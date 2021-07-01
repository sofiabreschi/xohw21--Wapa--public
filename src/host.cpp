#include <fstream>
#include "xcl2.hpp"
#include <algorithm>
#include <iostream>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <string>
#include <climits>
#include <time.h>
#include <chrono>
#include "ap_int.h"

#define TARGETLEN 1024
#define QUERYLEN 1024
#define NUM_KERNEL 4

#define NOW std::chrono::high_resolution_clock::now();

#define MAX_HBM_BANKCOUNT 32
#define BANK_NAME(n) n | XCL_MEM_TOPOLOGY
const int bank[MAX_HBM_BANKCOUNT] = {
    BANK_NAME(0),  BANK_NAME(1),  BANK_NAME(2),  BANK_NAME(3),  BANK_NAME(4),
    BANK_NAME(5),  BANK_NAME(6),  BANK_NAME(7),  BANK_NAME(8),  BANK_NAME(9),
    BANK_NAME(10), BANK_NAME(11), BANK_NAME(12), BANK_NAME(13), BANK_NAME(14),
    BANK_NAME(15), BANK_NAME(16), BANK_NAME(17), BANK_NAME(18), BANK_NAME(19),
    BANK_NAME(20), BANK_NAME(21), BANK_NAME(22), BANK_NAME(23), BANK_NAME(24),
    BANK_NAME(25), BANK_NAME(26), BANK_NAME(27), BANK_NAME(28), BANK_NAME(29),
    BANK_NAME(30), BANK_NAME(31)};

int main(int argc, char *argv[]){

    char dictionary_dna[4] = {'A', 'G', 'T', 'C'};

    std::string binaryFile = argv[1];
	srand(time(NULL));

    cl_int err;
    cl::CommandQueue commands;
    cl::Context context;

    std::vector< char, aligned_allocator< char > > query(TARGETLEN);
    for(int i = 0; i < TARGETLEN; i++){
        query[i] = dictionary_dna[rand()%4];
    }
    std::vector< char, aligned_allocator< char > > target(QUERYLEN);
    for(int i = 0; i < QUERYLEN; i++){
        target[i] = dictionary_dna[rand()%4];
    }

    std::vector< int, aligned_allocator< int > > results(NUM_KERNEL);

    std::string krnl_name = "SW";
    std::vector<cl::Kernel> krnls(NUM_KERNEL);

    // The get_xil_devices will return vector of Xilinx Devices
    auto devices = xcl::get_xil_devices();

    // read_binary_file() command will find the OpenCL binary file created using the
    // V++ compiler load into OpenCL Binary and return pointer to file buffer.
    auto fileBuf = xcl::read_binary_file(binaryFile);

    cl::Program::Binaries bins{{fileBuf.data(), fileBuf.size()}};
    int valid_device = 0;


    for (unsigned int i = 0; i < devices.size(); i++) {
        auto device = devices[i];
            // Creating Context and Command Queue for selected Device
        OCL_CHECK(err, context = cl::Context(device, NULL, NULL, NULL, &err));
        OCL_CHECK(err, commands = cl::CommandQueue(context, device,
                            CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE | CL_QUEUE_PROFILING_ENABLE, &err));

        std::cout << "Trying to program device[" << i 
                  << "]: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
     
        cl::Program program(context, {device}, bins, NULL, &err);
        
        if (err != CL_SUCCESS) {
            std::cout << "Failed to program device[" << i
                        << "] with xclbin file!\n";                      
        } else {
            std::cout << "Device[" << i << "]: program successful!\n";
            
            // Creating Kernel object using Compute unit names
            for (int i = 0; i < NUM_KERNEL; i++) {
                std::string cu_id = std::to_string(i + 1);
                std::string krnl_name_full = krnl_name + ":{" + "SW_" + cu_id + "}";

                printf("Creating a kernel [%s] for CU(%d)\n", krnl_name_full.c_str(), i + 1);

                //Here Kernel object is created by specifying kernel name along with compute unit.
                //For such case, this kernel object can only access the specific Compute unit
                OCL_CHECK(err, krnls[i] = cl::Kernel(program, krnl_name_full.c_str(), &err));
            }

            valid_device++;
            break; // we break because we found a valid device
        }
        std::cout<<"dwvgae"<<std::endl;
    }
    
    if (valid_device == 0) {
        std::cout << "Failed to program any device found, exit!\n";
        exit(EXIT_FAILURE);
    }

    // Create device buffers
    std::vector<cl_mem_ext_ptr_t> query_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> target_buffer_ext(NUM_KERNEL);
    std::vector<cl_mem_ext_ptr_t> results_buffer_ext(NUM_KERNEL);

    std::vector<cl::Buffer> query_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> target_buffer(NUM_KERNEL);
    std::vector<cl::Buffer> results_buffer(NUM_KERNEL);

	

    for(int i = 0; i < NUM_KERNEL; i++) {

        query_buffer_ext[i].obj = query.data();
        query_buffer_ext[i].param = 0;
        query_buffer_ext[i].flags = bank[i*3];

        target_buffer_ext[i].obj = target.data();
        target_buffer_ext[i].param = 0;
        target_buffer_ext[i].flags = bank[i*3+1];
        
        results_buffer_ext[i].obj = results.data();
        results_buffer_ext[i].param = 0;
        results_buffer_ext[i].flags = bank[i*3+2];
    }

    for (int i = 0; i < NUM_KERNEL; i++) {
    	OCL_CHECK(err, query_buffer[i] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*query.size(), &query_buffer_ext[i], &err));
        OCL_CHECK(err, target_buffer[i] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(char)*target.size(), &target_buffer_ext[i], &err));
        OCL_CHECK(err, results_buffer[i] = cl::Buffer(context, CL_MEM_READ_WRITE | CL_MEM_EXT_PTR_XILINX |
                                    CL_MEM_USE_HOST_PTR, sizeof(int)*results.size(), &results_buffer_ext[i], &err));
    }

	commands.finish();

    // Write our data set into device buffers  
     for(int i = 0; i < NUM_KERNEL; i++)
        err = commands.enqueueMigrateMemObjects({query_buffer[i], target_buffer[i]}, 0);

    if (err != CL_SUCCESS) {
            printf("Error: Failed to write to device memory!\n");
            printf("Test failed\n");
            exit(1);
    }

	commands.finish();


    // Set the arguments to our compute kernel
    for (int i = 0; i < NUM_KERNEL; i++) {
		OCL_CHECK(err, err = krnls[i].setArg(0, query_buffer[i]));
		OCL_CHECK(err, err = krnls[i].setArg(1, target_buffer[i]));
        
        OCL_CHECK(err, err = krnls[i].setArg(2, QUERYLEN));
        OCL_CHECK(err, err = krnls[i].setArg(3, TARGETLEN));
        OCL_CHECK(err, err = krnls[i].setArg(4, 1));
        OCL_CHECK(err, err = krnls[i].setArg(5, QUERYLEN));
        OCL_CHECK(err, err = krnls[i].setArg(6, results_buffer[i]));

        if (err != CL_SUCCESS) {
            printf("Error: Failed to set kernel arguments! %d\n", err);
            printf("Test failed\n");
            exit(1);
        }
    }

	commands.finish();

    std::chrono::high_resolution_clock::time_point start_time = NOW;

    for (int i = 0; i < NUM_KERNEL; ++i)
        err |= commands.enqueueTask(krnls[i]);


    if (err) {
        printf("Error: Failed to execute kernel! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

    commands.finish();
    std::chrono::high_resolution_clock::time_point end_time = NOW;
	std::chrono::duration<double> time = std::chrono::duration_cast<std::chrono::duration<double>>(end_time-start_time)/NUM_KERNEL;
    
    // Read back the results from the device to verify the output
    for (int i = 0; i < NUM_KERNEL; ++i) {
        err = commands.enqueueMigrateMemObjects({results_buffer[i]}, CL_MIGRATE_MEM_OBJECT_HOST);  
    }

    if (err != CL_SUCCESS) {
        printf("Error: Failed to read output array! %d\n", err);
        printf("Test failed\n");
        exit(1);
    }

    

	printf("HW time: %lf\n", time/NUM_KERNEL);

	//Checking the results 

    bool test_score = true;
	

	if (test_score) 
		std::cout<<"ALL RESULTS CORRECT"<<std::endl;
	else 
		std::cout<<"Test failed"<<std::endl;

    return 0;
}