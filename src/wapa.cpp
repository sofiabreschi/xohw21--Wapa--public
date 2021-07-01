//#include "bwamem.h"
#define DIMREF 1024
#define UNROLL_FACTOR 32
#define MINUS_INF -32000


#include <ap_int.h>
#include <stdio.h>


void computeAntidiag(unsigned short *antiDiag1_M,
							unsigned short *antiDiag2_M,
							unsigned short *antiDiag3_M,
							unsigned short *antiDiag2_I,
							unsigned short *antiDiag3_I,
							unsigned short *antiDiag2_D,
							unsigned short *antiDiag3_D,
							char *query,
							char *target,
							int qlen,
							int tlen,
							int minCol,
							int maxCol,
							int minRow,
							int maxRow,
							int antiDiagNo,
							int8_t open_gap,
							int8_t extend_gap,
							int8_t match,
							int8_t mismatch,
							int st,
							int en,
							int st_old,
							int en_old)
	{



		int8_t open_extend = open_gap + extend_gap;
		
		// printf("minRow %d\n", minRow);
		to_unroll:for(int i = 0; i < DIMREF; i+=UNROLL_FACTOR){
#pragma HLS DEPENDENCE variable=antiDiag1_M false inter
//#pragma HLS DEPENDENCE variable=antiDiag1_M false inter
#pragma HLS DEPENDENCE variable=antiDiag2_D false inter
			for(int j = 0; j < UNROLL_FACTOR; j++){
				if(st <= i + j && i + j < en){

					int target_position = minRow + i + j - minCol;
					int query_position = i + j ;
					int penalty_position = tlen - UNROLL_FACTOR -1;

					int tmp_h_M = 0, tmp_D = 0, tmp_v_M = 0, tmp_I = 0;


					int32_t M, D, I, tmpM;

					if(i + j  == 0){
						I = antiDiagNo*open_gap;
					}else if(i + j ==st && en > en_old && st_old == st){
						I = MINUS_INF;
					}else{
						I = antiDiag2_I[i + j -1];
					}
					if(i == maxCol-1 && antiDiagNo < qlen){
						M = (i + j )==0 ? 0 : antiDiagNo*open_gap;
						D = antiDiagNo*open_gap;
					}else if(i + j  == en-1 && en > en_old){
						M = antiDiag1_M[i + j -1];
						D = MINUS_INF;
					}else{
						M = (i + j )==0 ? antiDiagNo*open_gap : antiDiag1_M[i + j -1];
						D = antiDiag2_D[i + j];
					}

					//associate match score

					M += query[query_position] == target[target_position] ? match : mismatch;

					//check best dependency

					M = M > D ? M : D;
					M = M > I ? M : I;

					//aggiorno h successivo
					// printf("d: %d ", d);
					tmpM = M;
					//update e
					M -= open_extend;
					D -= extend_gap;
					D = M > D ? M : D;

					//update f
					I -= extend_gap;
					I = M > I ? M : I;

					//store new values
					antiDiag3_M[i + j] = tmpM;
					antiDiag3_D[i + j] = M > D ? M : D;
					antiDiag3_I[i + j] = M > I ? M : I;
				}
			}
		}
	}

extern "C"{
void SW (ap_uint<512> *query_global,
		ap_uint<512> *target_global,
		const int qlen,
		const int tlen, 
		const int n_couples,
		int bandwidth, 
		int *result){
   

#pragma HLS INTERFACE m_axi port = query_global offset = slave bundle = gmem0
#pragma HLS INTERFACE m_axi port = target_global offset = slave bundle = gmem1
#pragma HLS INTERFACE m_axi port = result offset = slave bundle = gmem2

#pragma HLS INTERFACE s_axilite port = query_global bundle = control
#pragma HLS INTERFACE s_axilite port = target_global bundle = control
#pragma HLS INTERFACE s_axilite port = result bundle = control

#pragma HLS INTERFACE s_axilite port = qlen bundle = control
#pragma HLS INTERFACE s_axilite port = tlen bundle = control
#pragma HLS INTERFACE s_axilite port = n_couples bundle = control
#pragma HLS INTERFACE s_axilite port = bandwidth bundle = control

#pragma HLS INTERFACE s_axilite port = return bundle = control

#pragma HLS INTERFACE ap_ctrl_chain port=return bundle=control

	char *query, *target;
	ap_uint<512> query_local[DIMREF/64], target_local[DIMREF/64];

	for(int i = 0; i < DIMREF*n_couples/64; i++){
		query_local[i] = query_global[i];
		target_local[i] = target_global[i];
	}

	query = (char*)query_local;
	target = (char*)target_local;

	unsigned short antiDiag1_M_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag1_M_buffer dim=1 factor=32 cyclic
	unsigned short antiDiag2_M_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag2_M_buffer dim=1 factor=32 cyclic
	unsigned short antiDiag3_M_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag3_M_buffer dim=1 factor=32 cyclic

	unsigned short antiDiag2_I_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag2_I_buffer dim=1 factor=32 cyclic
	unsigned short antiDiag3_I_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag3_I_buffer dim=1 factor=32 cyclic

	unsigned short antiDiag2_D_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag2_D_buffer dim=1 factor=32 cyclic
	unsigned short antiDiag3_D_buffer[DIMREF];
#pragma HLS ARRAY_PARTITION variable=antiDiag3_D_buffer dim=1 factor=32 cyclic


	unsigned short *antiDiag1_M = antiDiag1_M_buffer;
	unsigned short *antiDiag2_M = antiDiag2_M_buffer
	unsigned short *antiDiag3_M = antiDiag3_M_buffer;
	//
	unsigned short *antiDiag2_I = antiDiag2_I_buffer;
	unsigned short *antiDiag3_I = antiDiag3_I_buffer;
	//
	unsigned short *antiDiag2_D = antiDiag2_D_buffer;
	unsigned short *antiDiag3_D = antiDiag3_D_buffer;

	int8_t match = 1;
	int8_t mismatch = 17;
	int8_t extend_gap = 1;
	int8_t open_gap = 6;
	
	for(int n = 0; n < n_couples; n++){
	
		int minCol = 0;
		int maxCol = 0;
		int minRow = tlen;
		int maxRow = 0;
		
		int en_old = 0, st_old = 0;
		
		for (int antiDiagNo = 0; antiDiagNo < tlen + qlen - 1; ++antiDiagNo) { // antidiagonals in the outer loop

			unsigned short *t_M = antiDiag1_M;
			antiDiag1_M = antiDiag2_M;
			antiDiag2_M = antiDiag3_M;
			antiDiag3_M = t_M;
			//
			unsigned short *t_I = antiDiag2_I;
			antiDiag2_I = antiDiag3_I;
			antiDiag3_I = t_I;
			//
			unsigned short *t_D = antiDiag2_D;
			antiDiag2_D = antiDiag3_D;
			antiDiag3_D = t_D;
			//
			unsigned short bandwidthl, bandwidthr;
			if (bandwidth < 0) bandwidth = tlen > qlen? tlen : qlen;
			bandwidthl = bandwidthr = bandwidth;
			//
			int st = 0, en = tlen - 1;
			// find the boundaries
			if (st < antiDiagNo - qlen + 1) st = antiDiagNo - qlen + 1;
			if (en > antiDiagNo) en = antiDiagNo;
			if (st < (antiDiagNo-bandwidthr+1)>>1) st = (antiDiagNo-bandwidthr+1)>>1;
			if (en > (antiDiagNo+bandwidthl)>>1) en = (antiDiagNo+bandwidthl)>>1;

			if(antiDiagNo < qlen)
				maxCol++;
			if (antiDiagNo < tlen){
				minRow--;
				maxRow++;
			}
			if(antiDiagNo >= tlen)
				minCol++;

			if(bandwidth >= qlen){
					st = minCol;
					en = maxCol-1;
			}

			computeAntidiag(antiDiag1_M,
							antiDiag2_M,
							antiDiag3_M,
							antiDiag2_I,
							antiDiag3_I,
							antiDiag2_D,
							antiDiag3_D,
							query,target,
							qlen,tlen,
							minCol,maxCol,
							minRow,maxRow,
							antiDiagNo,
							open_gap,extend_gap,
							match,mismatch,
							st,en,st_old,en_old);


			en_old = en;
			st_old = st;

			if(antiDiagNo == qlen + tlen - 2){
				result[n] = antiDiag3_M[minCol];
			}
		}
	}

}
}