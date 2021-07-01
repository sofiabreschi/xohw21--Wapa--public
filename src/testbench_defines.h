#ifndef bwamem_h
#define bwamem_h

#include <stdio.h>
#include <assert.h>
#include <limits.h>
#include <math.h>
/*#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

#include "kstring.h"
#include "bwamem.h"
#include "bntseq.h"
#include "ksw.h"
#include "kvec.h"
#include "ksort.h"
#include "utils.h"

#ifdef USE_MALLOC_WRAPPERS
#  include "malloc_wrap.h"
#endif
*/
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <stdint.h>
#include <fstream>
#include "ap_int.h"
/*#include "bwamem.h"
#include "bwt.h"
#include "bntseq.h"
#include "bwa.h"
*/
using namespace std;

void SW (ap_uint<512> a[], ap_uint<512> b[], const int qlen, const int tlen, const int bandwidth, int score, int *result);
//void reverseStr(char[]);
//char random_char();
//int random_int();
//void SW (char[] , char[] ,const int, int );
//char* init_array(char*,const int );
//int** init_mat(int**,const int);
//char** init_mat_char(char**, int);
//void print1(const char[], const int);
//void print2(const char[], const int);
//void free_mat(int**,const int);
//void free_mat_char(char**, const int);
//void free_array(char*);
//void reverseStr(char[], const int);
#endif /* bwamem_h */
