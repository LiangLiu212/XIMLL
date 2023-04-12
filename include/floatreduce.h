#ifndef FLOAT_REDUCE_H__
#define FLOAT_REDUCE_H__
#include "utils.h"
#define MAX_BLOCK_SZ 1024
#define checkCudaErrors(val) check( (val), #val, __FILE__, __LINE__)

double gpu_float_sum_reduce(double* d_in, unsigned int d_in_len);

#endif // !REDUCE_H__


