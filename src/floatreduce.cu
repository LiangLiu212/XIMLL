#include "floatreduce.h"

#include <cmath>

#include "cuda_runtime.h"
#include "device_launch_parameters.h"

// Bandwidth: (((2^27) + 1) unsigned ints * 4 bytes/unsigned int)/(38.716 * 10^-3 s)
//  13.867 GB/s = 96.297% -> excellent memory bandwidth
// Reasonable point to stop working on this implementation's optimization
// Algorithm is not compute-intensive, so acheiving >75% of theoretical bandwidth is goal
// Main strategies used:
// - Process as much data as possible (in terms of algorithm correctness) in shared memory
// - Use sequential addressing to get rid of bank conflicts
__global__
void block_float_sum_reduce(double*  d_block_sums, 
	const double*  d_in,
	const unsigned int d_in_len)
{
	extern __shared__ double s_out[];

	unsigned int max_elems_per_block = blockDim.x;
	unsigned int glbl_tid = (blockDim.x) * blockIdx.x + threadIdx.x;
	unsigned int tid = threadIdx.x;
	
	// Zero out shared memory
	// Especially important when padding shmem for
	//  non-power of 2 sized input
//	s_out[threadIdx.x] = 0;
//	s_out[threadIdx.x + blockDim.x] = 0;

//	if (glbl_tid < d_in_len){
			s_out[tid] = d_in[glbl_tid];
//	}

	__syncthreads();

	// Actually do the reduction
	for (unsigned int s = 1; s < blockDim.x; s *= 2) {
			unsigned int index  = 2 *s *tid;
			if ( index < blockDim.x && index + s <  blockDim.x && index + s < d_in_len) {
					s_out[index] += s_out[index + s ];
			}
			__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0)
			d_block_sums[blockIdx.x] = s_out[0];
	__syncthreads();
}


double gpu_float_sum_reduce(double* d_in, unsigned int d_in_len)
{
		double total_sum = 0;

		// Set up number of threads and blocks
		// If input size is not power of two, the remainder will still need a whole block
		// Thus, number of blocks must be the least number of 2048-blocks greater than the input size
		unsigned int block_sz = MAX_BLOCK_SZ; // Halve the block size due to reduce3() and further 
		//  optimizations from there
		// our block_sum_reduce()
		unsigned int max_elems_per_block = block_sz ; // due to binary tree nature of algorithm
		// NVIDIA's reduceX()
		//unsigned int max_elems_per_block = block_sz;

		unsigned int grid_sz = 0;
		if (d_in_len <= max_elems_per_block)
		{
				grid_sz = (unsigned int)std::ceil(float(d_in_len) / float(max_elems_per_block));
		}
		else
		{
				grid_sz = (d_in_len ) / max_elems_per_block;
				if (d_in_len % max_elems_per_block != 0)
						grid_sz++;
		}

		// Allocate memory for array of total sums produced by each block
		// Array length must be the same as number of blocks / grid size
		double* d_block_sums;
		checkCudaErrors(cudaMalloc(&d_block_sums, sizeof(double) * grid_sz));
		checkCudaErrors(cudaMemset(d_block_sums, 0, sizeof(double) * grid_sz));

		// Sum data allocated for each block
		block_float_sum_reduce<<<grid_sz, block_sz, sizeof(double) * max_elems_per_block>>>(d_block_sums, d_in, d_in_len);

	//	std::cout << grid_sz << "	" << block_sz << "	" << max_elems_per_block << "  " << d_in_len << std::endl;
		//	reduce4<<<grid_sz, block_sz, sizeof(double) * block_sz>>>(d_block_sums, d_in, d_in_len);
		//print_d_array(d_block_sums, grid_sz);

		// Sum each block's total sums (to get global total sum)
		// Use basic implementation if number of total sums is <= 2048
		// Else, recurse on this same function
		if (grid_sz <= max_elems_per_block)
		{
				double* d_total_sum;
				checkCudaErrors(cudaMalloc(&d_total_sum, sizeof(double)));
				checkCudaErrors(cudaMemset(d_total_sum, 0, sizeof(double)));
				block_float_sum_reduce<<<1, block_sz, sizeof(double) * max_elems_per_block>>>(d_total_sum, d_block_sums, grid_sz);
		//		std::cout << 1 << "	" << block_sz << "	" << max_elems_per_block << "  " << grid_sz << std::endl;
				//	reduce4<<<1, block_sz, sizeof(double) * block_sz>>>(d_total_sum, d_block_sums, grid_sz);
				checkCudaErrors(cudaMemcpy(&total_sum, d_total_sum, sizeof(double), cudaMemcpyDeviceToHost));
				checkCudaErrors(cudaFree(d_total_sum));
		}
		else
		{
				double* d_in_block_sums;
				checkCudaErrors(cudaMalloc(&d_in_block_sums, sizeof(double) * grid_sz));
				checkCudaErrors(cudaMemcpy(d_in_block_sums, d_block_sums, sizeof(double) * grid_sz, cudaMemcpyDeviceToDevice));
				total_sum = gpu_float_sum_reduce(d_in_block_sums, grid_sz);
				checkCudaErrors(cudaFree(d_in_block_sums));
		}

		checkCudaErrors(cudaFree(d_block_sums));
		return total_sum;
}
