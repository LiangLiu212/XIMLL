#include "Amplitude.cuh"

__device__ double device_AADecay12(const double alpha, const double xiphi,const double theta, const double phi, const int index){

		double gamma = sqrt(1-alpha*alpha)*cos(xiphi);
		double beta = sqrt(1-alpha*alpha)*sin(xiphi);
		double c=0;
		switch(index){
				case 0: c=1; break;
				case 3: c=alpha; break;
				case 4: c=alpha*cos(phi)*sin(theta); break;
				case 5: c=gamma*cos(phi)*cos(theta) - beta*sin(phi); break;
				case 6: c=-(beta*cos(phi)*cos(theta)) - gamma*sin(phi); break;
				case 7: c=cos(phi)*sin(theta); break;
				case 8: c=alpha*sin(phi)*sin(theta); break;
				case 9: c=beta*cos(phi) + gamma*cos(theta)*sin(phi); break;
				case 10: c=gamma*cos(phi) - beta*cos(theta)*sin(phi); break;
				case 11: c=sin(phi)*sin(theta); break;
				case 12: c=alpha*cos(theta); break;
				case 13: c=-(gamma*sin(theta)); break;
				case 14: c=beta*sin(theta); break;
				case 15: c=cos(theta); break;
		}

		//printf("device_AADecay12 %d ~ %f\n", index, c);
		return c;
}


__device__ double  device_AAProd1212(double alpha_Jpsi, double delta_phi, double Jpsi_the, const int index){
		double v  = sqrt(1-alpha_Jpsi*alpha_Jpsi);
		double c=0;
		switch (index){
				case 0: c=(1+alpha_Jpsi*cos(Jpsi_the)*cos(Jpsi_the)); break;
				case 2: c=v*sin(Jpsi_the)*cos(Jpsi_the)*sin(delta_phi); break;
				case 5: c=sin(Jpsi_the)*sin(Jpsi_the); break;
				case 7: c=v*sin(Jpsi_the)*cos(Jpsi_the)* cos(delta_phi); break;
				case 8: c=-v*sin(Jpsi_the)*cos(Jpsi_the)*sin(delta_phi); break;
				case 10: c=alpha_Jpsi*sin(Jpsi_the)*sin(Jpsi_the); break;
				case 13: c=-v*sin(Jpsi_the)*cos(Jpsi_the)* cos(delta_phi); break;
				case 15: c=-alpha_Jpsi - cos(Jpsi_the)*cos(Jpsi_the); break;
		}
		//printf("device_AADecay1212 %d ~ %f\n", index, c);
		return c;
}


__device__ double Amp(double *Matrix, const int Index){
		double amp = 0;
		int i1 = 0;
		int i2 = 16;
		int i3 = 32;
		int i4 = 48;
		int i5 = 64;
		for(int mu=0; mu<4;mu++){// Xi loop
				for(int nu=0;nu<4;nu++){// Xibar loop
						for(int k=0;k<4;k++){
								for(int j=0;j<4;j++){
										amp += Matrix[Index*80 + i1 + mu*4 + nu]*  // G_numu
												Matrix[Index*80 + i2 + mu*4 + k]*  // A_numu for Xi decay
												Matrix[Index*80 + i3 + nu*4 + j]*  // A_numu for Xibar decay
												Matrix[Index*80 + i4 + k*4 + 0]*	// A_numu for lambda decay
												Matrix[Index*80 + i5 + j*4 + 0];  	// A_numu for lambdabar decay
								}
						}
				}
		}
		return amp;
}

__global__ void gpu_Amp(
				double *g_xithe,
				double *g_lthe,
				double *g_lphi,
				double *g_lbthe,
				double *g_lbphi,
				double *g_pthe,
				double *g_pphi,
				double *g_apthe,
				double *g_apphi,
				double *amp,
				const int g_NN,
				const AA_parameter g_para,
				const int g_flag, double *Matrix, const double norm){

		int index = blockIdx.x * blockDim.x + threadIdx.x;
		if( index <  g_NN){
				int iEvt = index / 80;

				int IMat = index % 80;
				int iMat = IMat / 16;
				int jMat = IMat % 16;
				if(g_flag  == 0){

						switch (iMat) {
								case 0: Matrix[index] = device_AAProd1212(g_para.alpha_jpsi, g_para.phi_jpsi, g_xithe[iEvt], jMat); break;
								case 1: Matrix[index] = device_AADecay12(g_para.alpha_xi, g_para.phi_xi, g_lthe[iEvt], g_lphi[iEvt], jMat); break;
								case 2: Matrix[index] = device_AADecay12(g_para.alpha_xibar, g_para.phi_xibar, g_lbthe[iEvt], g_lbphi[iEvt], jMat); break;
								case 3: Matrix[index] = device_AADecay12(g_para.alpha1_lambda, 0, g_pthe[iEvt], g_pphi[iEvt], jMat); break;
								case 4: Matrix[index] = device_AADecay12(g_para.alpha1_lambdabar, 0, g_apthe[iEvt], g_apphi[iEvt], jMat); break;
						}

				}
				else if(g_flag  == 1){

						switch (iMat) {
								case 0: Matrix[index] = device_AAProd1212(g_para.alpha_jpsi, g_para.phi_jpsi, g_xithe[iEvt], jMat); break;
								case 1: Matrix[index] = device_AADecay12(g_para.alpha_xi, g_para.phi_xi, g_lthe[iEvt], g_lphi[iEvt], jMat); break;
								case 2: Matrix[index] = device_AADecay12(g_para.alpha_xibar, g_para.phi_xibar, g_lbthe[iEvt], g_lbphi[iEvt], jMat); break;
								case 3: Matrix[index] = device_AADecay12(g_para.alpha2_lambda, 0, g_pthe[iEvt], g_pphi[iEvt], jMat); break;
								case 4: Matrix[index] = device_AADecay12(g_para.alpha2_lambdabar, 0, g_apthe[iEvt], g_apphi[iEvt], jMat); break;
						}
				}
				__syncthreads();

				if(IMat == 0){
						amp[iEvt] = log( Amp(Matrix, iEvt) / norm);
						//printf("%d ~ %f\n", iEvt, amp[iEvt]);
				}
		}
		__syncthreads();
}


// Bandwidth: (((2^27) + 1) unsigned ints * 4 bytes/unsigned int)/(38.716 * 10^-3 s)
//  13.867 GB/s = 96.297% -> excellent memory bandwidth
// Reasonable point to stop working on this implementation's optimization
// Algorithm is not compute-intensive, so acheiving >75% of theoretical bandwidth is goal
// Main strategies used:
// - Process as much data as possible (in terms of algorithm correctness) in shared memory
// - Use sequential addressing to get rid of bank conflicts
__global__
void block_sum_reduce(double* const d_block_sums, 
	double* const d_in,
	const unsigned int d_in_len)
{
	extern __shared__ unsigned int s_out[];

	unsigned int max_elems_per_block = blockDim.x * 2;
	unsigned int glbl_tid = blockDim.x * blockIdx.x + threadIdx.x;
	unsigned int tid = threadIdx.x;
	
	// Zero out shared memory
	// Especially important when padding shmem for
	//  non-power of 2 sized input
	s_out[threadIdx.x] = 0;
	s_out[threadIdx.x + blockDim.x] = 0;

	__syncthreads();

	// Copy d_in to shared memory per block
	if (glbl_tid < d_in_len)
	{
		s_out[threadIdx.x] = d_in[glbl_tid];
		if (glbl_tid + blockDim.x < d_in_len)
			s_out[threadIdx.x + blockDim.x] = d_in[glbl_tid + blockDim.x];
	}
	__syncthreads();

	// Actually do the reduction
	for (unsigned int s = blockDim.x; s > 0; s >>= 1) {
		if (tid < s) {
			s_out[tid] += s_out[tid + s];
		}
		__syncthreads();
	}

	// write result for this block to global mem
	if (tid == 0)
		d_block_sums[blockIdx.x] = s_out[0];
}







