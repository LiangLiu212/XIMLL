#ifndef AMPLITUDE_CUH
#define AMPLITUDE_CUH
#include "stdio.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <time.h>

struct AA_parameter{
	double alpha_jpsi;
	double phi_jpsi;
	double alpha_xi;
	double phi_xi;
	double alpha_xibar;
	double phi_xibar;
	double alpha1_lambda;
	double alpha1_lambdabar;
	double alpha2_lambda;
	double alpha2_lambdabar;
};

__device__ double device_AADecay12(const double alpha, const double xiphi,const double theta, const double phi, const int index);

__device__ double device_AAProd1212(double alpha_Jpsi, double delta_phi, double Jpsi_the, const int index);

__device__ double Amp(double *Matrix, const int Index);

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
				const int g_flag, double *Matrix);

#endif

