#ifndef GPU_ANGDISXIXI_H
#define GPU_ANGDISXIXI_H
#include <cmath>
#include <iostream>
#include "device_launch_parameters.h"
#include "cuda_runtime.h"


struct AA_Matrix{
	double thC[4][4];
	double tHa[4][4];
	double tHb[4][4];
	double tHc[4][4];
	double tHd[4][4];
};

__device__ void  gpu_AAProd1212(double alpha_Jpsi, double delta_phi, double Jpsi_the, double thC[4][4]){
		Double_t v  = sqrt(1-alpha_Jpsi*alpha_Jpsi);
	//	Double_t st = __sinf(Jpsi_the);
	//	Double_t ct = __cosf(Jpsi_the);
	//	Double_t cp = __cosf(delta_phi);
	//	Double_t sp = __sinf(delta_phi);
		for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
						Double_t c=0;
						if(k==0&&l==0)c=(1+alpha_Jpsi*__cosf(Jpsi_the)*__cosf(Jpsi_the));
						if(k==0&&l==2)c=v*__sinf(Jpsi_the)*__cosf(Jpsi_the)*__sinf(delta_phi);
						if(k==1&&l==1)c=__sinf(Jpsi_the)*__sinf(Jpsi_the);
						if(k==1&&l==3)c=v*__sinf(Jpsi_the)*__cosf(Jpsi_the)* __cosf(delta_phi);
						if(k==2&&l==0)c=-v*__sinf(Jpsi_the)*__cosf(Jpsi_the)*__sinf(delta_phi);
						if(k==2&&l==2)c=alpha_Jpsi*__sinf(Jpsi_the)*__sinf(Jpsi_the);
						if(k==3&&l==1)c=-v*__sinf(Jpsi_the)*__cosf(Jpsi_the)* __cosf(delta_phi);
						if(k==3&&l==3)c=-alpha_Jpsi - __cosf(Jpsi_the)*__cosf(Jpsi_the);
						thC[k][l]=c;
				}
		}
}


__device__ void gpu_AADecay12(const double alpha, const double xiphi, const double theta, const double phi, double tHa[4][4]) {

	//	double cos_theta = __cosf(theta);
	//	double sin_theta = __sinf(theta);
	//	double cos_phi = __cosf(phi);
	//	double sin_phi = __sinf(phi);
		// decay 1/2 -> 1/2 + 0
		double gamma = sqrt(1-alpha*alpha)*__cosf(xiphi);
		double beta = sqrt(1-alpha*alpha)*__sinf(xiphi);

		for(int k=0;k<4;k++){
				for(int l=0;l<4;l++){
						Double_t c=0;

						if(k==0&& l == 0) c=1;

						if(k==0&& l == 3) c=alpha;

						if(k==1&& l == 0) c=alpha*__cosf(phi)*__sinf(theta);

						if(k==1&& l == 1) c=gamma*__cosf(phi)*__cosf(theta) - beta*__sinf(phi);

						if(k==1&& l == 2) c=-(beta*__cosf(phi)*__cosf(theta)) - gamma*__sinf(phi);

						if(k==1&& l == 3) c=__cosf(phi)*__sinf(theta);

						if(k==2&& l == 0) c=alpha*__sinf(phi)*__sinf(theta);

						if(k==2&& l == 1) c=beta*__cosf(phi) + gamma*__cosf(theta)*__sinf(phi);

						if(k==2&& l == 2) c=gamma*__cosf(phi) - beta*__cosf(theta)*__sinf(phi);

						if(k==2&& l == 3) c=__sinf(phi)*__sinf(theta);

						if(k==3&& l == 0) c=alpha*__cosf(theta);

						if(k==3&& l == 1) c=-(gamma*__sinf(theta));

						if(k==3&& l == 2) c=beta*__sinf(theta);

						if(k==3&& l == 3) c=__cosf(theta);

						tHa[k][l]=c;
				}
		}


}



__device__ void Amp(double theta, double LThe,double LPhi,double LbThe, double LbPhi,double pThe,double pPhi, double apThe, double apPhi, double *pp, AA_Matrix *g_munu){

/*
		for(int i=0;i<4;i++){
				for(int j=0;j<4;j++){
						g_munu->thC[i][j] = 0;
						g_munu->tHa[i][j] = 0;
						g_munu->tHb[i][j] = 0;
						g_munu->tHc[i][j] = 0;
						g_munu->tHd[i][j] = 0;
				}
		}
*/
		//	printf("C i j : %f\n", Prob[0][0] );
	//	double  p_Jpsi_alpha = pp[0];
	//	double  p_Jpsi_phi = pp[1];
	//	double  p_Xi_alpha = pp[2];
	//	double  p_Xi_phi = pp[3];
	//	double  p_Xibar_alpha = pp[4];
	//	double  p_Xibar_phi = pp[5];
	//	double  p_L_alpha = pp[6];
	//	double  p_Lbar_alpha = pp[7];
		gpu_AAProd1212( pp[0], pp[1], theta, g_munu->thC);
		gpu_AADecay12(pp[2], pp[3], LThe, LPhi, g_munu->tHa);
		gpu_AADecay12(pp[4], pp[5], LbThe, LbPhi, g_munu->tHb);
		gpu_AADecay12(pp[6], 0, pThe, pPhi, g_munu->tHc);
		gpu_AADecay12(pp[7], 0, apThe, apPhi, g_munu->tHd);
	/*	double tep = 0;

		for(int mu=0; mu<4;mu++){// Xi loop
				for(int nu=0;nu<4;nu++){// Xibar loop
						for(int k=0;k<4;k++){
								for(int j=0;j<4;j++){
										tep +=  g_munu->thC[mu][nu]*
												g_munu->tHa[mu][k]*g_munu->tHb[nu][j]*
												g_munu->tHc[k][0]*g_munu->tHd[j][0];
								}
						}
				}
		}
		*/
	//	printf("tep : %f\n", p_Jpsi_alpha);
}

#endif // GPU_ANGDISXIXI_H
