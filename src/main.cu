#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "AngDisXiXi.hh"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVectorT.h"
#include "TStopwatch.h"
#include "TMath.h"
#include <TMinuit.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include "readData.h"
#include "gpu_AngDisXiXi.hh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#define _SLOOW
const Int_t NUM = 10000000;
#define THREADS_PER_BLOCK 1024
Int_t NN[4][6];
AngDisXiXi *angdis[4][2];
double **angdata[4][6];
double **gpu_angdata[4][6];
static int years;
int fit_flag = 0;
int fit_step = 0;

std::vector<int> i_year;


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



__global__ void gpu_aa(double *g_xithe,
					   double *g_lthe,
					   double *g_lphi,
					   double *g_lbthe,
					   double *g_lbphi,
					   double *g_pthe,
					   double *g_pphi,
					   double *g_apthe,
					   double *g_apphi,
					   const int g_NN,
					   const AA_parameter g_para,
					   AA_Matrix *g_munu,
					   const int g_flag, double *g_eval){
		int index = blockIdx.x * blockDim.x + threadIdx.x;
		if(index < g_NN){
				__shared__ double pp1[8], pp2[8];
				pp1[0] = g_para.alpha_jpsi;
				pp1[1] = g_para.phi_jpsi;
				pp1[2] = g_para.alpha_xi;
				pp1[3] = g_para.phi_xi;
				pp1[4] = g_para.alpha_xibar;
				pp1[5] = g_para.phi_xibar;
				pp1[6] = g_para.alpha1_lambda;
				pp1[7] = g_para.alpha1_lambdabar;
				pp2[0] = g_para.alpha_jpsi;
				pp2[1] = g_para.phi_jpsi;
				pp2[2] = g_para.alpha_xi;
				pp2[3] = g_para.phi_xi;
				pp2[4] = g_para.alpha_xibar;
				pp2[5] = g_para.phi_xibar;
				pp2[6] = g_para.alpha2_lambda;
				pp2[7] = g_para.alpha2_lambdabar;
				__syncthreads();
				if(g_flag < 2){
						Amp(g_xithe[index], g_lthe[index], g_lphi[index], g_lbthe[index], 
										g_lbphi[index], g_pthe[index], g_pphi[index], g_apthe[index], g_apphi[index], 
										pp1, &g_munu[index]);
				}
				else {
						Amp(g_xithe[index], g_lthe[index], g_lphi[index], g_lbthe[index], 
										g_lbphi[index], g_pthe[index], g_pphi[index], g_apthe[index], g_apphi[index], 
										pp2, &g_munu[index]);

				}
		//		double tep = 0;
				g_eval[index] = 0;
				for(int mu=0; mu<4;mu++){// Xi loop
						for(int nu=0;nu<4;nu++){// Xibar loop
								for(int k=0;k<4;k++){
										for(int j=0;j<4;j++){
												g_eval[index] +=  g_munu[index].thC[mu][nu]*
														g_munu[index].tHa[mu][k]*g_munu[index].tHb[nu][j]*
														g_munu[index].tHc[k][0]*g_munu[index].tHd[j][0];
										}
								}
						}
				}
				
	//			g_eval[index] = 1.0;


				//	printf("C numu  0 0 : %f\n", g_munu[index].thC[0][0]);
				//	printf("C numu  0 1 : %f\n", g_munu[index].thC[0][1]);
				//	printf("C numu  0 2 : %f\n", g_munu[index].thC[0][2]);
				//	printf("C numu  0 3 : %f\n", g_munu[index].thC[0][3]);


				//		#include "amplitude.cxx"

				//	printf("from gpu :  %f \n", g_eval[index]);
		}
}

//=====================================================================
void ReadData(int flag[4], const int index, const int MM)
{
		//	int years = 2;
		const char *type[6] = {"DATA", "BKG", "PHSP", "DATA", "BKG", "PHSP"};
		//	int flag[4] = {2009, 2012, 2018, 2019};
		const char *file[6] = {
				"/data/liul/workarea/XIXI/Rec3/mdiyRecpm2012.root", 
				"/data/liul/workarea/XIXI/Rec3/mdiyRecpm2012.root",
				"/data/liul/workarea/XIXI/Rec3/phspRecpm2012.root",
				"/data/liul/workarea/XIXI/Rec3/mdiyRecpp2012.root",
				"/data/liul/workarea/XIXI/Rec3/mdiyRecpp2012.root", 
				"/data/liul/workarea/XIXI/Rec3/phspRecpp2012.root" };

		for(int i = 0; i < years; i++){
				for(int j = 0; j < 6; j ++){
						angdata[i][j] = new double * [9];
						for(int l = 0; l < 9; l++){
								*(angdata[i][j] + l) =  new double [NUM];
						}
				}
		}

		for(int i = 0; i < years; i++)
				for(int j = 0; j < 2; j ++){
								if(angdis[i][j]){
						}
						else{
								angdis[i][j] = new AngDisXiXi();
						}
				}

		for(int i = 0; i < years; i++){ 		// read data
				for(int j = 0; j < 6; j ++){
						int l = j / 3;
						NN[i][j] =  readData(file[j], angdis[i][l], angdata[i][j], flag[i], type[j], index, MM);
				}
		}
		for(int i = 0; i < years; i++){
				for(int j = 0; j < 6; j ++){
						cout << "N[" << i << "][" << j << "] : " << NN[i][j] << endl;
				}
		}
		double **temp_angdata[years][6]; // define a temporary array 
		for(int i = 0; i < years; i++){
				for(int j = 0; j < 6; j ++){
						temp_angdata[i][j] = new double * [9];
						for(int l = 0; l < 9; l++){
								*(temp_angdata[i][j] + l) =  new double [NN[i][j]];
								for(int k = 0; k< NN[i][j]; k ++){
										*(*(temp_angdata[i][j] + l) + k) = *(*(angdata[i][j] + l) + k);
								}
						}
				}
		}



		for(int i = 0; i < years; i++){    // copy data from cpu to gpu
				for(int j = 0; j < 6; j ++){
						int size1 = NN[i][j] *sizeof(double);
						gpu_angdata[i][j] = new double * [9];
						for(int l = 0; l < 9; l++){
								cudaMalloc( (void **) &(*(gpu_angdata[i][j] + l)), size1 );
								cudaMemcpy( *(gpu_angdata[i][j] + l), *(temp_angdata[i][j] + l), size1, cudaMemcpyHostToDevice );
								//	for ( int k = 0; k < NN[i][j]; k++ ){
								//			cout << "read file : " << *(*((temp_angdata[i][j]) + 0) + k) << endl;
								//	}
						}
				}
		}
}


void fcnMLLG(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pp, Int_t iflag)
{ 
		double pp1[8], pp2[8];
		for(int i = 0; i < 6; i++){
				pp1[i] = pp[i];
				pp2[i] = pp[i];
		}
		pp1[6] = pp[6]; pp1[7] = pp[7];
		pp2[6] = pp[8]; pp2[7] = pp[9];

		AA_parameter aa_para;
		aa_para.alpha_jpsi = pp[0];
		aa_para.phi_jpsi = pp[1];
		aa_para.alpha_xi = pp[2];
		aa_para.phi_xi = pp[3];
		aa_para.alpha_xibar = pp[4];
		aa_para.phi_xibar = pp[5];
		aa_para.alpha1_lambda = pp[6];
		aa_para.alpha1_lambdabar = pp[7];
		aa_para.alpha2_lambda = pp[8];
		aa_para.alpha2_lambdabar = pp[9];

		double *host_eval;
		double *gpu_eval;
		AA_Matrix *gpu_munu;
		AA_Matrix *host_munu;

		cudaError_t cudaStatus;
		double loglike[4][4];
		int idx[4] = {0, 1, 3, 4};  // data1 bgk1 data2 bkg2
		for(int i = 0; i < years; i ++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
				for (int j = 0; j < 4; j++){

						host_munu = new AA_Matrix [NN[i][idx[j]]];
						host_eval = new double [NN[i][idx[j]]];
						int mat_size = NN[i][idx[j]]*sizeof(*gpu_munu);
						cudaMalloc( (void **) &gpu_munu,  mat_size);
						int size = NN[i][idx[j]] * sizeof(*gpu_eval);
						cudaMalloc( (void **) &gpu_eval, size);
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 003!" << endl;
								exit(1);
						}


						gpu_aa <<< (NN[i][idx[j]] + (THREADS_PER_BLOCK-1)) / THREADS_PER_BLOCK, THREADS_PER_BLOCK >>> ( 
										*(gpu_angdata[i][idx[j]] + 0), 
										*(gpu_angdata[i][idx[j]] + 1), 
										*(gpu_angdata[i][idx[j]] + 2), 
										*(gpu_angdata[i][idx[j]] + 3), 
										*(gpu_angdata[i][idx[j]] + 4), 
										*(gpu_angdata[i][idx[j]] + 5), 
										*(gpu_angdata[i][idx[j]] + 6), 
										*(gpu_angdata[i][idx[j]] + 7), 
										*(gpu_angdata[i][idx[j]] + 8), 
										NN[i][idx[j]], 
										aa_para, gpu_munu, j, gpu_eval);
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 004!" << endl;
								exit(1);
						}

						cudaDeviceSynchronize();
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 001!" << endl;
								exit(1);
						}
						cudaMemcpy( host_munu, gpu_munu, mat_size, cudaMemcpyDeviceToHost );
						cudaMemcpy( host_eval, gpu_eval, size, cudaMemcpyDeviceToHost );
						cudaFree( gpu_munu );
						cudaFree( gpu_eval );
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 002!" << endl;
								exit(1);
						}
						loglike[i][j] = 0;
						for(int evt = 0; evt < NN[i][idx[j]]; evt++){
								//	cout << "host C munu 0: " << 	*(host_eval + evt) << endl;
								if(*(host_eval + evt) <= 0){ f=0; cout << "data : " << *(host_eval + evt) << endl;  return; }
								loglike[i][j] += TMath::Log(*(host_eval + evt));
								//	cout << "host C munu 1: " << 	(host_munu + evt)->tHa[0][1] << endl;
								//	cout << "host C munu 2: " << 	(host_munu + evt)->tHa[0][2] << endl;
								//	cout << "host C munu 3: " << 	(host_munu + evt)->tHa[0][3] << endl;
						}
						delete [] host_munu;
						host_munu = NULL;
						delete [] host_eval;
						host_eval = NULL;
				}
		}
	//	exit(1);

		double norm[4][2];
		for (int i = 0; i < years; i++){
				for (int j = 0; j < 2; j++){
						norm[i][j] = 0;
						norm[i][j] = angdis[i][j]->CalcToIntegral();
						norm[i][j]/=Double_t(NN[i][3*j + 2]);
				}
		}

		double N_BKG[4][2] = {{477.15, 448.882}, {2293.86, 2427.39}, {8827.38, 9803.34}, {8234.57, 9115.42}};
		int idx_year = 0;

		double llk = 0;
		double l1 = 0;
		double l2 = 0;

		for (int i = 0; i < years; i++){
				if(i_year[i] == 2009) idx_year = 0;
				if(i_year[i] == 2012) idx_year = 1;
				if(i_year[i] == 2018) idx_year = 2;
				if(i_year[i] == 2019) idx_year = 3;
			//	cout <<  i_year[i] <<  "  Background  : "  << N_BKG[idx_year][0] << " " << N_BKG[idx_year][1] << endl;
				if(fit_flag == 1){
						l1 = - loglike[i][0] + N_BKG[idx_year][0]*loglike[i][1]/Double_t(NN[i][1]) + (Double_t(NN[i][0]) - N_BKG[idx_year][0])*TMath::Log(norm[i][0]);
						l2 = - loglike[i][2] + N_BKG[idx_year][1]*loglike[i][3]/Double_t(NN[i][4]) + (Double_t(NN[i][3]) - N_BKG[idx_year][1])*TMath::Log(norm[i][1]);
				}
				else if(fit_flag == 2){
						l1 = - loglike[i][0] + (Double_t(NN[i][0]))*TMath::Log(norm[i][0]);
						l2 = - loglike[i][2] + (Double_t(NN[i][3]))*TMath::Log(norm[i][1]);
				}
				else{
						cerr << "error fit_flag!" << endl;
				}
				llk += (l1 + l2);
		}

		f = llk;
		if(fit_step%100 == 0){
		std::cout << "Loglike: " << f << std::endl; 
		for( int i = 0; i<10 ; i++ ) cout<<pp[i]<<" ";
		cout << endl;
		}
		fit_step++;
}
//=====================================================================
// input [1] =  0; [2] =  type; [3] = step; [4] = output file
void XiXiMLL(int argc, char** argv, const int index, const int MM){


		ofstream out;
		TString outfile_name = argv[1];
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);
		int year[4];
		fit_flag = atoi(argv[2]);
		if(argc < 4 || argc > 7){
				cerr << "wrong arguments!" << endl;
				exit(1);
		}
		cout << "OK 11111111111" << endl;
		i_year.clear();
		if(argc >= 4 && argc <= 7){
				years = argc - 3;
				for(int i = 0; i < years; i++){
						year[i] = atoi(argv[i+3]);
						if(year[i] != 2009 &&  year[i] != 2012 && year[i] != 2018 && year[i] != 2019 ){
								cerr << "wrong data sets : " << year[i] << endl;
								exit(1);
						}
						i_year.push_back(year[i]);
				}
		}
		ReadData(year, index, MM);
//		for(int i = 0; i < years; i++){
//				for(int j = 0; j < 2; j++){
//						angdis[i][j]->PrintInt();
//				}
//		}
		cout << "OK 11111111113" << endl;
		// instantiating the values to be measured 
		//
		// initial values for fit
		double Jpsi_alpha       = 0.586;  	  // alpha_J/Psi 
		double Jpsi_phi       =  1.121;		  //-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		double xi_alpha     = -0.3756;  		  // alpha (Sgm->p pi0)
		double xi_phi   = 0.012;   			  // alpha (Sgm->pbar pi0)
		double xib_alpha     = 0.3756;   
		double xib_phi   = -0.012;  
		double L1_alpha     = 0.692;   
		double L2_alpha     = -0.751;   
		double L3_alpha     = 0.751;   
		double L4_alpha     = -0.692;   
		cout << "OK" << endl;
		// cout << argv[1] << endl;
		// fit nr is used to tell which analysis cuts that are used
		TMinuit *minuit=new TMinuit(10);
		Int_t ierflag=0; 
		Double_t arglist[100];
		cout << "OK 11111111111" << endl;
		minuit->SetFCN(fcnMLLG);
		cout << "OK 11111111111" << endl;
		arglist[0]= 0;
		minuit->mnexcm("SET PRINT",arglist,1,ierflag);
		arglist[0]= 0.5;
		minuit->mnexcm("SET ERR",arglist,1,ierflag);
		minuit->mnparm(0, "alpha_jpsi" ,Jpsi_alpha, 0.001, -1., 1., ierflag);
		minuit->mnparm(1, "dphi_jpsi", Jpsi_phi, 0.001, -TMath::Pi(), TMath::Pi(), ierflag);
		cout << "OK 11111111111" << endl;
		minuit->mnparm(2, "xi_alpha" , xi_alpha, 0.001, -1., 0., ierflag);
		minuit->mnparm(3, "xi_phi" , xi_phi, 0.001,  -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(4, "xib_alpha" , xib_alpha, 0.001, 0, 1., ierflag);
		minuit->mnparm(5, "xib_phi" , xib_phi, 0.001,  -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(6, "L1_alpha" , L1_alpha, 0.001, 0., 1., ierflag);
		minuit->mnparm(7, "L2_alpha" , L2_alpha, 0.001, -1., 0., ierflag);
		minuit->mnparm(8, "L3_alpha" , L3_alpha, 0.001, 1., 0., ierflag);
		minuit->mnparm(9, "L4_alpha" , L4_alpha, 0.001, -1., 0., ierflag);
		cout << "OK 11111111113" << endl;
		//minuit->mnparm(3, "A_CP" , A_CP, 0.001, 0.,0., ierflag); 
		// 	 		minuit->FixParameter(7);
		//		 minuit->FixParameter(9);
		//		 minuit->FixParameter(12);
		//		 minuit->FixParameter(13);
		cout << "OK 11111111114" << endl;
		minuit->mnexcm("MINI",arglist,0,ierflag); //minimization using the migrag
		cout << "OK 11111111115" << endl;
		//limits both 0 implies no limit 
		minuit->mnexcm("MINOS",arglist,0,ierflag);
		cout << "OK 11111111115" << endl;
		minuit->mnmatu(1);
		cout << "OK 11111111115" << endl;
		Double_t fmin, fedm, errdef;
		Int_t   npari, nparx, istat; 
		minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
		double res[10], err_res[10];
		for(int p = 0; p < 10; p++)
				minuit->GetParameter(p, res[p], err_res[p]);
		out << fmin << "," << istat << "," << NN[0][0] << ","<< NN[0][1] << ","; 
		out << res[0]<< "," << err_res[0] << "," << res[1]<< "," << err_res[1]<< ","; 
		out << res[2]<< "," << err_res[2] << "," << res[3]<< "," << err_res[3]<< ",";
		out << res[4]<< "," << err_res[4] << "," << res[5]<< "," << err_res[5]<< ",";
		out << res[6]<< "," << err_res[6] << "," << res[7]<< "," << err_res[7]<< ",";
		out << res[8]<< "," << err_res[8] << "," << res[9]<< "," << err_res[9]<< endl;
		out.close();
	//	return 0;
}

int main(int argc, char **argv){
	for(int i  = 0; i < 30; i++){
		XiXiMLL(argc, argv, i, 30);
	}
}

