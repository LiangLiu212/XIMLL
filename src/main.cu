#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "AngDisXiXi.hh"
#include "Amplitude.cuh"
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
#include <time.h>
#include <fstream>
#include <string>
#include <cstdlib>
#include "readData.h"
#include "gpu_AngDisXiXi.hh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdlib>
#include <getopt.h>
#include "rootfile.h"
#include <map>
#define _SLOOW
const Int_t NUM = 10000000;
#define THREADS_PER_BLOCK 6
#define MATRIX_SIZE 80
Int_t NN[4][12];
AngDisXiXi *angdis[4][2];
double **angdata[4][12];
double **gpu_angdata[4][12];
double *gpu_Matrix[4][12];
double *gpu_amp[4][12];
double *out_amp[4][12];
static int years;
static int fit_flag = 0;
static int fit_step = 0;
static int nsample = 0;
static int readdata_index = 0;

int verbose_flag;
std::vector<int> i_year;
std::vector<int> idx;


void InitialMemory(rootfile *rf){

		nsample = rf->size()/years;

		for(int i = 0; i < years; i++){
				for(int j = 0; j < nsample; j ++){
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
}

//=====================================================================
void ReadData(rootfile *rf, const int index, const int MM)
{
		//	int years = 2;
		for(int i = 0; i < years; i++){ 		// read data
				for(int j = 0; j < nsample; j++){
						int l = -1;
						int n = i*nsample + j;
						if(!rf->type(n).CompareTo("xixipm")){
								l = 0;
						}
						else if(!rf->type(n).CompareTo("xixipp")){
								l = 1;
						}
						NN[i][j] =  readData(rf->file(n), angdis[i][l], angdata[i][j], rf->year(n), rf->sample(n), index, MM);
						if(rf->sample(n).CompareTo("phsp")){
								idx.push_back(n);
						}
				}
		}
		for(int i = 0; i < years; i++){
				for(int j = 0; j < nsample; j ++){
						cout << "N[" << i << "][" << j << "] : " << NN[i][j] << endl;
				}
		}

		double **temp_angdata[years][12]; // define a temporary array 
		for(int i = 0; i < years; i++){
				for(int j = 0; j < nsample; j ++){
						temp_angdata[i][j] = new double * [9];
						for(int l = 0; l < 9; l++){
								*(temp_angdata[i][j] + l) =  new double [NN[i][j] + THREADS_PER_BLOCK];
								for(int k = 0; k< NN[i][j]; k ++){
										*(*(temp_angdata[i][j] + l) + k) = *(*(angdata[i][j] + l) + k);
								}
						}
				}
		}
		if(readdata_index == 0){

				for(int i = 0; i < years; i++){    // copy data from cpu to gpu
						for(int j = 0; j < nsample; j ++){
								int size1 = (NN[i][j] + THREADS_PER_BLOCK) *sizeof(double);
								gpu_angdata[i][j] = new double * [9];
								for(int l = 0; l < 9; l++){
										cudaMalloc( (void **) &(*(gpu_angdata[i][j] + l)), size1 );
										cudaMemcpy( *(gpu_angdata[i][j] + l), *(temp_angdata[i][j] + l), size1, cudaMemcpyHostToDevice );
										delete [] *(temp_angdata[i][j] + l);
								}
								cudaMalloc( (void **) &gpu_Matrix[i][j], size1 * MATRIX_SIZE  );
								cudaMalloc( (void **) &gpu_amp[i][j], size1 * MATRIX_SIZE  );
								out_amp[i][j] = new double [NN[i][j] + THREADS_PER_BLOCK];
						}
				}
		}
		else{
				for(int i = 0; i < years; i++){    // copy data from cpu to gpu
						for(int j = 0; j < nsample; j ++){
								int size1 = (NN[i][j] + THREADS_PER_BLOCK) *sizeof(double);
								for(int l = 0; l < 9; l++){
										cudaMemcpy( *(gpu_angdata[i][j] + l), *(temp_angdata[i][j] + l), size1, cudaMemcpyHostToDevice );
										delete [] *(temp_angdata[i][j] + l);
								}
						}
				}
		}
		readdata_index++;
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


		cudaError_t cudaStatus;
		clock_t start,end;
		double loglike[4][12];
		for(int i = 0; i < years; i ++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
				for (int j = 0; j < (nsample - 2); j++){
						int flag = j / ((nsample - 2) / 2);
						start = clock();
						gpu_Amp <<< (NN[i][idx[j]] * MATRIX_SIZE + MATRIX_SIZE * THREADS_PER_BLOCK ) / (MATRIX_SIZE * THREADS_PER_BLOCK), MATRIX_SIZE * THREADS_PER_BLOCK >>> ( 
										*(gpu_angdata[i][idx[j]] + 0), 
										*(gpu_angdata[i][idx[j]] + 1), 
										*(gpu_angdata[i][idx[j]] + 2), 
										*(gpu_angdata[i][idx[j]] + 3), 
										*(gpu_angdata[i][idx[j]] + 4), 
										*(gpu_angdata[i][idx[j]] + 5), 
										*(gpu_angdata[i][idx[j]] + 6), 
										*(gpu_angdata[i][idx[j]] + 7), 
										*(gpu_angdata[i][idx[j]] + 8),
										gpu_amp[i][idx[j]],
										(NN[i][idx[j]] + THREADS_PER_BLOCK)*80, 
										aa_para, flag, gpu_Matrix[i][idx[j]]);
						cudaDeviceSynchronize(); // wait until prior kernel is finished
						end = clock();
					//	double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
					//	cout << "GPU 3: running kernel " << time3 << " seconds" << endl;
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 004!" << endl;
								exit(1);
						}

						int mat_size = (NN[i][idx[j]] + THREADS_PER_BLOCK) *sizeof(double);
						cudaMemcpy( out_amp[i][idx[j]], gpu_amp[i][idx[j]], mat_size, cudaMemcpyDeviceToHost );
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 002!" << endl;
								exit(1);
						}
						loglike[i][j] = 0;
						for(int evt = 0; evt < NN[i][idx[j]]; evt++){
								//	cout << "host C munu 0: " << 	*(host_eval + evt) << endl;
								if(*(out_amp[i][idx[j]] + evt) <= 0){ f=0; cout << "data : " << *(out_amp[i][idx[j]] + evt) << endl;  return; }
								loglike[i][j] += TMath::Log(*(out_amp[i][idx[j]] + evt));
						}
				}
		}
		//	exit(1);
		double norm[4][2];
		for (int i = 0; i < years; i++){
				for (int j = 0; j < 2; j++){
						norm[i][j] = 0;
						norm[i][j] = angdis[i][j]->CalcToIntegral();
				}
		}

		double N_BKG[4][2] = {{827.43333, 653.96667}, {4443, 3498}, {8827.38, 9803.34}, {8234.57, 9115.42}};
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
				if(nsample == 6){
				//		cout << loglike[i][0] << "  " << loglike[i][1] << "  " << loglike[i][2] << "  " << loglike[i][3] << endl;
						l1 = - loglike[i][0] + N_BKG[idx_year][0]*loglike[i][1]/Double_t(NN[i][1]) + (Double_t(NN[i][0]) - N_BKG[idx_year][0])*TMath::Log(norm[i][0]/Double_t(NN[i][2]));
						l2 = - loglike[i][2] + N_BKG[idx_year][1]*loglike[i][3]/Double_t(NN[i][4]) + (Double_t(NN[i][3]) - N_BKG[idx_year][1])*TMath::Log(norm[i][1]/Double_t(NN[i][5]));
				}
				else if(nsample == 4){
						l1 = - loglike[i][0] + (Double_t(NN[i][0]))*TMath::Log(norm[i][0]/NN[i][1]);
						l2 = - loglike[i][2] + (Double_t(NN[i][3]))*TMath::Log(norm[i][1]/NN[i][3]);
				}
				else if(nsample == 8){
						l1 = - loglike[i][0] - loglike[i][1] + N_BKG[idx_year][0]*loglike[i][2]/Double_t(NN[i][2]) + (Double_t(NN[i][0]) + NN[i][1] - N_BKG[idx_year][0])*TMath::Log(norm[i][0]/Double_t(NN[i][3]));
						l2 = - loglike[i][3] - loglike[i][4] + N_BKG[idx_year][1]*loglike[i][5]/Double_t(NN[i][6]) + (Double_t(NN[i][4]) + NN[i][5] - N_BKG[idx_year][1])*TMath::Log(norm[i][1]/Double_t(NN[i][7]));
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
void XiXiMLL(rootfile *rf, int index, int MM){


		ofstream out;
		TString outfile_name = "out.txt";
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);

		ReadData(rf, index, MM);



		cout << "OK 11111111113" << endl;
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

		int c;
		int m_command;
		vector<TString> m_year;
		TString m_bkg;
		TString m_inclusive;
		TString m_version;

		fit_flag = 2;

		while (1)
		{
				static struct option long_options[] =
				{
						/* These options set a flag. */
						{"data", no_argument,       0, 'd'},
						{"iocheck",   no_argument,       0, 0},
						{"inclusive",   no_argument,       0, 3},
						/* These options don’t set a flag.
						   We distinguish them by their indices. */
						{"bkg1",     no_argument,       0, 1},
						{"bkg2",  no_argument,       0, 2},
						{"version",  required_argument, 0, 'v'},
						{"mix",  no_argument, 0, 'm'},
						{"year",    required_argument, 0, 'y'},
						{0, 0, 0, 0}
				};
				/* getopt_long stores the option index here. */
				int option_index = 0;

				c = getopt_long (argc, argv, "dv:my:",
								long_options, &option_index);

				/* Detect the end of the options. */
				if (c == -1)
						break;

				switch (c)
				{
						case 'd':
								m_command = 0;
								break;
						case 0:
								m_command = 1;
								puts ("option -a\n");
								break;

						case 3: 
								m_inclusive = "inclusive";
								break;
						case 1:
								m_bkg = "bkg1";
								fit_flag = 1;
								puts ("option -b\n");
								break;

						case 2:
								m_bkg = "bkg2";
								fit_flag = 1;
								break;

						case 'v':
								m_version = optarg;
								printf ("option -d with value `%s'\n", optarg);
								break;
						case 'y':
								m_year.push_back(optarg);
								i_year.push_back(atoi(optarg));

								printf ("option -i with value `%s'\n", optarg);
								while (optind < argc && argv[optind][0] != '-'){
										m_year.push_back(argv[optind]);
										i_year.push_back(atoi(argv[optind]));
										optind++;
								}

								printf ("option -d with value `%s'\n", optarg);
								break;

						case 'm':
								break;

						case '?':

								/* getopt_long already printed an error message. */
								break;

						default:
								abort ();
				}
		}

		/* Instead of reporting ‘--verbose’
		   and ‘--brief’ as they are encountered,
		   we report the final status resulting from them. */
		if (verbose_flag)
				puts ("verbose flag is set");

		/* Print any remaining command line arguments (not options). */
		if (optind < argc)
		{
				printf ("non-option ARGV-elements: ");
				while (optind < argc)
						printf ("%s ", argv[optind++]);
				putchar ('\n');
		}



		map<TString, TString> files;
		rootfile *rf = new rootfile();
		switch (m_command){
				case 0: {
						}
				case 1: {
								TString path = "/data/liul/workarea/XIXI/fit/boost";
								TString infile;
								TString type;
								years = m_year.size();
								for(int i = 0; i < m_year.size(); i++){
										infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
										type = "mdiy";
										rf->Setfile(infile);
										rf->Setyear(m_year[i]);
										rf->Setsample(type);
										rf->Settype(Form("xixipm"));
										files.insert(make_pair(infile, type));
										if(!m_inclusive.CompareTo("inclusive")){
												infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "inclusive/inclusive.root";
												type = "inclusive";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipm"));
										}
										if(!m_bkg.CompareTo("bkg1")){
												infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg1/bkg30x.root";
												type = "bkg1";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipm"));
										}
										else if(!m_bkg.CompareTo("bkg2")){
												infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg2/bkg30x.root";
												type = "bkg2";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipm"));
										}
										infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "phsp/phsp30x.root";
										type = "phsp";
										files.insert(make_pair(infile, type));
										rf->Setfile(infile);
										rf->Setyear(m_year[i]);
										rf->Setsample(type);
										rf->Settype(Form("xixipm"));

										infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
										type = "mdiy";
										files.insert(make_pair(infile, type));
										rf->Setfile(infile);
										rf->Setyear(m_year[i]);
										rf->Setsample(type);
										rf->Settype(Form("xixipp"));
										if(!m_inclusive.CompareTo("inclusive")){
												infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "inclusive/inclusive.root";
												type = "inclusive";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipm"));
										}
										if(!m_bkg.CompareTo("bkg1")){
												infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + m_bkg +  "/bkg30x.root";
												type = "bkg1";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipp"));
										}
										else if(!m_bkg.CompareTo("bkg2")){
												infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + m_bkg+ "/bkg30x.root";
												type = "bkg2";
												files.insert(make_pair(infile, type));
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(type);
												rf->Settype(Form("xixipp"));
										}
										infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "phsp/phsp30x.root";
										type = "phsp";
										files.insert(make_pair(infile, type));
										rf->Setfile(infile);
										rf->Setyear(m_year[i]);
										rf->Setsample(type);
										rf->Settype(Form("xixipp"));
								}

						}
		}
		for(int i  = 0; i < rf->size(); i++){
				cout << rf->file(i) << " => " << rf->year(i) << '\n';
		}
		cout << endl;
		InitialMemory(rf);

		for(int i  = 0; i < 30; i++){
				XiXiMLL(rf, i, 30);
		}
		return 0;
}

