#include "TRandom.h"
#include "TH1F.h"
//#include "TF1.h"
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
#include "myMinuit.h"
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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdlib>
#include <getopt.h>
#include "rootfile.cuh"
#include <map>
#include "Math/IFunction.h"
#include "Minuit2/Minuit2Minimizer.h"
#include "Math/Functor.h"
#include <string>
#include <iostream>
#define _SLOOW
const Int_t NUM = 10000000;
#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 6
#define MATRIX_SIZE 80
#endif
int verbose_flag;
rootfile *rf;

bool isIO;

bool fixAcpXi;
bool fixphicpXi;
bool fixAcp0;
bool fixAcpm;

bool flagcpu = false;

#include "Minimizer.cxx"
#include "Minimizeriso.cxx"

void fcnMLLG(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pp, Int_t iflag)
{
		double llf = 0.0;
		llf = rf->cpufcnmll(pp);
		if(llf == 0){ cout << "ERROR" << endl; return;}
		f =  llf;
}
// using minuit2 to minimize the likelihood function
double Rosenbrock(const double * par){
		double pp[10]={0};
		for(int i = 0; i < 10; i++){
				pp[i] = par[i];
		}
		double llf = 0.0;
		if(flagcpu){
				llf = rf->cpufcnmll(pp);
		}
		else{
				llf = rf->fcnmll(pp);
		}
		if(llf == 0){ cout << "ERROR" << endl; exit(0);}
		return llf;

}


//=====================================================================
// using minuit2 to minimize the likelihood function
void Minuit2XiXiMLL(int index, int MM, TString outfile_name = "out.txt"){

		ofstream out;
		outfile_name = "Minuit2" + outfile_name;
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);

		cout << "OK 11111111113" << endl;
		double Jpsi_alpha       = 0.586;  	  // alpha_J/Psi 
		double Jpsi_phi       =  1.121;		  //-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		double xi_alpha     = -0.3756;  		  // alpha (Sgm->p pi0)
		double xi_phi   = 0.012;   			  // alpha (Sgm->pbar pi0)
		double xib_alpha     = 0.3756;   
		double xib_phi   = -0.012;  
		double L1_alpha     = 0.672;   
		double L2_alpha     = -0.751;   
		double L3_alpha     = 0.751;   
		double L4_alpha     = -0.672;   
		cout << "OK" << endl;
		// cout << argv[1] << endl;
		// fit nr is used to tell which analysis cuts that are used

		const char * algoName = "migrad";
		const int printlevel = 1;
	//	ROOT::Math::Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer(algoName);
		ROOT::Minuit2::Minuit2Minimizer* min = new ROOT::Minuit2::Minuit2Minimizer(algoName);
		min->SetMaxFunctionCalls(1000000);
		min->SetTolerance(0.001);
		min->SetPrintLevel(printlevel);
		ROOT::Math::Functor f(&Rosenbrock,10);
		min->SetFunction(f);
		min->SetErrorDef(0.5);

		double step[10] = {0};
		for(int i = 0; i < 10;i++){
				step[i] = 0.001;
		}
		min->SetVariable(0,"alpha_jpsi", Jpsi_alpha, step[0]);
		min->SetVariable(1,"dphi_jpsi", Jpsi_phi, step[1]);
		min->SetVariable(2,"xi_alpha", xi_alpha, step[2]);
		min->SetVariable(3,"xi_phi", xi_phi, step[3]);
		min->SetVariable(4,"xib_alpha", xib_alpha, step[4]);
		min->SetVariable(5,"xib_phi", xib_phi, step[5]);
		min->SetVariable(6,"L1_alpha", L1_alpha, step[6]);
		min->SetVariable(7,"L2_alpha", L2_alpha, step[7]);
		min->SetVariable(8,"L3_alpha", L3_alpha, step[8]);
		min->SetVariable(9,"L4_alpha", L4_alpha, step[9]);
		min->SetVariableLimits(0, -1., 1.);
		min->SetVariableLimits(1, -TMath::Pi(), TMath::Pi());
		min->SetVariableLimits(2, -1., 0.);
		min->SetVariableLimits(3, -TMath::Pi(), TMath::Pi());
		min->SetVariableLimits(4, 0, 1.);
		min->SetVariableLimits(5, -TMath::Pi(), TMath::Pi());
		min->SetVariableLimits(6, 0., 1.);
		min->SetVariableLimits(7, -1., 0.);
		min->SetVariableLimits(8, 0., 1.);
		min->SetVariableLimits(9, -1., 0.);
		min->Minimize();

		
		cout << "Cov Matrix" << endl;
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						cout << Form("%10.3e", min->CovMatrix(i, j)) << "  ";
				}
				cout << endl;
		}
	
		cout << "Corr. Matrix" << endl;
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						cout << Form("%10.3f", min->Correlation(i, j)) << "  ";
				}
				cout << endl;
		}



//#define minuit2_test
#ifndef minuit2_test

		cout << "OK 11111111113" << endl;
		double res[10], err_res[10];
		double low[10], high[10];
		for(int i = 0; i < 10; i++){
				low[i] = 0;
				high[i] = 0;
				min->GetMinosError(i, low[i], high[i]);
				res[i] = min->State().Parameter(i).Value();
				err_res[i] = (fabs(low[i]) + fabs(high[i]))/2.0;
		}

		double fmin = min->State().Fval();
		int istat = min->Status();
		out << fmin << "," << istat << ","; 
		out << res[0]<< "," << err_res[0] << "," << res[1]<< "," << err_res[1]<< ","; 
		out << res[2]<< "," << err_res[2] << "," << res[3]<< "," << err_res[3]<< ",";
		out << res[4]<< "," << err_res[4] << "," << res[5]<< "," << err_res[5]<< ",";
		out << res[6]<< "," << err_res[6] << "," << res[7]<< "," << err_res[7]<< ",";
		out << res[8]<< "," << err_res[8] << "," << res[9]<< "," << err_res[9]<< endl;
	
		for(int i = 0; i < 10; i++){
				for(int j = 0; j< 10; j++){
						out << min->CovMatrix(i, j) << ",";
				}
				out << endl;
		}
		out.close();
		//	return 0;
#endif
}

//=====================================================================
// input [1] =  0; [2] =  type; [3] = step; [4] = output file
void XiXiMLL(int index, int MM, TString outfile_name = "out.txt"){

		ofstream out;
		cout << outfile_name << endl;
		outfile_name = "Minuit" + outfile_name;
		out.open(outfile_name, ios::out | ios::app);

		cout << "OK 11111111113" << endl;
		double Jpsi_alpha       = 0.586;  	  // alpha_J/Psi 
		double Jpsi_phi       =  1.121;		  //-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		double xi_alpha     = -0.3756;  		  // alpha (Sgm->p pi0)
		double xi_phi   = 0.012;   			  // alpha (Sgm->pbar pi0)
		double xib_alpha     = 0.3756;   
		double xib_phi   = -0.012;  
		double L1_alpha     = 0.672;   
		double L2_alpha     = -0.751;   
		double L3_alpha     = 0.751;   
		double L4_alpha     = -0.672;   
		cout << "OK" << endl;
		// cout << argv[1] << endl;
		// fit nr is used to tell which analysis cuts that are used
		myMinuit *minuit=new myMinuit(10);
	//	if(!isIO){
	//			minuit->setRandomSeed(3423);
	//	}
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
		out << fmin << "," << istat << ","; 
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
		int m_command = -1;
		vector<TString> m_year;
		vector<TString> m_channel;
		vector<TString> m_namesample; 		// 
		vector<TString> m_sysbkg; // test the systematic uncertainties of etac charge bkg1
		m_year.clear();
		m_channel.clear();
		m_namesample.clear();
		m_sysbkg.clear();

		TString m_version;
		TString m_normalization = "mdiy";
		TString m_outfile = "out.txt";
		int iJob1 = 0;
		int iJob2 = 30;
		bool fbkg3 = false;

		fixAcpXi = false;
		fixphicpXi = false;
		fixAcp0 = false;
		fixAcpm = false;

		rf = new rootfile();

		while (1)
		{
				static struct option long_options[] =
				{
						/* These options set a flag. */
						{"data", no_argument,       0, 'd'},

						{"iocheck",   no_argument,       0, 'i'},

						/* These options don’t set a flag.
						   We distinguish them by their indices. */
						{"bkg1",     required_argument,    0, 1001},
						{"bkg2",  no_argument,       0, 1002},
						{"bkg3",  no_argument,       0, 1003},
						{"etac",  required_argument,       0, 1004},
						{"charge",  required_argument,       0, 1005},
						{"sideband",  no_argument,       0, 1006},

						{"rdmseed",  required_argument, 0, 2001},

						{"fixAcpXi",  no_argument,      0, 3001},
						{"fixphicpXi",  no_argument,    0, 3002},
						{"fixAcp0",  no_argument,       0, 3003},
						{"fixAcpm",  no_argument,       0, 3004},
						{"cpu",  no_argument,       0, 3005},

						{"version",  required_argument, 0, 'v'},

						{"year",    required_argument, 0, 'y'},

						{"jobs",    required_argument, 0, 'j'},

						{"outfile",    required_argument, 0, 'o'},
						{"normalization",    required_argument, 0, 'n'},
						{"cutXiDL",    required_argument, 0, 1},
						{"cutLmdDL",    required_argument, 0, 2},
						{"cutXiCosTheta",    required_argument, 0, 3},
						{"cutmXi",    required_argument, 0, 4},
						{"cutmLmd",    required_argument, 0, 5},
						{"cutchi2kmf",    required_argument, 0, 6},
						{"cutchi2Xi",    required_argument, 0, 7},
						{"cutchi2Lmd",    required_argument, 0, 8},
						{"cutmn1",    required_argument, 0, 9},
						{"cutmn2",    required_argument, 0, 10},
						{"cutncos",    required_argument, 0, 11},
						{"cutpcos",    required_argument, 0, 12},
						{"cutnbarcos",    required_argument, 0, 13},
						{"cutpbarcos",    required_argument, 0, 14},
						{"fast",      no_argument, 0, 'f'},
						{"help",      no_argument, 0, 'h'},
						{0, 0, 0, 0}
				};
				/* getopt_long stores the option index here. */
				int option_index = 0;

				c = getopt_long (argc, argv, "dv:my:t:n:j:o:h",
								long_options, &option_index);


				/* Detect the end of the options. */
				if (c == -1)
						break;

				switch (c)
				{
						case 'd':
								m_command = 0;
								m_namesample.push_back("data");
								m_namesample.push_back("mdiy");
								m_namesample.push_back("phsp");
								m_channel.push_back("xixipm");
								m_channel.push_back("xixipp");
								isIO = false;
								break;
						case 'i':
								m_namesample.push_back("data");
								m_namesample.push_back("mdiy");
								m_namesample.push_back("phsp");
								m_channel.push_back("xixipm");
								m_channel.push_back("xixipp");
								m_command = 1;
								isIO = true;
								break;

						case 1001:
								m_namesample.push_back("bkg1");
								puts ("option -b\n");
								m_sysbkg.push_back(optarg);
								break;

						case 1002:
								m_namesample.push_back("bkg2");
								break;
						case 1003:
								m_namesample.push_back("bkg3");
								break;
						case 1004:
								m_namesample.push_back("etac");
								m_sysbkg.push_back(optarg);
								break;
						case 1005:
								m_namesample.push_back("charge");
								m_sysbkg.push_back(optarg);
								break;
						case 1006:
								m_namesample.push_back("sideband");
								break;
						case 2001:
								rf->SetSeed(atoi(optarg));
								cout << "SedSeed : " << atoi(optarg) << endl;
								break;

						case 3001:
								fixAcpXi = true;
								break;
						case 3002:
								fixphicpXi = true;
								break;
						case 3003:
								fixAcp0 = true;
								break;
						case 3004:
								fixAcpm = true;
								break;
						case 3005:
								flagcpu = true;
								break;


						case 'v':
								m_version = optarg;
								printf ("option -d with value `%s'\n", optarg);
								break;
						case 'y':
								m_year.push_back(optarg);

								printf ("option -i with value `%s'\n", optarg);
								while (optind < argc && argv[optind][0] != '-'){
								printf ("option -i with value `%s'\n", argv[optind]);
										m_year.push_back(argv[optind]);
										optind++;
								}

								printf ("option -d with value `%s'\n", optarg);
								break;
						case 1: 
								cout << "setCutXiDL" << endl;
								rf->setCutXiDL(atof(optarg)); break;
						case 2: 
								cout << "setCutLmdDL" << endl;
								cout << atof(optarg) << endl;
								rf->setCutLmdDL(atof(optarg)); break;
						case 3: 
								cout << "setCutXiCosTheta" << endl;
								rf->setCutXiCosTheta(atof(optarg)); break;
						case 4: 
								cout << "setCutmXi" << endl;
								rf->setCutmXi(atof(optarg)); break;
						case 5: 
								cout << "setCutmLmd" << endl;
								rf->setCutmLmd(atof(optarg)); break;
						case 6: 
								cout << "setCutchi2kmf" << endl;
								rf->setCutchi2kmf(atof(optarg)); break;
						case 7: 
								cout << "setCutchi2Xi" << endl;
								rf->setCutchi2Xi(atof(optarg)); break;
						case 8: 
								cout << "setCutchi2Lmd" << endl;
								rf->setCutchi2Lmd(atof(optarg)); break;
						case 9: 
								cout << "setCutmn1" << endl;
								rf->setCutmn1(atof(optarg)); break;
						case 10: 
								cout << "setCutmn2" << endl;
								rf->setCutmn2(atof(optarg)); break;
						case 11: 
								cout << "setCutncos" << endl;
								rf->setCutncos(atof(optarg)); break;
						case 12: 
								cout << "setCutpcos" << endl;
								rf->setCutpcos(atof(optarg)); break;
						case 13: 
								cout << "setCutnbarcos" << endl;
								rf->setCutnbarcos(atof(optarg)); break;
						case 14: 
								cout << "setCutpbarcos" << endl;
								rf->setCutpbarcos(atof(optarg)); break;
						case 'n':
								m_normalization = optarg;
								cout << "normalization : " << optarg << endl;
								break;
						case 'j':
								iJob1 = atoi(optarg);
								while (optind < argc && argv[optind][0] != '-'){
										iJob2 = atoi(argv[optind]);
										optind++;
								}
								cout << "index of Jobs : " << iJob1 << "	" << iJob2 << endl;
								break;
						case 'o': 
								m_outfile = optarg;
								break;
						case 'h':
								cout <<"--data <no argument> Fit experiment data. " << endl;
										cout << "--iocheck <no argument> Perform the Input-Output check. " << endl;
										cout << "--inclusive <no argument> Using inclusive mc to estimate the background. " << endl;
										cout << "--bkg1 <no argument> Using mdiy mc to estimate the photon contamination. " << endl;
										cout << "--bkg2 <no argument> Using sideband method to estimate the photon contamination. " << endl;
										cout << "--version <v1, truev1> Specify the version. " << endl;
										cout << "--mix <no argument> " << endl;
										cout << "--year <2009, 2012, 2018, 2019> Specify the year of experiment data taken. " << endl;
										cout << "--typea " << endl;
										cout << "--jobs <num1> <num2> the range of the job's index in IO check. " << endl;
										cout << "--outfile <file name> Specify the outfile " << endl;
										cout << "--normalization <phsp, mdiy> Using mdiy or phsp to do the normalization. " << endl;
										cout << "--cutXiDL <cut>" << endl;
										cout << "--cutLmdDL <cut>" << endl;
										cout << "--cutXiCosTheta <cut>" << endl;
										cout << "--cutmXi <cut>" << endl;
										cout << "--cutmLmd <cut>" << endl;
										cout << "--cutchi2kmf <cut>" << endl;
										cout << "--cutchi2Xi <cut>" << endl;
										cout << "--cutchi2Lmd <cut>" << endl;
										cout << "--cutmn1 <cut>" << endl;
										cout << "--cutmn2 <cut>" << endl;
										cout<<		"--fast <no argument> Read normalization from binary files." << endl;


								break;

						case '?':
								cout <<"--data <no argument> Fit experiment data. " << endl;
										cout << "--iocheck <no argument> Perform the Input-Output check. " << endl;
										cout << "--inclusive <no argument> Using inclusive mc to estimate the background. " << endl;
										cout << "--bkg1 <no argument> Using mdiy mc to estimate the photon contamination. " << endl;
										cout << "--bkg2 <no argument> Using sideband method to estimate the photon contamination. " << endl;
										cout << "--version <v1, truev1> Specify the version. " << endl;
										cout << "--mix <no argument> " << endl;
										cout << "--year <2009, 2012, 2018, 2019> Specify the year of experiment data taken. " << endl;
										cout << "--typea " << endl;
										cout << "--jobs <num1> <num2> the range of the job's index in IO check. " << endl;
										cout << "--outfile <file name> Specify the outfile " << endl;
										cout << "--normalization <phsp, mdiy> Using mdiy or phsp to do the normalization. " << endl;
										cout << "--cutXiDL <cut>" << endl;
										cout << "--cutLmdDL <cut>" << endl;
										cout << "--cutXiCosTheta <cut>" << endl;
										cout << "--cutmXi <cut>" << endl;
										cout << "--cutmLmd <cut>" << endl;
										cout << "--cutchi2kmf <cut>" << endl;
										cout << "--cutchi2Xi <cut>" << endl;
										cout << "--cutchi2Lmd <cut>" << endl;
										cout << "--cutmn1 <cut>" << endl;
										cout << "--cutmn2 <cut>" << endl;
										cout<<		"--fast <no argument> Read normalization from binary files." << endl;

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

				while (optind < argc)
						printf ("%s ", argv[optind++]);
				putchar ('\n');
		}

		rf->Norm(m_normalization);
		rf->IOcheck(isIO);
		TString path = "/data/liul/workarea/XIXI/fit/boost";
		switch (m_command){
				case 0: {
								TString infile;
								for(int i = 0; i < m_year.size(); i++){
										for(int ch = 0; ch < m_channel.size(); ch++){
												rf->SetVersion(m_version);
												rf->SetYear(m_year[i]);
												rf->SetChannel(m_channel[ch]);
												for(int j = 0; j < m_namesample.size(); j++){
														rf->SetSample(m_namesample[j]);
												}
										}
								}

								cout << endl;
								rf->InitialMemory();
								rf->ReadData();
								rf->MassFit();
								rf->SetBKGSysTest(m_sysbkg);
							//	XiXiMLL(1, 30, m_outfile);
								clock_t start,end;
								start = clock();
								Minuit2XiXiMLL(1, 30, m_outfile);
								end = clock();
								double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
								cout << "CPU 3:  " << time3 << " seconds" << endl;
								rf->Print();
								rf->FreeMemory();
								break;
						}
				case 1: {
								TString infile;
								for(int i = 0; i < m_year.size(); i++){
										for(int ch = 0; ch < m_channel.size(); ch++){
												rf->SetVersion(m_version);
												rf->SetYear(m_year[i]);
												rf->SetChannel(m_channel[ch]);
												for(int j = 0; j < m_namesample.size(); j++){
														rf->SetSample(m_namesample[j]);
												}
										}
								}

								cout << endl;
								for(int i = iJob1; i < iJob2; i++){
										rf->InitialMemory();
										if(flagcpu){
												cout << "CPU 3 S:   seconds" << endl;
										}
										else{
												cout << "GPU 3 S:   seconds" << endl;
										}
										rf->Trial(i);
										rf->ReadData();
										rf->MassFit();
									//	XiXiMLL(1, 30, m_outfile);
										clock_t start,end;
										start = clock();
										Minuit2XiXiMLL(1, 30, m_outfile);
										end = clock();
										double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
										if(flagcpu){
												cout << "CPU 3:  " << time3 << " seconds" << endl;
										}
										else{
												cout << "CPU 3:  " << time3 << " seconds" << endl;
										}
										rf->Print();
										rf->FreeMemory();
								}
								break;
						}
				case -1:{
								cout << "test !" << endl;
								break;
						}
		}
		/*
		   for(int i  = 0; i < rf->size(); i++){
		   cout << rf->file(i) << " => " << rf->year(i) << '\n';
		   }
		   cout << endl;
		   rf->InitialMemory();
		   for(int i  = 0; i < 30; i++){
		   rf->IOReadData(i, 30);
		   rf->MassFit();
		//		XiXiMLL(i, 30);
		rf->FreeMemory();
		}
		 */

		return 0;
}

