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
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cstdlib>
#include <getopt.h>
#include "rootfile.cuh"
#include <map>
#define _SLOOW
const Int_t NUM = 10000000;
#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 6
#define MATRIX_SIZE 80
#endif
int verbose_flag;
rootfile *rf;

void fcnMLLG(Int_t &npar, Double_t *gin, Double_t &f, Double_t *pp, Int_t iflag)
{
		double llf = rf->IOfcnmll(pp);
		if(llf == 0){ cout << "ERROR" << endl; return;}
		f =  llf;
}
//=====================================================================
// input [1] =  0; [2] =  type; [3] = step; [4] = output file
void XiXiMLL(int index, int MM){

		ofstream out;
		TString outfile_name = "out.txt";
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);

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
		int m_command;
		vector<TString> m_year;
		vector<TString> m_namesample; 		// 
		TString m_version;


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
						{"type",    required_argument, 0, 't'},
						{0, 0, 0, 0}
				};
				/* getopt_long stores the option index here. */
				int option_index = 0;

				c = getopt_long (argc, argv, "dv:my:t:",
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
								break;
						case 0:
								m_command = 1;
								m_namesample.push_back("mdiy");
								m_namesample.push_back("phsp");
								puts ("option -a\n");
								break;

						case 3: 
								m_namesample.push_back("inclusive");
								break;
						case 1:
								m_namesample.push_back("bkg1");
								puts ("option -b\n");
								break;

						case 2:
								m_namesample.push_back("bkg2");
								break;

						case 'v':
								m_version = optarg;
								printf ("option -d with value `%s'\n", optarg);
								break;
						case 'y':
								m_year.push_back(optarg);

								printf ("option -i with value `%s'\n", optarg);
								while (optind < argc && argv[optind][0] != '-'){
										m_year.push_back(argv[optind]);
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
			
				while (optind < argc)
						printf ("%s ", argv[optind++]);
				putchar ('\n');
		}



		rf = new rootfile();
		rf->SetNyear(m_year.size());

		TString path = "/data/liul/workarea/XIXI/fit/boost";
		switch (m_command){
				case 0: {
								TString infile;
								for(int i = 0; i < m_year.size(); i++){
										for(int j = 0; j < m_namesample.size(); j++){
												if(!m_namesample[j].CompareTo("data")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "data/data.root";
												}
												else if(!m_namesample[j].CompareTo("mdiy")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
												}
												else if(!m_namesample[j].CompareTo("phsp")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "phsp/phsp30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg1")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg1/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg2")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg2/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("inclusive")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "inclusive/inclusive.root";
												}
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(m_namesample[j]);
												rf->Settype(Form("xixipm"));
												rf->Setversion(m_version);
										}

										for(int j = 0; j < m_namesample.size(); j++){
												if(!m_namesample[j].CompareTo("data")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "data/data.root";
												}
												else if(!m_namesample[j].CompareTo("mdiy")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
												}
												else if(!m_namesample[j].CompareTo("phsp")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "phsp/phsp30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg1")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "bkg1/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg2")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "bkg2/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("inclusive")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "inclusive/inclusive.root";
												}
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(m_namesample[j]);
												rf->Settype(Form("xixipp"));
												rf->Setversion(m_version);
										}
								}
								break;



						}
				case 1: {
								TString infile;
								for(int i = 0; i < m_year.size(); i++){
										for(int j = 0; j < m_namesample.size(); j++){
												if(!m_namesample[j].CompareTo("data")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "data/data.root";
												}
												else if(!m_namesample[j].CompareTo("mdiy")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
												}
												else if(!m_namesample[j].CompareTo("phsp")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "phsp/phsp30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg1")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg1/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg2")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "bkg2/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("inclusive")){
														infile = path + "/" + m_year[i] + "/" + "xixipm" + "/" + m_version + "/" + "inclusive/inclusive.root";
												}
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(m_namesample[j]);
												rf->Settype(Form("xixipm"));
												rf->Setversion(m_version);
										}

										for(int j = 0; j < m_namesample.size(); j++){
												if(!m_namesample[j].CompareTo("data")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "data/data.root";
												}
												else if(!m_namesample[j].CompareTo("mdiy")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "mdiy/mdiy30x.root";
												}
												else if(!m_namesample[j].CompareTo("phsp")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "phsp/phsp30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg1")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "bkg1/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("bkg2")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "bkg2/bkg30x.root";
												}
												else if(!m_namesample[j].CompareTo("inclusive")){
														infile = path + "/" + m_year[i] + "/" + "xixipp" + "/" + m_version + "/" + "inclusive/inclusive.root";
												}
												rf->Setfile(infile);
												rf->Setyear(m_year[i]);
												rf->Setsample(m_namesample[j]);
												rf->Settype(Form("xixipp"));
												rf->Setversion(m_version);
										}
								}

								for(int i  = 0; i < rf->size(); i++){
										cout << rf->file(i) << " => " << rf->year(i) << '\n';
								}
								cout << endl;
								rf->InitialMemory();
								for(int i  = 0; i < 30; i++){
										rf->IOReadData(i, 30);
									//	rf->MassFit();
										XiXiMLL(i, 30);
										rf->FreeMemory();
								}
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

