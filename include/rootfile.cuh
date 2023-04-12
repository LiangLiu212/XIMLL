#ifndef ROOT_FILE_CUH
#define ROOT_FILE_CUH

#include "TRandom.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TString.h"
//#include "TF1.h"
#include "TChain.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TVectorT.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TStopwatch.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TDirectory.h"
#include <TMinuit.h>
#include <vector>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <iostream>
#include "AngDisXiXi.hh"
#include "Amplitude.cuh"
#include <TApplication.h>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TArrow.h"
#include "RooFit.h"
#include <iostream>
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooPlot.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooGaussian.h"
#include "RooHistPdf.h"
#include "RooDataSet.h"
#include "RooKeysPdf.h"
#include "RooFFTConvPdf.h"
#include "RooArgusBG.h"
#include "RooArgusGauss.h"
#include "RooCmdArg.h"
#include "RooExtendPdf.h"
#include "RooFitResult.h"
#include "RooFit.h"
#include "RooArgusPoly.h"
#include "RooDataHist.h"
using namespace RooFit;
#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 6
#define MATRIX_SIZE 80
#endif

using namespace std;
class rootfile {
		public:
				rootfile() {
						NUM = 10000000;
						N_Sample = -1;
						cut_LmdDL = 0.0;
						cut_XiDL = 0.0;
						cut_XiCosTheta = 0.84;
						cut_mXi = 0.011;
						cut_mLmd1 = 0.011;
						cut_chi2kmf = 200.;
						cut_chi2Xi = 500.;
						cut_chi2Lmd = 500.;
						cut_mn2 = 0.955;
						cut_mn1 = 0.925;
						cut_ncos = 0;
						cut_pcos = 0;
						cut_nbarcos = 0;
						cut_pbarcos = 0;
						cut_deltancos = 1.0;
						cut_deltapcos = 1.0;
						cut_deltanbarcos = 1.0;
						cut_deltapbarcos = 1.0;
						fit_step=0;
						flag_fast = false;
						m_seed = 321;
						double tmpNBKG[4][2] = {{861.385, 657.545}, {4425.2, 3421.11}, {16359.6, 12720.5}, {16589.0, 12382.8}};
						for(int i = 0; i < 4; i++)
								for(int j =0; j < 2; j++)
										NBKG[i][j] = tmpNBKG[i][j];
				}
				~rootfile();
				void setCutXiDL(const double cut) {cut_XiDL = cut;}
				void setCutLmdDL(const double cut) {cut_LmdDL = cut;}
				void setCutXiCosTheta(const double cut) {cut_XiCosTheta = cut;}
				void setCutmXi(const double cut) {cut_mXi = cut;}
				void setCutmLmd(const double cut) {cut_mLmd1 = cut;}
				void setCutchi2kmf(const double cut) {cut_chi2kmf = cut;}
				void setCutchi2Xi(const double cut) {cut_chi2Xi = cut;}
				void setCutchi2Lmd(const double cut) {cut_chi2Lmd = cut;}
				void setCutmn1(const double cut) {cut_mn1 = cut;}
				void setCutmn2(const double cut) {cut_mn2 = cut;}
				void setCutncos(const double cut) {cut_ncos = cut; cut_deltancos = 0.1;}
				void setCutpcos(const double cut) {cut_pcos = cut; cut_deltapcos = 0.1;}
				void setCutnbarcos(const double cut) {cut_nbarcos = cut; cut_deltanbarcos = 0.1;}
				void setCutpbarcos(const double cut) {cut_pbarcos = cut; cut_deltapbarcos = 0.1;}

				void SetSeed(const int rdm) {m_seed = rdm;}

				void SetVersion(const TString str){
						N_Sample++;
						D_Sample[N_Sample].m_sample.clear();
						D_Sample[N_Sample].m_version = str;
				}
				void SetYear(const TString str){
						D_Sample[N_Sample].m_year = str;
				}
				void SetChannel(const TString str){
						D_Sample[N_Sample].m_channel = str;
				}


				void SetSample(TString str) { 
						cout << str << endl;
						D_Sample[N_Sample].m_sample.push_back(str);
				}

				void MassFit();
				void massFit(const int index);
				void massFitIO(const int index);
				void ReadData();
				void readData(int index, int jndex);
				void readDataIO(int index, int jndex);
				int RunHigh(const TString year);
				int RunLow(const TString year);
				bool Selection(const TString year, const TString channel, const TString sample, const int index);
				void InitialMemory();
				void FreeMemory();
				double fcnmll(double *pp);
				double cpufcnmll(double *pp);
				void Norm(const TString str){m_norm = str;};
				void Trial(const int n){i_trial = n;}
				void IOcheck(bool b){m_isIO = b;}
				void setBKG();
				double setBKGMassFit(double nn, TString year, TString channel, TString sample);
				void Print();
				void SetBKGSysTest(vector<TString> vstr);
		private:
				struct DataSample{
						TString m_year;
						TString m_channel;
						TString m_version;
						vector<TString> m_sample;

						vector<double> datamassn;
						vector<double> mdiymassn;
						vector<double> bkg1massn;
						vector<double> bkgcharged;
						vector<double> bkgetac;
						vector<double> bkgsideband;
						vector<int> NN;
						double n_bkg;
						double n_signal;
						double n_etac;
						double n_charge;
						double n_bkgerr;
						double n_signalerr;
						double n_etacerr;
						double n_chargeerr;
				};

				struct DataSample D_Sample[20];
				int fit_step;
				int i_trial;
				bool m_isIO;
				int m_seed;

				TString m_norm;
				int m_nyear;
				int N_Sample;


				TRandom *r1;

				TFile *fbkgcorr;
				TFile *fpicorr;
				TFile *fpi0corr;

				TH2D *hcorrbkg[4][2];
				TH2D *hpicorr[4][2];
				TH2D *hpi0corr[4][2];


				Int_t NUM;
				AngDisXiXi *angdis[4][2];
				Int_t NN[4][20];
				vector<double> datamassn;
				vector<double> mdiymassn;
				double **angdata[4][2][10];
				double **gpu_angdata[4][2][10];
				double *gpu_Matrix[4][2][20];
				double *gpu_amp[4][2][20];
				double *out_amp[4][2][20];
				double NBKG[4][2];
				int nBkg[4][2];
				double fra[4][2];
				double alpha[4][2];
				double fraerr[4][2];
				double alphaerr[4][2];
				double *t_sumAmp;

				double cut_LmdDL, cut_XiDL, cut_XiCosTheta, cut_mXi2, cut_mXi;
				double cut_mLmd1, cut_chi2kmf, cut_chi2Xi, cut_chi2Lmd, cut_mn2, cut_mn1;
				double cut_ncos, cut_nbarcos, cut_pcos, cut_pbarcos;
				double cut_deltancos, cut_deltapcos, cut_deltanbarcos, cut_deltapbarcos;
				Double_t m_LmdDL, m_XiDL, m_XiCosTheta, m_mXi2, m_mXi1, m_mLmd1, m_mn;
				Double_t m_chi2kmf, m_chi2Xi, m_chi2Lmd, m_angle_gam1, m_angle_gam2, m_lmd_p, m_lmd_cos;
				Double_t the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi;
				Double_t m_pion1_1_cos, m_pion1_1_pt;
				Double_t m_pion1_2_cos, m_pion1_2_pt;
				Double_t m_pion2_1_cos, m_pion2_1_pt;
				Double_t m_pion0_cos, m_pion0_rho;
				int m_runNo;
				bool flag_fast;
				int IndexChannel(TString str){
						if(!str.CompareTo("xixipm")){
								return 0;
						}
						else if(!str.CompareTo("xixipp")){
								return 1;
						}
				}
				int IndexYear(TString str){
						if(!str.CompareTo("2009")){
								return 0;
						}
						else if(!str.CompareTo("2012")){
								return 1;
						}
						else if(!str.CompareTo("2018")){
								return 2;
						}
						else if(!str.CompareTo("2019")){
								return 3;
						}
				}
				int IndexSample(TString str){
						for(int i = 0; i < D_Sample[0].m_sample.size(); i++){
								if(!str.CompareTo(D_Sample[0].m_sample[i])){
										return i;
								}
						}
						return -1;
				/*		if(!str.CompareTo("data")){
								return 0;
						}
						else if(!str.CompareTo("mdiy")){
								return 1;
						}else if(!str.CompareTo("phsp")){
								return 2;
						}else if(!str.CompareTo("bkg1")){
								return 3;
						}else if(!str.CompareTo("bkg2")){
								return 4;
						}else if(!str.CompareTo("bkg3")){
								return 5;
						}else if(!str.CompareTo("eta_c")){
								return 6;
						}else if(!str.CompareTo("charge")){
								return 7;
						}else if(!str.CompareTo("sideband")){
								return 8;
						}*/
				}

			
				Int_t CorrFactor(TString year, TString channel, TString sample);
				Double_t calculate_int(Double_t par0, Double_t par1, Double_t intlow, Double_t intup);
};
#endif
