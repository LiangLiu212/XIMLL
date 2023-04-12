#ifndef READ_FILE_H
#define READ_FILE_H

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


#ifndef THREADS_PER_BLOCK
#define THREADS_PER_BLOCK 6
#define MATRIX_SIZE 80
#endif

using namespace std;
class readfile {
		public:
				readfile() {
						m_file.clear();
						m_year.clear();
						m_type.clear();
						m_sample.clear();
						NUM = 10000000;
						cut_LmdDL = 0.0;
						cut_XiDL = 0.0;
						cut_XiCosTheta = 0.84;
						cut_mXi = 0.011;
						cut_mLmd1 = 0.011;
						cut_chi2kmf = 200.;
						cut_chi2Xi = 999.;
						cut_chi2Lmd = 999.;
						cut_mn2 = 0.955;
						cut_mn1 = 0.925;
						fit_step=0;
						flag_fast = false;
						double tmpNBKG[4][2] = {{861.385, 657.545}, {4425.2, 3421.11}, {16359.6, 12720.5}, {16589.0, 12382.8}};
						for(int i = 0; i < 4; i++)
								for(int j =0; j < 2; j++)
										NBKG[i][j] = tmpNBKG[i][j];
				}
				virtual void setCutXiDL(const double cut) {cut_XiDL = cut;}
				virtual void setCutLmdDL(const double cut) {cut_LmdDL = cut;}
				virtual void setCutXiCosTheta(const double cut) {cut_XiCosTheta = cut;}
				virtual void setCutmXi(const double cut) {cut_mXi = cut;}
				virtual void setCutmLmd(const double cut) {cut_mLmd1 = cut;}
				virtual void setCutchi2kmf(const double cut) {cut_chi2kmf = cut;}
				virtual void setCutchi2Xi(const double cut) {cut_chi2Xi = cut;}
				virtual void setCutchi2Lmd(const double cut) {cut_chi2Lmd = cut;}
				virtual void setCutmn1(const double cut) {cut_mn1 = cut;}
				virtual void setCutmn2(const double cut) {cut_mn2 = cut;}

				virtual void Setfile(TString str) { m_file.push_back(str);}
				virtual void Setyear(TString str) { m_year.push_back(str);}
				virtual void Settype(TString str) { m_type.push_back(str);}
				virtual void Setsample(TString str) { m_sample.push_back(str);}
				virtual void Setversion(TString str) { m_version.push_back(str);}
				virtual void SetNyear(const int nyear) {m_nyear = nyear;}
				virtual void SetFast(const bool fast){flag_fast = fast;}
				virtual TString year(int n) { return m_year[n];}
				virtual int iyear(const int n){
						if(!m_year[n].CompareTo("2009")){
								return 0;
						}
						else if(!m_year[n].CompareTo("2012")){
								return 1;
						}
						else if(!m_year[n].CompareTo("2018")){
								return 2;
						}
						else if(!m_year[n].CompareTo("2019")){
								return 3;
						}
				}
				virtual TString file(int n) { return m_file[n];}
				virtual TString type(int n) { return m_type[n];}
				virtual TString sample(int n) {return m_sample[n];}
				virtual TString version(int n) {return m_version[n];}
				virtual int Nyear() {return m_nyear;}
				virtual int size() {return m_file.size();}
				virtual void SetNorm(const TString norm){
						m_norm = norm;
				}
				virtual TString Norm(){
						return m_norm;
				}

				virtual void InitialMemory();
				virtual void FreeMemory();
				virtual void IOFreeMemory();
				virtual void ReadData(const int index, const int MM);
				virtual void IOReadData(const int index, const int MM);
				virtual void MassFit(const int index, const TString m_outfile);
				virtual double IOfcnmll(double *pp);
				virtual double fcnmll(double *pp);
				virtual void readBKG(const int index);
		public:
				int fit_step;
				vector<TString> m_file;
				vector<TString> m_year;
				vector<TString> m_type;
				vector<TString> m_sample;
				vector<TString> m_version;
				TString m_norm;
				int m_nyear;
				int nsample;
				Int_t NUM;
				AngDisXiXi *angdis[4][2];
				Int_t NN[4][20];
				vector<double> datamassn;
				vector<double> mdiymassn;
				double **angdata[4][20];
				double **gpu_angdata[4][20];
				double *gpu_Matrix[4][20];
				double *gpu_amp[4][20];
				double *out_amp[4][20];
				double NBKG[4][2];
				int nBkg[4][2];
				virtual int IOreadData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM);
				virtual int readData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM);
				virtual double massFit(const int iyear, const int itype);
				double cut_LmdDL, cut_XiDL, cut_XiCosTheta, cut_mXi2, cut_mXi;
				double cut_mLmd1, cut_chi2kmf, cut_chi2Xi, cut_chi2Lmd, cut_mn2, cut_mn1;
				bool flag_fast;
};
#endif
