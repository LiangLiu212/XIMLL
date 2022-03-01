#ifndef ROOT_FILE_H
#define ROOT_FILE_H

#include "TRandom.h"
#include "TApplication.h"
#include "TH1F.h"
#include "TString.h"
#include "TF1.h"
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
class rootfile {
		public:
				rootfile() {
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

				}
				void Setfile(TString str) { m_file.push_back(str);}
				void Setyear(TString str) { m_year.push_back(str);}
				void Settype(TString str) { m_type.push_back(str);}
				void Setsample(TString str) { m_sample.push_back(str);}
				void Setversion(TString str) { m_version.push_back(str);}
				void SetNyear(const int nyear) {m_nyear = nyear;}
				TString year(int n) { return m_year[n];}
				TString file(int n) { return m_file[n];}
				TString type(int n) { return m_type[n];}
				TString sample(int n) {return m_sample[n];}
				TString version(int n) {return m_version[n];}
				int Nyear() {return m_nyear;}
				int size() {return m_file.size();}

				void InitialMemory();
				void FreeMemory();
				void ReadData(const int index, const int MM);
				void IOReadData(const int index, const int MM);
				void MassFit();
				double IOfcnmll(double *pp);
				double fcnmll(double *pp);
		private:
				int fit_step;
				vector<TString> m_file;
				vector<TString> m_year;
				vector<TString> m_type;
				vector<TString> m_sample;
				vector<TString> m_version;
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
				int IOreadData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM);
				int readData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM);
				double massFit(const int iyear, const int itype);
				double cut_LmdDL, cut_XiDL, cut_XiCosTheta, cut_mXi2, cut_mXi;
				double cut_mLmd1, cut_chi2kmf, cut_chi2Xi, cut_chi2Lmd, cut_mn2, cut_mn1;
};
#endif
