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
using namespace std;


class rootfile {
		public:
		rootfile() {
				m_file.clear();
				m_year.clear();
				m_type.clear();
				m_sample.clear();
		}
		void Setfile(TString str) { m_file.push_back(str);}
		void Setyear(TString str) { m_year.push_back(str);}
		void Settype(TString str) { m_type.push_back(str);}
		void Setsample(TString str) { m_sample.push_back(str);}
		TString year(int n) { return m_year[n];}
		TString file(int n) { return m_file[n];}
		TString type(int n) { return m_type[n];}
		TString sample(int n) {return m_sample[n];}
		int size() {return m_file.size();}

		private:
				vector<TString> m_file;
				vector<TString> m_year;
				vector<TString> m_type;
				vector<TString> m_sample;
};

#endif
