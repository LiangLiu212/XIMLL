#include "TFile.h"
#include "TTree.h"
#include "readData.h"

int readData(const char *filename, AngDisXiXi *angdis, double **para, int flag, const char *type){

		bool phsp = false;
		if(strcmp(type, "PHSP") == 0) phsp = true;

		int runNo_low = 0;
		int runNo_high = 0;

		if(flag == 2009){
				runNo_low = 9800;
				runNo_high =  11000;
		}
		else if(flag == 2012){
				runNo_low = 27100;
				runNo_high = 28400;
		}
		else if(flag == 2018){
				runNo_low = 52840;
				runNo_high = 56646;
		}
		else if(flag == 2019){
				runNo_low = 56778 ;
				runNo_high = 59115;
		}

		Int_t NN1;
		std::vector<double> gD1Xithe;
		std::vector<double> gD1Lthe;
		std::vector<double> gD1Lphi;
		std::vector<double> gD1Lbthe;
		std::vector<double> gD1Lbphi;
		std::vector<double> gD1pthe;
		std::vector<double> gD1pphi;
		std::vector<double> gD1apthe;
		std::vector<double> gD1apphi;

		gD1Xithe.clear();
		gD1Lthe.clear();
		gD1Lphi.clear();
		gD1Lbthe.clear();
		gD1Lbphi.clear();
		gD1pthe.clear();
		gD1pphi.clear();
		gD1apthe.clear();
		gD1apphi.clear();

		Double_t the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi;
		int runNo;

		TFile *f1 = new TFile(filename, "read");
		TTree *t1 = (TTree*)f1->Get("xixi");
		t1->SetBranchAddress("the", &the);
		t1->SetBranchAddress("Lthe", &Lthe);
		t1->SetBranchAddress("Lphi", &Lphi);
		t1->SetBranchAddress("Lbthe", &Lbthe);
		t1->SetBranchAddress("Lbphi", &Lbphi);
		t1->SetBranchAddress("pthe", &pthe);
		t1->SetBranchAddress("pphi", &pphi);
		t1->SetBranchAddress("apthe", &apthe);
		t1->SetBranchAddress("apphi", &apphi);
		t1->SetBranchAddress("runNo", &runNo);
		int nn = 0;
		int NEvt = t1->GetEntries();
//		NEvt = 10000;
		for(int i = 0; i <  NEvt; i++){
				t1->GetEntry(i);
				if(abs(runNo) < runNo_low || abs(runNo) > runNo_high) continue;
				nn++;
				gD1Xithe.push_back(the);
				gD1Lthe.push_back(Lthe);
				gD1Lphi.push_back(Lphi);
				gD1Lbthe.push_back(Lbthe);
				gD1Lbphi.push_back(Lbphi);
				gD1pthe.push_back(pthe);
				gD1pphi.push_back(pphi);
				gD1apthe.push_back(apthe);
				gD1apphi.push_back(apphi);
				if(phsp){
						angdis->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				}
		}
		cout << filename << ", " << type << ", " << flag << ", number : " <<   nn << endl;

		NN1 = nn;

		for(int i = 0; i < NN1; i++){
				*(*(para+0)+i) = gD1Xithe[i];
				*(*(para+1)+i) = gD1Lthe[i];
				*(*(para+2)+i) = gD1Lphi[i];
				*(*(para+3)+i) = gD1Lbthe[i];
				*(*(para+4)+i) = gD1Lbphi[i];
				*(*(para+5)+i) = gD1pthe[i];
				*(*(para+6)+i) = gD1pphi[i];
				*(*(para+7)+i) = gD1apthe[i];
				*(*(para+8)+i) = gD1apphi[i];
		}
		f1->Close();
		return NN1;
}
