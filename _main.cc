#include "TRandom.h"
#include "TH1F.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "include/AngDisXiXi.hh"
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
#define _SLOOW
Int_t NN1,NN2, NN3, NN4, NN5, NN6;
Int_t fit_nr;
AngDisXiXi *angdis1 = 0;
AngDisXiXi *angdis2 = 0;
std::vector<double> gD1Xithe;
std::vector<double> gD1Lthe;
std::vector<double> gD1Lphi;
std::vector<double> gD1Lbthe;
std::vector<double> gD1Lbphi;
std::vector<double> gD1pthe;
std::vector<double> gD1pphi;
std::vector<double> gD1apthe;
std::vector<double> gD1apphi;

std::vector<double> gM1Xithe;
std::vector<double> gM1Lthe;
std::vector<double> gM1Lphi;
std::vector<double> gM1Lbthe;
std::vector<double> gM1Lbphi;
std::vector<double> gM1pthe;
std::vector<double> gM1pphi;
std::vector<double> gM1apthe;
std::vector<double> gM1apphi;

std::vector<double> gB1Xithe;
std::vector<double> gB1Lthe;
std::vector<double> gB1Lphi;
std::vector<double> gB1Lbthe;
std::vector<double> gB1Lbphi;
std::vector<double> gB1pthe;
std::vector<double> gB1pphi;
std::vector<double> gB1apthe;
std::vector<double> gB1apphi;

std::vector<double> gD2Xithe;
std::vector<double> gD2Lthe;
std::vector<double> gD2Lphi;
std::vector<double> gD2Lbthe;
std::vector<double> gD2Lbphi;
std::vector<double> gD2pthe;
std::vector<double> gD2pphi;
std::vector<double> gD2apthe;
std::vector<double> gD2apphi;

std::vector<double> gM2Xithe;
std::vector<double> gM2Lthe;
std::vector<double> gM2Lphi;
std::vector<double> gM2Lbthe;
std::vector<double> gM2Lbphi;
std::vector<double> gM2pthe;
std::vector<double> gM2pphi;
std::vector<double> gM2apthe;
std::vector<double> gM2apphi;

std::vector<double> gB2Xithe;
std::vector<double> gB2Lthe;
std::vector<double> gB2Lphi;
std::vector<double> gB2Lbthe;
std::vector<double> gB2Lbphi;
std::vector<double> gB2pthe;
std::vector<double> gB2pphi;
std::vector<double> gB2apthe;
std::vector<double> gB2apphi;

//=====================================================================
void ReadData()
{
		gM1Xithe.clear();
		gM1Lthe.clear();
		gM1Lphi.clear();
		gM1Lbthe.clear();
		gM1Lbphi.clear();
		gM1pthe.clear();
		gM1pphi.clear();
		gM1apthe.clear();
		gM1apphi.clear();

		gD1Xithe.clear();
		gD1Lthe.clear();
		gD1Lphi.clear();
		gD1Lbthe.clear();
		gD1Lbphi.clear();
		gD1pthe.clear();
		gD1pphi.clear();
		gD1apthe.clear();
		gD1apphi.clear();

		gB1Xithe.clear();
		gB1Lthe.clear();
		gB1Lphi.clear();
		gB1Lbthe.clear();
		gB1Lbphi.clear();
		gB1pthe.clear();
		gB1pphi.clear();
		gB1apthe.clear();
		gB1apphi.clear();

		gM2Xithe.clear();
		gM2Lthe.clear();
		gM2Lphi.clear();
		gM2Lbthe.clear();
		gM2Lbphi.clear();
		gM2pthe.clear();
		gM2pphi.clear();
		gM2apthe.clear();
		gM2apphi.clear();

		gD2Xithe.clear();
		gD2Lthe.clear();
		gD2Lphi.clear();
		gD2Lbthe.clear();
		gD2Lbphi.clear();
		gD2pthe.clear();
		gD2pphi.clear();
		gD2apthe.clear();
		gD2apphi.clear();

		gB2Xithe.clear();
		gB2Lthe.clear();
		gB2Lphi.clear();
		gB2Lbthe.clear();
		gB2Lbphi.clear();
		gB2pthe.clear();
		gB2pphi.clear();
		gB2apthe.clear();
		gB2apphi.clear();

		Double_t the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi;
		Double_t n_mass, chisq;
		Int_t nn = 0;

		angdis1->InitialInt();
		angdis2->InitialInt();
		TFile *f1 = new TFile("data/XIXI_PHSP1.root", "read");
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
		for(int i = 0; i <  t1->GetEntries(); i++){
				t1->GetEntry(i);

				nn++;
				gM1Xithe.push_back(the);
				gM1Lthe.push_back(Lthe);
				gM1Lphi.push_back(Lphi);
				gM1Lbthe.push_back(Lbthe);
				gM1Lbphi.push_back(Lbphi);
				gM1pthe.push_back(pthe);
				gM1pphi.push_back(pphi);
				gM1apthe.push_back(apthe);
				gM1apphi.push_back(apphi);
				angdis1->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
		}

		cout << nn-1 << endl;

		NN1 = nn-1;
		nn = 0;

		TFile *f2 = new TFile("data/XIXI_BKG1.root", "read");
		TTree *t2 = (TTree*)f2->Get("xixi");
		t2->SetBranchAddress("the", &the);
		t2->SetBranchAddress("Lthe", &Lthe);
		t2->SetBranchAddress("Lphi", &Lphi);
		t2->SetBranchAddress("Lbthe", &Lbthe);
		t2->SetBranchAddress("Lbphi", &Lbphi);
		t2->SetBranchAddress("pthe", &pthe);
		t2->SetBranchAddress("pphi", &pphi);
		t2->SetBranchAddress("apthe", &apthe);
		t2->SetBranchAddress("apphi", &apphi);
		for(int i = 0; i <  t2->GetEntries(); i++){
				t2->GetEntry(i);
				
				nn++;
				gB1Xithe.push_back(the);
				gB1Lthe.push_back(Lthe);
				gB1Lphi.push_back(Lphi);
				gB1Lbthe.push_back(Lbthe);
				gB1Lbphi.push_back(Lbphi);
				gB1pthe.push_back(pthe);
				gB1pphi.push_back(pphi);
				gB1apthe.push_back(apthe);
				gB1apphi.push_back(apphi);
		}

		cout << nn-1 << endl;
		NN2 = nn -1;
		nn = 0;

		TFile *f3 = new TFile("data/XIXI_DATA1.root", "read");
		TTree *t3 = (TTree*)f3->Get("xixi");
		t3->SetBranchAddress("the", &the);
		t3->SetBranchAddress("Lthe", &Lthe);
		t3->SetBranchAddress("Lphi", &Lphi);
		t3->SetBranchAddress("Lbthe", &Lbthe);
		t3->SetBranchAddress("Lbphi", &Lbphi);
		t3->SetBranchAddress("pthe", &pthe);
		t3->SetBranchAddress("pphi", &pphi);
		t3->SetBranchAddress("apthe", &apthe);
		t3->SetBranchAddress("apphi", &apphi);
		for(int i = 0; i <  t3->GetEntries(); i++){
				t3->GetEntry(i);

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
		}
		cout << nn-1 << endl;
		NN3 = nn -1;
		nn=0;

		TFile *f4 = new TFile("data/XIXI_PHSP2.root", "read");
		TTree *t4 = (TTree*)f4->Get("xixi");
		t4->SetBranchAddress("the", &the);
		t4->SetBranchAddress("Lthe", &Lthe);
		t4->SetBranchAddress("Lphi", &Lphi);
		t4->SetBranchAddress("Lbthe", &Lbthe);
		t4->SetBranchAddress("Lbphi", &Lbphi);
		t4->SetBranchAddress("pthe", &pthe);
		t4->SetBranchAddress("pphi", &pphi);
		t4->SetBranchAddress("apthe", &apthe);
		t4->SetBranchAddress("apphi", &apphi);
		for(int i = 0; i <  t4->GetEntries(); i++){
				t4->GetEntry(i);

				nn++;
				gM2Xithe.push_back(the);
				gM2Lthe.push_back(Lthe);
				gM2Lphi.push_back(Lphi);
				gM2Lbthe.push_back(Lbthe);
				gM2Lbphi.push_back(Lbphi);
				gM2pthe.push_back(pthe);
				gM2pphi.push_back(pphi);
				gM2apthe.push_back(apthe);
				gM2apphi.push_back(apphi);
				angdis2->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);

		}
		cout << nn-1 << endl;


		NN4 = nn -1;
		nn = 0;
		TFile *f5 = new TFile("data/XIXI_BKG2.root", "read");
		TTree *t5 = (TTree*)f5->Get("xixi");
		t5->SetBranchAddress("the", &the);
		t5->SetBranchAddress("Lthe", &Lthe);
		t5->SetBranchAddress("Lphi", &Lphi);
		t5->SetBranchAddress("Lbthe", &Lbthe);
		t5->SetBranchAddress("Lbphi", &Lbphi);
		t5->SetBranchAddress("pthe", &pthe);
		t5->SetBranchAddress("pphi", &pphi);
		t5->SetBranchAddress("apthe", &apthe);
		t5->SetBranchAddress("apphi", &apphi);
		for(int i = 0; i <  t5->GetEntries(); i++){
				t5->GetEntry(i);
				nn++;
				gB2Xithe.push_back(the);
				gB2Lthe.push_back(Lthe);
				gB2Lphi.push_back(Lphi);
				gB2Lbthe.push_back(Lbthe);
				gB2Lbphi.push_back(Lbphi);
				gB2pthe.push_back(pthe);
				gB2pphi.push_back(pphi);
				gB2apthe.push_back(apthe);
				gB2apphi.push_back(apphi);
		}
		cout << nn-1 << endl;

		NN5 = nn -1;
		nn =0;

		TFile *f6 = new TFile("data/XIXI_DATA2.root", "read");
		TTree *t6 = (TTree*)f6->Get("xixi");
		t6->SetBranchAddress("the", &the);
		t6->SetBranchAddress("Lthe", &Lthe);
		t6->SetBranchAddress("Lphi", &Lphi);
		t6->SetBranchAddress("Lbthe", &Lbthe);
		t6->SetBranchAddress("Lbphi", &Lbphi);
		t6->SetBranchAddress("pthe", &pthe);
		t6->SetBranchAddress("pphi", &pphi);
		t6->SetBranchAddress("apthe", &apthe);
		t6->SetBranchAddress("apphi", &apphi);
		for(int i = 0; i <  t6->GetEntries(); i++){
				t6->GetEntry(i);

				nn++;
				gD2Xithe.push_back(the);
				gD2Lthe.push_back(Lthe);
				gD2Lphi.push_back(Lphi);
				gD2Lbthe.push_back(Lbthe);
				gD2Lbphi.push_back(Lbphi);
				gD2pthe.push_back(pthe);
				gD2pphi.push_back(pphi);
				gD2apthe.push_back(apthe);
				gD2apphi.push_back(apphi);
		}
		cout << nn-1 << endl;
		NN6 = nn -1;

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

	 	angdis1->SetParameter(pp1);
	 	angdis2->SetParameter(pp2);
		Double_t LogLike = 0;
		Double_t LogLikeB = 0;
		Double_t LogLikeB1 = 0;
		Double_t LogLikeB2 = 0;
		Double_t delta = 0;
		Double_t eval  = 0;

		for( int i = 0; i < NN3; i++ ){
				Double_t the  = gD1Xithe[i];
				Double_t Lthe 	 = gD1Lthe[i];
				Double_t Lphi 	 = gD1Lphi[i];
				Double_t Lbthe 	 = gD1Lbthe[i];
				Double_t Lbphi 	 = gD1Lbphi[i];
				Double_t pthe 	 = gD1pthe[i];
				Double_t pphi 	 = gD1pphi[i];
				Double_t apthe 	 = gD1apthe[i];
				Double_t apphi 	 = gD1apphi[i];

				eval = angdis1->Amp(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				if(eval <= 0){ f=0; cout << "data : " << the << endl;  return; }
				LogLike -= TMath::Log(eval);
		}


		for( int i = 0; i < NN6; i++ ){
				Double_t the  = gD2Xithe[i];
				Double_t Lthe 	 = gD2Lthe[i];
				Double_t Lphi 	 = gD2Lphi[i];
				Double_t Lbthe 	 = gD2Lbthe[i];
				Double_t Lbphi 	 = gD2Lbphi[i];
				Double_t pthe 	 = gD2pthe[i];
				Double_t pphi 	 = gD2pphi[i];
				Double_t apthe 	 = gD2apthe[i];
				Double_t apphi 	 = gD2apphi[i];

				eval = angdis2->Amp(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				if(eval <= 0){ f=0; cout << "data : " << the << endl;  return; }
				LogLike -= TMath::Log(eval);
		}




		for( int i = 0; i < NN2; i++ ){
				Double_t the  = gB1Xithe[i];
				Double_t Lthe 	 = gB1Lthe[i];
				Double_t Lphi 	 = gB1Lphi[i];
				Double_t Lbthe 	 = gB1Lbthe[i];
				Double_t Lbphi 	 = gB1Lbphi[i];
				Double_t pthe 	 = gB1pthe[i];
				Double_t pphi 	 = gB1pphi[i];
				Double_t apthe 	 = gB1apthe[i];
				Double_t apphi 	 = gB1apphi[i];

				eval = angdis1->Amp(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				if(eval <= 0){ f=0; cout << "data : " << the << endl;  return; }
				LogLikeB1 += TMath::Log(eval);
		}

		for( int i = 0; i < NN5; i++ ){
				Double_t the  = gB2Xithe[i];
				Double_t Lthe 	 = gB2Lthe[i];
				Double_t Lphi 	 = gB2Lphi[i];
				Double_t Lbthe 	 = gB2Lbthe[i];
				Double_t Lbphi 	 = gB2Lbphi[i];
				Double_t pthe 	 = gB2pthe[i];
				Double_t pphi 	 = gB2pphi[i];
				Double_t apthe 	 = gB2apthe[i];
				Double_t apphi 	 = gB2apphi[i];

				eval = angdis2->Amp(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				if(eval <= 0){ f=0; cout << "data : " << the << endl;  return; }
				LogLikeB2 += TMath::Log(eval);
		}
		
		double Nbkg1 = 19901.2;
		double Nbkg2 = 21900;
		
		LogLikeB1 = Nbkg1*LogLikeB1/Double_t(NN2);
		LogLikeB2 = Nbkg2*LogLikeB2/Double_t(NN5);
		LogLikeB = LogLikeB1 + LogLikeB2;

		Double_t norm1=0;
		norm1 = angdis1->CalcToIntegral();
		cout << " norm1 : " << norm1 << endl;
		norm1/=Double_t(NN1);
		Double_t norm2=0;
		norm2 = angdis2->CalcToIntegral();
		cout << " norm2 : " << norm2 << endl;
		norm2/=Double_t(NN4);

		cout << norm1 << " : " << norm2 << endl;
		for( int i = 0; i<10 ; i++ ) cout<<pp[i]<<" ";
		cout << endl;
		f = LogLike + LogLikeB + (NN3 - Nbkg1)*TMath::Log(norm1) + (NN6 - Nbkg2)*TMath::Log(norm2);

		cout << LogLike << " " << LogLikeB << " " << (NN3 - Nbkg1)*TMath::Log(norm1) << " " << (NN6 - Nbkg2)*TMath::Log(norm2) << endl;
		std::cout << "Loglike: " << f << std::endl; 
}
//=====================================================================
// input [1] =  0; [2] =  type; [3] = step; [4] = output file
int main(int argc, char** argv){
		ofstream out;
		TString outfile_name = argv[1];
		cout << outfile_name << endl;
		out.open(outfile_name, ios::out | ios::app);
		// instantiating the values to be measured 
		//
		Double_t pp[8];	
		for( int i = 0; i < 8; i++ ) pp[i]=1;
		// starting values for fit
		Double_t Jpsi_alpha       = 0.586;    // alpha_J/Psi 
		Double_t Jpsi_phi       =  1.121;//-TMath::Pi()/4.; // relative phase, Dphi_J/Psi
		Double_t xi_alpha     = -0.3756;   // alpha (Sgm->p pi0)
		Double_t xi_phi   = 0.012;   // alpha (Sgm->pbar pi0)
		Double_t xib_alpha     = 0.3756;   
		Double_t xib_phi   = -0.012;  
		Double_t L1_alpha     = 0.692;   
		Double_t L2_alpha     = -0.751;   
		Double_t L3_alpha     = 0.751;   
		Double_t L4_alpha     = -0.692;   
		angdis1 = new AngDisXiXi();
		angdis2 = new AngDisXiXi();
		ReadData();
		cout << "OK 11111111111" << endl;
		angdis1->PrintInt();
		angdis2->PrintInt();
		cout << "OK 11111111113" << endl;

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
		minuit->mnparm(0, "alpha_jpsi" ,Jpsi_alpha, 0.01, -1., 1., ierflag);
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
		out << fmin << "," << istat << "," << NN1 << ","<< NN2 << ","; 
		out << res[0]<< "," << err_res[0] << "," << res[1]<< "," << err_res[1]<< ","; 
		out << res[2]<< "," << err_res[2] << "," << res[3]<< "," << err_res[3]<< ",";
		out << res[4]<< "," << err_res[4] << "," << res[5]<< "," << err_res[5]<< ",";
		out << res[6]<< "," << err_res[6] << "," << res[7]<< "," << err_res[7]<< ",";
		out << res[8]<< "," << err_res[8] << "," << res[9]<< "," << err_res[9]<< endl;
		out.close();
		return 1;
}
