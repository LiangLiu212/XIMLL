#include "cpufit.h"



#ifdef FOUR_CORES
void *cpufit::hFCN0(void *ptr)
{
  gLL[0] = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = *(*(angdata[iyear][0]+0)+i);
    Double_t thL1   = *(*(angdata[iyear][0]+1)+i);
    Double_t pL1    = *(*(angdata[iyear][0]+2)+i);
    Double_t thL2   = *(*(angdata[iyear][0]+3)+i);
    Double_t pL2    = *(*(angdata[iyear][0]+4)+i);
    Double_t thp1   = *(*(angdata[iyear][0]+5)+i);
    Double_t pp1    = *(*(angdata[iyear][0]+6)+i);
    Double_t thp2   = *(*(angdata[iyear][0]+7)+i);
    Double_t pp2    = *(*(angdata[iyear][0]+8)+i);

    eval = Pangdis[0]->W(thXi,thL1,pL1,thp1,pp1);
    //eval = gAngDis->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[0] += eval;
  }
}
void *cpufit::hFCN1(void *ptr)
{
  gLL[1] = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=NN2/8;
  n0=n1;n1=2*NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = *(*(angdata[iyear][0]+0)+i);
    Double_t thL1   = *(*(angdata[iyear][0]+1)+i);
    Double_t pL1    = *(*(angdata[iyear][0]+2)+i);
    Double_t thL2   = *(*(angdata[iyear][0]+3)+i);
    Double_t pL2    = *(*(angdata[iyear][0]+4)+i);
    Double_t thp1   = *(*(angdata[iyear][0]+5)+i);
    Double_t pp1    = *(*(angdata[iyear][0]+6)+i);
    Double_t thp2   = *(*(angdata[iyear][0]+7)+i);
    Double_t pp2    = *(*(angdata[iyear][0]+8)+i);
    eval = Pangdis[1]->W(thXi,thL1,pL1,thp1,pp1);
    //eval = gAngDis1->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[1] += eval;
  }
}
void *cpufit::hFCN2(void *ptr)
{
  gLL[2] = 0;
  Double_t delta   = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=2*NN2/8;
  n0=n1;n1=3*NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = *(*(angdata[iyear][0]+0)+i);
    Double_t thL1   = *(*(angdata[iyear][0]+1)+i);
    Double_t pL1    = *(*(angdata[iyear][0]+2)+i);
    Double_t thL2   = *(*(angdata[iyear][0]+3)+i);
    Double_t pL2    = *(*(angdata[iyear][0]+4)+i);
    Double_t thp1   = *(*(angdata[iyear][0]+5)+i);
    Double_t pp1    = *(*(angdata[iyear][0]+6)+i);
    Double_t thp2   = *(*(angdata[iyear][0]+7)+i);
    Double_t pp2    = *(*(angdata[iyear][0]+8)+i);

    eval = Pangdis[2]->W(thXi,thL1,pL1,thp1,pp1);
    //eval = gAngDis2->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[2] += eval;
  }
}
void *cpufit::hFCN3(void *ptr)
{
  gLL[3] = 0;
  Double_t delta   = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=3*NN2/8;
  n0=n1;n1=4*NN2/8;
  for( int i = n0; i < n1; i++ ){

    
    Double_t thXi   = gMXiTh[i];
    Double_t thL1   = gMLamTh[i];
    Double_t pL1    = gMLamPh[i];
    Double_t thp1   = gMPrTh[i];
    Double_t pp1    = gMPrPh[i];
    Double_t thL2   = gMLambarTh[i];
    Double_t pL2    = gMLambarPh[i];
    Double_t thp2   = gMPrbarTh[i];
    Double_t pp2    = gMPrbarPh[i];
    
    eval = Pangdis[3]->W(thXi,thL1,pL1,thp1,pp1);
    //eval = gAngDis3->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[3] += eval;
  }
}
void *cpufit::hFCN4(void *ptr)
{
  gLL[4] = 0;
  Double_t eval    = 0;
  Int_t n0=4*NN2/2;
  Int_t n1=5*NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = gMXiTh[i];
    Double_t thL1  	= gMLamTh[i];
    Double_t pL1 	  = gMLamPh[i];
    Double_t thp1  	= gMPrTh[i];
    Double_t pp1 	  = gMPrPh[i];
    Double_t thL2  	= gMLambarTh[i];
    Double_t pL2 	  = gMLambarPh[i];
    Double_t thp2  	= gMPrbarTh[i];
    Double_t pp2 	  = gMPrbarPh[i];

    //eval = gAngDis->W(thXi,thL1,pL1,thp1,pp1);
    eval = Nangdis[][0][0]->Wbar(thXi,thL2,pL2,thp2,pp2);
	  gLL[4] += eval;
  }
}
void *cpufit::hFCN5(void *ptr)
{
  gLL[5] = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=5*NN2/8;
  n0=n1;n1=6*NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = gMXiTh[i];
    Double_t thL1  	= gMLamTh[i];
    Double_t pL1 	  = gMLamPh[i];
    Double_t thp1  	= gMPrTh[i];
    Double_t pp1 	  = gMPrPh[i];
    Double_t thL2  	= gMLambarTh[i];
    Double_t pL2 	  = gMLambarPh[i];
    Double_t thp2  	= gMPrbarTh[i];
    Double_t pp2 	  = gMPrbarPh[i];
    
    //eval = gAngDis1->W(thXi,thL1,pL1,thp1,pp1);
    eval = gAngDisbar1->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[5] += eval;
  }
}
void *cpufit::hFCN6(void *ptr)
{
  gLL[6] = 0;
  Double_t delta   = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=6*NN2/8;
  n0=n1;n1=7*NN2/8;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = gMXiTh[i];
    Double_t thL1  	= gMLamTh[i];
    Double_t pL1 	  = gMLamPh[i];
    Double_t thp1  	= gMPrTh[i];
    Double_t pp1 	  = gMPrPh[i];
    Double_t thL2  	= gMLambarTh[i];
    Double_t pL2 	  = gMLambarPh[i];
    Double_t thp2  	= gMPrbarTh[i];
    Double_t pp2 	  = gMPrbarPh[i];
    
    //eval = gAngDis2->W(thXi,thL1,pL1,thp1,pp1);
    eval = gAngDisbar2->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[6] += eval;
  }
}
void *cpufit::hFCN7(void *ptr)
{
  gLL[7] = 0;
  Double_t delta   = 0;
  Double_t eval    = 0;
  Int_t n0=0;
  Int_t n1=7*NN2/8;
  n0=n1;n1=NN2;
  for( int i = n0; i < n1; i++ ){
    Double_t thXi   = gMXiTh[i];
    Double_t thL1  	= gMLamTh[i];
    Double_t pL1 	  = gMLamPh[i];
    Double_t thp1  	= gMPrTh[i];
    Double_t pp1 	  = gMPrPh[i];
    Double_t thL2  	= gMLambarTh[i];
    Double_t pL2 	  = gMLambarPh[i];
    Double_t thp2  	= gMPrbarTh[i];
    Double_t pp2 	  = gMPrbarPh[i];
    
    //eval = gAngDis3->W(thXi,thL1,pL1,thp1,pp1);
    eval = gAngDisbar3->Wbar(thXi,thL2,pL2,thp2,pp2);
    gLL[7] += eval;
  }
}
#endif




void cpufit::InitialMemory(){
		nsample = size()/Nyear();
		for(int i = 0; i < Nyear(); i++){
				for(int j = 0; j < nsample; j ++){
						angdata[i][j] = new double * [10];
						for(int l = 0; l < 10; l++){
								*(angdata[i][j] + l) =  new double [NUM];
						}
				}
		}
		for(int i = 0; i < Nyear(); i++){
				for(int j = 0; j < 2; j ++){
						for(int k = 0; k < 8; k++){
								Pangdis[i][j][k] = new AngDis(alpha_jpsi, dphi_jpsi, polarisation, alpha_xi , phi_xi , alpha_lam);
								Nangdis[i][j][k] = new AngDisbar(alpha_jpsi, dphi_jpsi, polarisation, alpha_xibar, phi_xibar, alpha_lambar);
						}
				}
		}
}



//=====================================================================
void readfile::IOReadData(const int index, const int MM)
{
		//	int years = 2;
		for(int i = 0; i < m_nyear; i++){ 		// read data
				double pp1[8], pp2[8];
				pp1[0] = 0.586;   pp2[0] = 0.586;
				pp1[1] = 1.213;   pp2[1] = 1.213; 
				pp1[2] = -0.375;  pp2[2] = -0.375;
				pp1[3] = 0.02;    pp2[3] = 0.02;
				pp1[4] = 0.375;   pp2[4] = 0.375;
				pp1[5] = -0.02;    pp2[5] = -0.02;
				pp1[6] = 0.692;   pp2[6] = 0.757;
				pp1[7] = -0.757;  pp2[7] = -0.692;

				for(int j = 0; j < nsample; j++){
						int l = -1;
						int n = i*nsample + j;
						if(!type(n).CompareTo("xixipm")){
								l = 0;
						}
						else if(!type(n).CompareTo("xixipp")){
								l = 1;
						}
						cout << i << "   "  << l  << endl;
						NN[i][j] = IOreadData(n, angdata[i][j], index, MM);
				}
		}
		for(int i = 0; i < m_nyear; i++){
				for(int j = 0; j < nsample; j ++){
						if(j == 0){
								cout << "N[" << i << "][" << j << "] : " << NN[i][j]  <<"  " << nBkg[i][j]<< endl;
						}
						else  if(j == nsample/2){
								cout << "N[" << i << "][" << j << "] : " << NN[i][j]  <<"  " << nBkg[i][1]<< endl;
						}
						else {
								cout << "N[" << i << "][" << j << "] : " << NN[i][j] << endl;
						}
				}
		}
}


int readfile::IOreadData(const int n, double **para, const int  index, const int MM){
		cout << n << "   " << m_file[n] << " " << m_sample[n] << " " <<  m_year[n]<<  "  " << type(n) << endl;


		bool data = false;
		bool mdiy = false;
		bool phsp = false;
		bool bkg1 = false;
		bool bkg2 = false;
		bool bkg3 = false;
		bool inc = false;
		if(!m_sample[n].CompareTo("data")) data = true;
		if(!m_sample[n].CompareTo("mdiy")) mdiy = true;
		if(!m_sample[n].CompareTo("phsp")) phsp = true;
		if(!m_sample[n].CompareTo("bkg1")) bkg1 = true;
		if(!m_sample[n].CompareTo("bkg2")) bkg2 = true;
		if(!m_sample[n].CompareTo("mdiybkg3")) bkg3 = true;
		if(!m_sample[n].CompareTo("inclusive")) inc = true;

		int runNo_low = 0;
		int runNo_high = 0;
		int iyear = -1;
		int itype = -1;
		if(!m_type[n].CompareTo("xixipm")){
				itype = 0;
		}
		else if(!m_type[n].CompareTo("xixipp")){
				itype = 1;
		}

		if(!m_year[n].CompareTo("2009")){
				iyear = 0;
				runNo_low = 9800;
				runNo_high =  11000;
		}
		else if(!m_year[n].CompareTo("2012")){
				iyear = 1;
				runNo_low = 27100;
				runNo_high = 28400;
		}
		else if(!m_year[n].CompareTo("2018")){
				iyear = 2;
				runNo_low = 52840;
				runNo_high = 56646;
		}
		else if(!m_year[n].CompareTo("2019")){
				iyear = 3;
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
		std::vector<double> gD1sW;

		gD1Xithe.clear();
		gD1Lthe.clear();
		gD1Lphi.clear();
		gD1Lbthe.clear();
		gD1Lbphi.clear();
		gD1pthe.clear();
		gD1pphi.clear();
		gD1apthe.clear();
		gD1apphi.clear();
		gD1sW.clear();

		Double_t the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi;
		Double_t m_LmdDL, m_XiDL, m_XiCosTheta, m_mXi2, m_mXi1, m_mLmd1, m_mn;
		Double_t m_chi2kmf, m_chi2Xi, m_chi2Lmd, m_angle_gam1, m_angle_gam2;
		Double_t sWeight;

		int runNo;

		TFile *f1 = new TFile(m_file[n], "read");
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
		t1->SetBranchAddress("LmdDL", &m_LmdDL);
		t1->SetBranchAddress("XiDL", &m_XiDL);
		t1->SetBranchAddress("XiCosTheta", &m_XiCosTheta);
		t1->SetBranchAddress("mXi2", &m_mXi2);
		t1->SetBranchAddress("mXi1", &m_mXi1);
		t1->SetBranchAddress("mLmd1", &m_mLmd1);
		t1->SetBranchAddress("mn", &m_mn);
		t1->SetBranchAddress("chi2kmf", &m_chi2kmf);
		t1->SetBranchAddress("chi2Xi", &m_chi2Xi);
		t1->SetBranchAddress("chi2Lmd", &m_chi2Lmd);
		t1->SetBranchAddress("angle_gam1", &m_angle_gam1);
		t1->SetBranchAddress("angle_gam2", &m_angle_gam2);
		int nn = 0;
		int NEvt = t1->GetEntries();
		int low = 0;
		int high = NEvt;
		int iEvt = NEvt/30;
		int iy  = n/nsample;
		int isample  = n%nsample;
		int isam = isample/(nsample/2);
		if(mdiy){
				low = index * iEvt;
				high = (index+1) * iEvt;
				nBkg[iy][isam] = 0;
		}
		int plow = -1;
		int phigh = -1;
		if(phsp){
				if(!m_norm.CompareTo("mdiy")){
						plow = index * iEvt;
						phigh = (index+1) * iEvt;
				}
		}

		int count1 = 0;
		int count2 = 0;
		for(int i = low; i < high; i++){
				t1->GetEntry(i);
				if(phsp){
						if(!m_norm.CompareTo("mdiy")){
								if(i >= plow && i < phigh) continue;
						}
				}
				if(abs(runNo) < runNo_low || abs(runNo) > runNo_high) continue;
				if(m_LmdDL < cut_LmdDL) continue;
				if(m_XiDL < cut_XiDL) continue;
				if(fabs(m_XiCosTheta) >  cut_XiCosTheta) continue;
				if(fabs(m_mXi2 - 1.32171) > cut_mXi) continue;
				if(fabs(m_mXi1 - 1.32171) > cut_mXi) continue;
				if(fabs(m_mLmd1 - 1.1157) > cut_mLmd1) continue;
				if(m_chi2kmf > cut_chi2kmf) continue;
				if(m_chi2Xi > cut_chi2Xi) continue;
				if(m_chi2Lmd > cut_chi2Lmd) continue;

				if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				if(m_angle_gam1 > 0.3) continue;
				if(m_angle_gam2 > 0.3) continue;

				gD1Xithe.push_back(the);
				gD1Lthe.push_back(Lthe);
				gD1Lphi.push_back(Lphi);
				gD1Lbthe.push_back(Lbthe);
				gD1Lbphi.push_back(Lbphi);
				gD1pthe.push_back(pthe);
				gD1pphi.push_back(pphi);
				gD1apthe.push_back(apthe);
				gD1apphi.push_back(apphi);
				nn++;
		}
		//	cout << infile << ", " << type << ", " << flag << ", number : " <<   nn << endl;
		NN1 = nn;
		for(int i = 0; i < NN1; i++){
				*(*(para+0)+i) = gD1Xithe[i];

				*(*(para+1)+i) = gD1Lthe[i];
				*(*(para+2)+i) = gD1Lphi[i];a

				*(*(para+3)+i) = gD1Lbthe[i];
				*(*(para+4)+i) = gD1Lbphi[i];

				*(*(para+5)+i) = gD1pthe[i];
				*(*(para+6)+i) = gD1pphi[i];

				*(*(para+7)+i) = gD1apthe[i];
				*(*(para+8)+i) = gD1apphi[i];a
		}
		f1->Close();
		cout << m_file[n] << endl;
		return NN1;
}

void cpufit::fcnMLLG(){

		for(int i = 0; i < Nyear(); i++){
				for(int k = 0; k < 8; k++){
						Pangdis[i][0][k]->Set(pp[0],pp[1],pp[2],pp[3],pp[4],pp[5]);
						Nangdis[i][0][k]->Setbar(pp[0],pp[1],pp[2],pp[6],pp[7],pp[8]);
				}
		}

		gLL[0] = 0;   gLL[1] = 0;   gLL[2] = 0;   gLL[3] = 0;   
		gLL[4] = 0;   gLL[5] = 0;   gLL[6] = 0;   gLL[7] = 0;  
		t[0]->Run();  t[1]->Run();  t[2]->Run();  t[3]->Run();   
		t[4]->Run();  t[5]->Run();  t[6]->Run();  t[7]->Run();   
		Double_t LogLike = 0;
		Double_t LogLikebar = 0;
		for(int i = 0; i < Nyear(); i++){
				for( int i = 0; i < NN1/2; i++ ){
						Double_t thXi   = *(*(angdata[i][j]+0)+i);
						Double_t thL1   = *(*(angdata[i][j]+1)+i);
						Double_t pL1    = *(*(angdata[i][j]+2)+i);
						Double_t thp1   = *(*(angdata[i][j]+5)+i);
						Double_t pp1    = *(*(angdata[i][j]+6)+i);
						Double_t eval   = gAngDisEXP->W(thXi,thL1,pL1,thp1,pp1);	
						LogLike -= TMath::Log(eval);
				}
				for( int i = 0; i < NN1/2; i++ ){
						Double_t thXi   = *(*(angdata[i][j]+0)+i);
						Double_t thL2   = *(*(angdata[i][j]+3)+i);
						Double_t pL2    = *(*(angdata[i][j]+4)+i);
						Double_t thp2   = *(*(angdata[i][j]+7)+i);
						Double_t pp2    = *(*(angdata[i][j]+8)+i);
						Double_t eval   = gAngDisEXPbar->Wbar(thXi,thL2,pL2,thp2,pp2);
						LogLikebar -= TMath::Log(eval);    

				}
		}

#ifdef FOUR_CORES	
		t[0]->Join(); t[1]->Join(); t[2]->Join(); t[3]->Join();	
		t[4]->Join(); t[5]->Join(); t[6]->Join(); t[7]->Join(); 
		Double_t norm    = (gLL[0]+gLL[1]+gLL[2]+gLL[3])/Double_t(NN2/2);
		Double_t normbar = (gLL[4]+gLL[5]+gLL[6]+gLL[7])/Double_t(NN2/2);

#endif	
		LogLike += Double_t(NN1/2)*TMath::Log(norm);
		LogLikebar += Double_t(NN1/2)*TMath::Log(normbar);
		f = LogLike + LogLikebar;
}

void cpufit::IOMLL(){
#ifdef FOUR_CORES 
		t[0] = new TThread("t0", hFCN0, (void*) 0);
		t[1] = new TThread("t1", hFCN1, (void*) 1);
		t[2] = new TThread("t2", hFCN2, (void*) 2);
		t[3] = new TThread("t3", hFCN3, (void*) 3);   
		t[4] = new TThread("t4", hFCN4, (void*) 4);
		t[5] = new TThread("t5", hFCN5, (void*) 5);
		t[6] = new TThread("t6", hFCN6, (void*) 6);
		t[7] = new TThread("t7", hFCN7, (void*) 7);   
#endif  

		TMinuit *minuit=new TMinuit(9);
		Int_t ierflag=0; 
		Double_t arglist[100];
		minuit->SetFCN(fcnMLLG);
		arglist[0]= 0;
		minuit->mnexcm("SET PRINT",arglist,1,ierflag);
		arglist[0]= 0.5;
		// these two parameters are always the same
		minuit->mnexcm("SET ERR",arglist,1,ierflag);
		minuit->mnparm(0, "alpha_jpsi" ,alpha_jpsi, 0.01, -1., 1., ierflag);
		minuit->mnparm(1, "dphi_jpsi", dphi_jpsi, 0.01, -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(2, "P_e" ,  polarisation, 0.01, -1., 1.5, ierflag);
		minuit->mnparm(3, "aXi" ,  alpha_xi,     0.001, -1., 0., ierflag);
		minuit->mnparm(4, "pXi" ,  phi_xi,       0.001, -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(5, "aL"  ,  alpha_lam,    0.001, -1., 1., ierflag);
		minuit->mnparm(6, "aXib" , alpha_xibar,  0.001, -1., 1., ierflag);
		minuit->mnparm(7, "pXib" , phi_xibar,    0.001, -TMath::Pi(), TMath::Pi(), ierflag);
		minuit->mnparm(8, "aLb"  , alpha_lambar, 0.001, -1., 1., ierflag);
		//minuit->FixParameter(3);
		//minuit->FixParameter(4);
		//minuit->FixParameter(5);
		minuit->mnexcm("MINI",arglist,0,ierflag);  //minimization using the migrag
		//limits both 0 implies no limit  
		minuit->mnexcm("MINOS",arglist,0,ierflag);
		minuit->mnmatu(1);

		double res[9], err_res[9];
		for(int p = 0; p < 9; p++)
				minuit->GetParameter(p, res[p], err_res[p]);

		Double_t 	fmin, fedm, errdef;
		Int_t 		npari, nparx, istat;	
		minuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);
		Double_t matrix[9][9];
		minuit->mnemat(&matrix[0][0],9);
}






