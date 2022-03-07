#include "rootfile.cuh"

void rootfile::InitialMemory(){
		nsample = size()/Nyear();
		for(int i = 0; i < Nyear(); i++){
				for(int j = 0; j < nsample; j ++){
						angdata[i][j] = new double * [9];
						for(int l = 0; l < 9; l++){
								*(angdata[i][j] + l) =  new double [NUM];
						}
				}
		}
		for(int i = 0; i < Nyear(); i++)
				for(int j = 0; j < 2; j ++){
								angdis[i][j] = new AngDisXiXi();
								angdis[i][j]->InitialInt();
				}
}


//=====================================================================
void rootfile::ReadData(const int index, const int MM)
{
		//	int years = 2;
		for(int i = 0; i < m_nyear; i++){ 		// read data
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
						NN[i][j] = readData(n, angdis[i][l], angdata[i][j], index, MM);
	//					if(sample(n).CompareTo("phsp")){
	//							idx.push_back(n);
	//					}
				}
		}
		for(int i = 0; i < m_nyear; i++){
				for(int j = 0; j < nsample; j ++){
						cout << "N[" << i << "][" << j << "] : " << NN[i][j] << endl;
				}
		}

		double **temp_angdata[m_nyear][20]; // define a temporary array 
		for(int i = 0; i < m_nyear; i++){
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

		for(int i = 0; i < m_nyear; i++){    // copy data from cpu to gpu
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

void rootfile::FreeMemory(){
		for(int i = 0; i < m_nyear; i++){
				for(int j = 0; j < nsample; j ++){
						for(int l = 0; l < 9; l++){
								cudaFree((*(gpu_angdata[i][j] + l)));
						}
						cudaFree(gpu_Matrix[i][j]);
						cudaFree(gpu_amp[i][j]);
						delete [] out_amp[i][j];
				}
		}
}


int rootfile::readData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM){
		cout << n << "   " << m_file[n] << " " << m_sample[n] << " " <<  m_year[n]<<  "  " << type(n) << endl;


		bool data = false;
		bool mdiy = false;
		bool phsp = false;
		bool bkg = false;
		bool inc = false;
		if(!m_sample[n].CompareTo("data")) data = true;
		if(!m_sample[n].CompareTo("mdiy")) mdiy = true;
		if(!m_sample[n].CompareTo("phsp")) phsp = true;
		if(!m_sample[n].CompareTo("bkg1")) bkg = true;
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
		Double_t m_LmdDL, m_XiDL, m_XiCosTheta, m_mXi2, m_mXi1, m_mLmd1, m_mn;
		Double_t m_chi2kmf, m_chi2Xi, m_chi2Lmd, m_angle_gam1, m_angle_gam2;

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
		//		NEvt = 10000;

		int count = 0;
		for(int i = low; i <  high; i++){
				t1->GetEntry(i);
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


				if(data){
						ang->setDataMass(count, m_mn);
						count++;
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(mdiy){
						if(m_angle_gam1 > 0.3) continue;
						if(m_angle_gam2 > 0.3) continue;
						ang->setMCMass(count, m_mn);
						count++;
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(phsp){
						if(m_angle_gam1 > 0.3) continue;
						if(m_angle_gam2 > 0.3) continue;
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(bkg){
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(inc){
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}

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
				if(phsp){
						ang->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				}
		}
		//	cout << infile << ", " << type << ", " << flag << ", number : " <<   nn << endl;

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
		cout << m_file[n] << endl;
		cout << "Finish" << endl;
		return NN1;
}

void rootfile::MassFit(){
		int iyear =-1;
		int itype = -1;
		cout << "Mass Fit" << endl;
		for(int i = 0; i < Nyear(); i++){
				for(int j =0; j < 2; j++){
				NBKG[i][j] = 0;
				NBKG[i][j]  = massFit(i, j);
				}
		}


}

double rootfile::massFit(const int iyear, const int itype){
		using namespace RooFit;
		RooRealVar mn("mn", "M(n#pi^{-}) (GeV/#font[52]{c}^{2})",0.90,0.98);
		RooDataSet signal("signal", "signal", mn);
		cout << angdis[iyear][itype]->getNmc() << endl;
		for(int i = 0; i < angdis[iyear][itype]->getNmc(); i++){
				if(i == 10000) break;
				mn = angdis[iyear][itype]->MCMass(i);
				signal.add(mn);
		}
		RooKeysPdf keysshape("keysshape", "keysshape", mn, signal, RooKeysPdf::MirrorBoth, 2);

		RooDataSet data("data", "data", mn);

		cout << angdis[iyear][itype]->getNdata() << endl;
		double dataNevt = angdis[iyear][itype]->getNdata();
		for(int i = 0; i < angdis[iyear][itype]->getNdata(); i++){
				mn = angdis[iyear][itype]->DataMass(i);
				data.add(mn);
		}


		RooRealVar mean2("mean2", "mean2", 0.005, -0.05, 0.05);	
		RooRealVar sigma2("sigma2", "sigma2",0.005, 0, 0.05);
		RooGaussian ga1("ga1", "ga1", mn, mean2, sigma2);
		mn.setBins(10000, "cache");
		RooFFTConvPdf shape("shape", "shape", mn, keysshape, ga1);


		RooRealVar m0("m0", "m0", 0.98, 0.96, 0.999);
		RooRealVar k("k", "k", -10, -50, -1.0);
		RooArgusBG argus("argus", "argus", mn, m0, k);

		RooRealVar fra1("fra1", "fra1", dataNevt*100.85, 0, 20000000);
		RooRealVar fra2("fra2", "fra2", dataNevt*10.4, 0, 20000000);
		RooAddPdf sum("sum", "sig+bak", RooArgList(shape, argus), RooArgList(fra1, fra2));

		RooFitResult *res = sum.fitTo(data, "r", Extended(), Save(1));
		mn.setRange("bkg1", 0.9, 0.98);
		mn.setRange("bkg2", cut_mn1, cut_mn2);
		RooAbsReal *intpolyX10 = argus.createIntegral(mn,  NormSet(mn), Range("bkg1"));
		RooAbsReal *intpolyX11 = argus.createIntegral(mn,  NormSet(mn), Range("bkg2"));
		double back10 = intpolyX10->getVal();
		double back11 = intpolyX11->getVal();
		cout << back10 << endl;
		cout << back11 << endl;
		cout << fra2.getVal() * back11 << endl;
		return fra2.getVal() * back11;


}

double rootfile::fcnmll(double *pp){
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
	//	clock_t start,end;
		double loglike[4][12];
		int years = Nyear();
		for(int i = 0; i < years; i ++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
				for (int j = 0; j < (nsample); j++){
						int flag = j / ((nsample) / 2);
						int n = i*nsample + j;
						if(!sample(n).CompareTo("phsp")) continue;
						if(!sample(n).CompareTo("mdiy")) continue;
					//	start = clock();
						gpu_Amp <<< (NN[i][j] * MATRIX_SIZE + MATRIX_SIZE * THREADS_PER_BLOCK ) / (MATRIX_SIZE * THREADS_PER_BLOCK), MATRIX_SIZE * THREADS_PER_BLOCK >>> ( 
										*(gpu_angdata[i][j] + 0), 
										*(gpu_angdata[i][j] + 1), 
										*(gpu_angdata[i][j] + 2), 
										*(gpu_angdata[i][j] + 3), 
										*(gpu_angdata[i][j] + 4), 
										*(gpu_angdata[i][j] + 5), 
										*(gpu_angdata[i][j] + 6), 
										*(gpu_angdata[i][j] + 7), 
										*(gpu_angdata[i][j] + 8),
										gpu_amp[i][j],
										(NN[i][j] + THREADS_PER_BLOCK)*80, 
										aa_para, flag, gpu_Matrix[i][j]);
						cudaDeviceSynchronize(); // wait until prior kernel is finished
					//	end = clock();
					//	double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
					//	cout << "GPU 3: running kernel " << time3 << " seconds" << endl;
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 004!" << endl;
								exit(1);
						}

						int mat_size = (NN[i][j] + THREADS_PER_BLOCK) *sizeof(double);
						cudaMemcpy( out_amp[i][j], gpu_amp[i][j], mat_size, cudaMemcpyDeviceToHost );
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 002!" << endl;
								exit(1);
						}
						loglike[i][j] = 0;
						if(!sample(n).CompareTo("phsp")){
								for(int evt = 0; evt < NN[i][j]; evt++){
										if(*(out_amp[i][j] + evt) <= 0){ cout << "data : " << *(out_amp[i][j] + evt) << endl;  return 0; }
										loglike[i][j] += *(out_amp[i][j] + evt);
								}

						}
						else{
								for(int evt = 0; evt < NN[i][j]; evt++){
										//	cout << "host C munu 0: " << 	*(host_eval + evt) << endl;
										if(*(out_amp[i][j] + evt) <= 0){ cout << "data : " << *(out_amp[i][j] + evt) << endl;  return 0; }
										loglike[i][j] += TMath::Log(*(out_amp[i][j] + evt));
								}
						}
				}
		}
		//	exit(1);
		double norm[4][2];
		for (int i = 0; i < years; i++){
				for (int j = 0; j < 2; j++){
						norm[i][j] = 0;
					//	angdis[i][j]->PrintInt();
						norm[i][j] = angdis[i][j]->CalcToIntegral();
				}
		}



		double llk = 0;
		double l1 = 0;
		double l2 = 0;

		for(int i = 0; i < years; i++){
	//			l1 = - loglike[i][0] + NBKG[i][0]*loglike[i][3]/Double_t(NN[i][3]) + (Double_t(NN[i][0]) - NBKG[i][0])*TMath::Log(loglike[i][2]/Double_t(NN[i][2]));
	//			l2 = - loglike[i][4] + NBKG[i][1]*loglike[i][7]/Double_t(NN[i][7]) + (Double_t(NN[i][4]) - NBKG[i][1])*TMath::Log(loglike[i][6]/Double_t(NN[i][6]));
				l1 = - loglike[i][0] + NBKG[i][0]*loglike[i][3]/Double_t(NN[i][3]) + (Double_t(NN[i][0]) - NBKG[i][0])*TMath::Log(norm[i][0]/Double_t(NN[i][2]));
				l2 = - loglike[i][4] + NBKG[i][1]*loglike[i][7]/Double_t(NN[i][7]) + (Double_t(NN[i][4]) - NBKG[i][1])*TMath::Log(norm[i][1]/Double_t(NN[i][6]));
				llk += (l1 + l2);
		}
/*
		if(fit_step%100 == 0){
				std::cout << "Loglike: " << llk << std::endl; 
				for( int i = 0; i<10 ; i++ ) cout<<pp[i]<<" ";
				cout << endl;
		}
		*/
		fit_step++;
		return llk;

}


//=====================================================================
void rootfile::IOReadData(const int index, const int MM)
{
		//	int years = 2;
		for(int i = 0; i < m_nyear; i++){ 		// read data
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
						NN[i][j] = IOreadData(n, angdis[i][l], angdata[i][j], index, MM);
	//					if(sample(n).CompareTo("phsp")){
	//							idx.push_back(n);
	//					}
				}
		}
		for(int i = 0; i < m_nyear; i++){
				for(int j = 0; j < nsample; j ++){
						cout << "N[" << i << "][" << j << "] : " << NN[i][j] << endl;
				}
		}

		double **temp_angdata[m_nyear][20]; // define a temporary array 
		for(int i = 0; i < m_nyear; i++){
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

		for(int i = 0; i < m_nyear; i++){    // copy data from cpu to gpu
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


int rootfile::IOreadData(const int n, AngDisXiXi *ang, double **para, const int  index, const int MM){
		cout << n << "   " << m_file[n] << " " << m_sample[n] << " " <<  m_year[n]<<  "  " << type(n) << endl;


		bool data = false;
		bool mdiy = false;
		bool phsp = false;
		bool bkg = false;
		bool inc = false;
		if(!m_sample[n].CompareTo("data")) data = true;
		if(!m_sample[n].CompareTo("mdiy")) mdiy = true;
		if(!m_sample[n].CompareTo("phsp")) phsp = true;
		if(!m_sample[n].CompareTo("bkg1")) bkg = true;
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
		Double_t m_LmdDL, m_XiDL, m_XiCosTheta, m_mXi2, m_mXi1, m_mLmd1, m_mn;
		Double_t m_chi2kmf, m_chi2Xi, m_chi2Lmd, m_angle_gam1, m_angle_gam2;

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
		if(mdiy){
				low = index * iEvt;
				high = (index+1) * iEvt;
		}
		//		NEvt = 10000;

		int count1 = 0;
		int count2 = 0;
		for(int i = low; i <  high; i++){
				t1->GetEntry(i);
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


				if(data){
						ang->setDataMass(count1, m_mn);
						count1++;
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(mdiy){
						ang->setDataMass(count1, m_mn);
						count1++;
						if(m_angle_gam1 > 0.3) continue;
						if(m_angle_gam2 > 0.3) continue;
						if(m_angle_gam1 < 0.3 && m_angle_gam2 < 0.3){
								ang->setMCMass(count2, m_mn);
								count2++;
						}
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(phsp){
						if(m_angle_gam1 > 0.3) continue;
						if(m_angle_gam2 > 0.3) continue;
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(bkg){
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}
				if(inc){
						if(m_mn < cut_mn1 || m_mn > cut_mn2) continue;
				}

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
				if(phsp){
						ang->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
				}
		}
		//	cout << infile << ", " << type << ", " << flag << ", number : " <<   nn << endl;

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
		cout << m_file[n] << endl;
		cout << "Finish" << endl;
		return NN1;
}




double rootfile::IOfcnmll(double *pp){
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
		//		clock_t start,end;
		double loglike[4][12];
		int years = Nyear();
		for(int i = 0; i < years; i ++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
				for (int j = 0; j < (nsample); j++){
						int flag = j / ((nsample) / 2);
						int n = i*nsample + j;
						if(!sample(n).CompareTo("phsp")) continue;
						//	if(!sample(n).CompareTo("mdiy")) continue;
						//	start = clock();
						gpu_Amp <<< (NN[i][j] * MATRIX_SIZE + MATRIX_SIZE * THREADS_PER_BLOCK ) / (MATRIX_SIZE * THREADS_PER_BLOCK), MATRIX_SIZE * THREADS_PER_BLOCK >>> ( 
										*(gpu_angdata[i][j] + 0), 
										*(gpu_angdata[i][j] + 1), 
										*(gpu_angdata[i][j] + 2), 
										*(gpu_angdata[i][j] + 3), 
										*(gpu_angdata[i][j] + 4), 
										*(gpu_angdata[i][j] + 5), 
										*(gpu_angdata[i][j] + 6), 
										*(gpu_angdata[i][j] + 7), 
										*(gpu_angdata[i][j] + 8),
										gpu_amp[i][j],
										(NN[i][j] + THREADS_PER_BLOCK)*80, 
										aa_para, flag, gpu_Matrix[i][j]);
						cudaDeviceSynchronize(); // wait until prior kernel is finished
						//	end = clock();
						//	double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
						//	cout << "GPU 3: running kernel " << time3 << " seconds" << endl;
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 004!" << endl;
								exit(1);
						}

						int mat_size = (NN[i][j] + THREADS_PER_BLOCK) *sizeof(double);
						cudaMemcpy( out_amp[i][j], gpu_amp[i][j], mat_size, cudaMemcpyDeviceToHost );
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 002!" << endl;
								exit(1);
						}
						loglike[i][j] = 0;
						if(!sample(n).CompareTo("phsp")){
								for(int evt = 0; evt < NN[i][j]; evt++){
										if(*(out_amp[i][j] + evt) <= 0){ cout << "data : " << *(out_amp[i][j] + evt) << endl;  return 0; }
										loglike[i][j] += *(out_amp[i][j] + evt);
								}

						}
						else{
								for(int evt = 0; evt < NN[i][j]; evt++){
										//	cout << "host C munu 0: " << 	*(host_eval + evt) << endl;
										if(*(out_amp[i][j] + evt) <= 0){ cout << "data : " << *(out_amp[i][j] + evt) << endl;  return 0; }
										loglike[i][j] += TMath::Log(*(out_amp[i][j] + evt));
								}
						}
				}
		}
		//	exit(1);
		double norm[4][2];
		for (int i = 0; i < years; i++){
				for (int j = 0; j < 2; j++){
						norm[i][j] = 0;
						//	angdis[i][j]->PrintInt();
						norm[i][j] = angdis[i][j]->CalcToIntegral();
				}
		}

		double llk = 0;
		double l1 = 0;
		double l2 = 0;

		for(int i = 0; i < years; i++){
				//	cout << loglike[i][0] << "		" << loglike[i][2] << "		" << loglike[i][3] << "		" << loglike[i][5] << endl;
				//	l1 = - loglike[i][0] + NBKG[i][0]*loglike[i][3]/Double_t(NN[i][3]) + (Double_t(NN[i][0]) - NBKG[i][0])*TMath::Log(loglike[i][2]/Double_t(NN[i][2]));
				//	l2 = - loglike[i][4] + NBKG[i][1]*loglike[i][7]/Double_t(NN[i][7]) + (Double_t(NN[i][4]) - NBKG[i][1])*TMath::Log(loglike[i][6]/Double_t(NN[i][6]));
				//	l1 = - loglike[i][0] + NBKG[i][0]*loglike[i][2]/Double_t(NN[i][2]) + (Double_t(NN[i][0]) - NBKG[i][0])*TMath::Log(norm[i][0]/Double_t(NN[i][1]));
				//	l2 = - loglike[i][3] + NBKG[i][1]*loglike[i][5]/Double_t(NN[i][5]) + (Double_t(NN[i][3]) - NBKG[i][1])*TMath::Log(norm[i][1]/Double_t(NN[i][4]));
				l1 = - loglike[i][0]  + (Double_t(NN[i][0]))*TMath::Log(norm[i][0]/Double_t(NN[i][1]));
				l2 = - loglike[i][3]  + (Double_t(NN[i][3]))*TMath::Log(norm[i][1]/Double_t(NN[i][4]));
				llk += (l1 + l2);
		}

		if(fit_step%100 == 0){
				std::cout << "Loglike: " << llk << std::endl; 
				for( int i = 0; i<10 ; i++ ) cout<<pp[i]<<" ";
				cout << endl;
		}
		fit_step++;
		return llk;

}

