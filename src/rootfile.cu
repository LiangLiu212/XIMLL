#include "rootfile.cuh"
#include "TROOT.h"
#include "TMath.h"
#include "device_launch_parameters.h"
#include "cuda_runtime.h"
#include "floatreduce.h"
using namespace ROOT;

rootfile::~rootfile(){
		delete r1;
}

void rootfile::ReadData(){
		for(int i = 0; i < (N_Sample +1); i++){
				D_Sample[i].datamassn.clear();
				D_Sample[i].mdiymassn.clear();
				D_Sample[i].bkg1massn.clear();
				D_Sample[i].NN.clear();
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(m_isIO){
								readDataIO(i, j);
						}
						else {
						readData(i, j);
						}
				}
		}
		cout << "Hello 1" << endl;
        double **temp_angdata[4][2][20]; // define a temporary array
		for(int i = 0; i < (N_Sample +1); i++){
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						int iyear = IndexYear(D_Sample[i].m_year);
						int ich = IndexChannel(D_Sample[i].m_channel);
						int isample = IndexSample(D_Sample[i].m_sample[j]);
						temp_angdata[iyear][ich][isample] = new double * [9];
						for(int l = 0; l < 9; l++){
								*(temp_angdata[iyear][ich][isample] + l) = new double [D_Sample[i].NN[j] + THREADS_PER_BLOCK];
								for(int k = 0; k< D_Sample[i].NN[j]; k ++){
										*(*(temp_angdata[iyear][ich][isample] + l)+k) = *(*(angdata[iyear][ich][isample] + l) + k);
								}
						}
				}
		}

		cout << "Hello 2" << endl;
        for(int i = 0; i < (N_Sample +1); i++){    // copy data from cpu to gpu
                for(int j = 0; j < D_Sample[i].m_sample.size(); j ++){
						int iyear = IndexYear(D_Sample[i].m_year);
						int ich = IndexChannel(D_Sample[i].m_channel);
						int isample = IndexSample(D_Sample[i].m_sample[j]);

                        if(!D_Sample[i].m_sample[j].CompareTo("phsp")) continue;
                        if(!D_Sample[i].m_sample[j].CompareTo("mdiy")) continue;
                        int size1 = (D_Sample[i].NN[j] + THREADS_PER_BLOCK) *sizeof(double);
                        gpu_angdata[iyear][ich][isample] = new double * [9];
                        for(int l = 0; l < 9; l++){
                                cudaMalloc( (void **) &(*(gpu_angdata[iyear][ich][isample] + l)), size1 );
                                cudaMemcpy( *(gpu_angdata[iyear][ich][isample] + l), *(temp_angdata[iyear][ich][isample] + l), size1, cudaMemcpyHostToDevice );
                                delete [] *(temp_angdata[iyear][ich][isample] + l);
                        }
                        cudaMalloc( (void **) &gpu_Matrix[iyear][ich][isample], size1 * MATRIX_SIZE  );
                        cudaMalloc( (void **) &gpu_amp[iyear][ich][isample], size1 * MATRIX_SIZE  );
                        out_amp[iyear][ich][isample] = new double [D_Sample[i].NN[j] + THREADS_PER_BLOCK];
                }
        }
		cudaMalloc((void **)&t_sumAmp, sizeof(double));
		setBKG();

}


void rootfile::FreeMemory(){
		for(int i = 0; i < (N_Sample +1); i++){
				for(int j = 0; j < D_Sample[i].m_sample.size(); j ++){
						int iyear = IndexYear(D_Sample[i].m_year);
						int ich = IndexChannel(D_Sample[i].m_channel);
						int isample = IndexSample(D_Sample[i].m_sample[j]);
                        if(!D_Sample[i].m_sample[j].CompareTo("phsp")) continue;
                        if(!D_Sample[i].m_sample[j].CompareTo("mdiy")) continue;
                        for(int l = 0; l < 9; l++){
                                cudaFree((*(gpu_angdata[iyear][ich][isample] + l)));
                        }
                        cudaFree(gpu_Matrix[iyear][ich][isample]);
                        cudaFree(gpu_amp[iyear][ich][isample]);
                        delete [] out_amp[iyear][ich][isample];
				}
		}
}

void rootfile::readData(int index, int jndex){
		TString path = "/data/liul/workarea/XIXI/fit/boost";
		TString year = D_Sample[index].m_year;
		TString channel = D_Sample[index].m_channel;
		TString version = D_Sample[index].m_version;
		TString sample = D_Sample[index].m_sample[jndex];
		TString insample = sample;
		if(!sample.CompareTo("phsp")){
				if(!m_norm.CompareTo("mdiy")){
				insample = "mdiy";
				}
				else if(!m_norm.CompareTo("phsp")){
						insample = "phsp";
				}
		}
		if(!sample.CompareTo("bkg1")){
				insample = "mdiy";
		}
		if(!sample.CompareTo("bkg2")){
				insample = "data";
		}
		if(!sample.CompareTo("sideband")){
				insample = "data";
		}
		TString infile  = path + "/" + year + "/" + channel + "/" + version + "/" + insample + "/boost.root";
		cout << infile << " ==> " << sample << endl;
//#define read_data
#ifndef read_data
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



		TFile *f1 = new TFile(infile, "read");
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
		t1->SetBranchAddress("runNo", &m_runNo);
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
		t1->SetBranchAddress("lmd_p", &m_lmd_p);
		t1->SetBranchAddress("lmd_cos", &m_lmd_cos);
		t1->SetBranchAddress("pion1_1_cos", &m_pion1_1_cos);
		t1->SetBranchAddress("pion1_1_pt", &m_pion1_1_pt);
		t1->SetBranchAddress("pion1_2_cos", &m_pion1_2_cos);
		t1->SetBranchAddress("pion1_2_pt", &m_pion1_2_pt);
		t1->SetBranchAddress("pion2_1_cos", &m_pion2_1_cos);
		t1->SetBranchAddress("pion2_1_pt", &m_pion2_1_pt);
		t1->SetBranchAddress("pion0_cos", &m_pion0_cos);
		t1->SetBranchAddress("pion0_rho", &m_pion0_rho);

		int nn = 0;
		int NEvt = t1->GetEntries();
		int low = 0;
		int high = NEvt;
		r1->SetSeed(3251);

		for(int i = low; i <  high; i++){
				t1->GetEntry(i);
				if(!Selection(year, channel, sample, index)) continue;
				if(!sample.CompareTo("phsp")){
						int factor = CorrFactor(year, channel, sample);
						for(int iCorr = 0; iCorr < factor; iCorr++){
								gD1Xithe.push_back(the);
								gD1Lthe.push_back(Lthe);
								gD1Lphi.push_back(Lphi);
								gD1Lbthe.push_back(Lbthe);
								gD1Lbphi.push_back(Lbphi);
								gD1pthe.push_back(pthe);
								gD1pphi.push_back(pphi);
								gD1apthe.push_back(apthe);
								gD1apphi.push_back(apphi);
								if(!m_norm.CompareTo("phsp")){
										angdis[IndexYear(year)][IndexChannel(channel)]->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
								}
								else if(!m_norm.CompareTo("mdiy")){
										angdis[IndexYear(year)][IndexChannel(channel)]->AddToIntegralmDIY(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
								}
								nn++;
						}
				}
				else if(!sample.CompareTo("bkg1")){
						int factor = CorrFactor(year, channel, sample);
						for(int iCorr = 0; iCorr < factor; iCorr++){
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
				}
				else{
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
		}
		D_Sample[index].NN.push_back(nn);

		int iyear = IndexYear(year);
		int ich = IndexChannel(channel);
		int isample = IndexSample(sample);
		cout << "Hello " << iyear << ich << isample << endl;

		for(int i = 0; i < nn; i++){
				*(*(angdata[iyear][ich][isample]+0)+i) = gD1Xithe[i];
				*(*(angdata[iyear][ich][isample]+1)+i) = gD1Lthe[i];
				*(*(angdata[iyear][ich][isample]+2)+i) = gD1Lphi[i];
				*(*(angdata[iyear][ich][isample]+3)+i) = gD1Lbthe[i];
				*(*(angdata[iyear][ich][isample]+4)+i) = gD1Lbphi[i];
				*(*(angdata[iyear][ich][isample]+5)+i) = gD1pthe[i];
				*(*(angdata[iyear][ich][isample]+6)+i) = gD1pphi[i];
				*(*(angdata[iyear][ich][isample]+7)+i) = gD1apthe[i];
				*(*(angdata[iyear][ich][isample]+8)+i) = gD1apphi[i];
		}


		cout << nn << endl;
		f1->Close();
#endif
}

int rootfile::RunHigh(const TString year){
		if(!year.CompareTo("2009")){
				return 11000;
		}
		else if(!year.CompareTo("2012")){
				return 28400;
		}
		else if(!year.CompareTo("2018")){
				return 56646;
		}
		else if(!year.CompareTo("2019")){
				return 59115;
		}
}


int rootfile::RunLow(const TString year){
		if(!year.CompareTo("2009")){
				return 9800;
		}
		else if(!year.CompareTo("2012")){
				return 27100;
		}
		else if(!year.CompareTo("2018")){
				return 52840;
		}
		else if(!year.CompareTo("2019")){
				return 56778;
		}
}

bool rootfile::Selection(const TString year, const TString channel, const TString sample, const int index){
		if(!sample.CompareTo("sideband")){
				if(abs(m_runNo) > RunLow(year) && abs(m_runNo) < RunHigh(year) && m_LmdDL > cut_LmdDL && m_XiDL > cut_XiDL && fabs(m_XiCosTheta) < cut_XiCosTheta){
						if(fabs(m_mXi2 - 1.32171) < (3*cut_mXi + 0.005) && fabs(m_mXi2 - 1.32171) > (cut_mXi + 0.005) && fabs(m_mXi1 - 1.32171) < (3*cut_mXi + 0.005) && fabs(m_mXi1 - 1.32171) > (cut_mXi + 0.005)  && fabs(m_mLmd1 - 1.1157) < cut_mLmd1 ){
								if(m_chi2kmf < cut_chi2kmf && m_chi2Xi < cut_chi2Xi && m_chi2Lmd < cut_chi2Lmd){
										if(!channel.CompareTo("xixipm")){
												if(fabs(cos(pthe) - cut_ncos) < cut_deltancos && fabs(cos(apthe) - cut_pbarcos) < cut_deltapbarcos){
														D_Sample[index].bkgsideband.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}															
												}
										}
										else if(!channel.CompareTo("xixipp")){
												if(fabs(cos(pthe) - cut_pcos) < cut_deltapcos && fabs(cos(apthe) - cut_nbarcos) < cut_deltanbarcos){
														D_Sample[index].bkgsideband.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}															
												}
										}
								}
						}
				}
		}

		else if(abs(m_runNo) > RunLow(year) && abs(m_runNo) < RunHigh(year) && m_LmdDL > cut_LmdDL && m_XiDL > cut_XiDL && fabs(m_XiCosTheta) < cut_XiCosTheta){
				if(fabs(m_mXi2 - 1.32171) < cut_mXi && fabs(m_mXi1 - 1.32171) < cut_mXi && fabs(m_mLmd1 - 1.1157) < cut_mLmd1 ){
						if(m_chi2kmf < cut_chi2kmf && m_chi2Xi < cut_chi2Xi && m_chi2Lmd < cut_chi2Lmd){
								if(!channel.CompareTo("xixipm")){ 
										if(fabs(cos(pthe) - cut_ncos) < cut_deltancos && fabs(cos(apthe) - cut_pbarcos) < cut_deltapbarcos){
												if(!sample.CompareTo("data")){
														D_Sample[index].datamassn.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
												else if(!sample.CompareTo("mdiy")){
														if(m_angle_gam1 < 0.3 && m_angle_gam2 < 0.3){
																D_Sample[index].mdiymassn.push_back(m_mn);
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("phsp")){
														if(m_angle_gam1 < 0.3 && m_angle_gam2 < 0.3){
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("bkg1")){
														if(m_angle_gam1 > 0.3 || m_angle_gam2 > 0.3){
																D_Sample[index].bkg1massn.push_back(m_mn);
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("bkg2")){
														if((m_mn > 0.905 && m_mn < 0.920) || (m_mn > 0.960 && m_mn < 0.975)){
																return true;
														}
												}
												else if(!sample.CompareTo("etac")){
														D_Sample[index].bkgetac.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
												else if(!sample.CompareTo("charge")){
														D_Sample[index].bkgcharged.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
										}
								}
								else if(!channel.CompareTo("xixipp")){ 
										if(fabs(cos(pthe) - cut_pcos) < cut_deltapcos && fabs(cos(apthe) - cut_nbarcos) < cut_deltanbarcos){
												if(!sample.CompareTo("data")){
														D_Sample[index].datamassn.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
												else if(!sample.CompareTo("mdiy")){
														if(m_angle_gam1 < 0.3 && m_angle_gam2 < 0.3){
																D_Sample[index].mdiymassn.push_back(m_mn);
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("phsp")){
														if(m_angle_gam1 < 0.3 && m_angle_gam2 < 0.3){
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("bkg1")){
														if(m_angle_gam1 > 0.3 || m_angle_gam2 > 0.3){
																D_Sample[index].bkg1massn.push_back(m_mn);
																if(m_mn > cut_mn1 && m_mn < cut_mn2){
																		return true;
																}
														}
												}
												else if(!sample.CompareTo("bkg2")){
														if((m_mn > 0.905 && m_mn < 0.920) || (m_mn > 0.960 && m_mn < 0.975)){
																return true;
														}
												}
												else if(!sample.CompareTo("etac")){
														D_Sample[index].bkgetac.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
												else if(!sample.CompareTo("charge")){
														D_Sample[index].bkgcharged.push_back(m_mn);
														if(m_mn > cut_mn1 && m_mn < cut_mn2){
																return true;
														}
												}
										}
								}
						}
				}
		}
		else {
				return false;
		}
}

void rootfile::MassFit(){
		for(int i = 0; i < (N_Sample +1); i++){
				if(m_isIO){
						massFitIO(i);
				}
				else{
						massFit(i);
				}
		}
		for(int i = 0; i < (N_Sample +1); i++){
				cout << "sample " << i << " : " << D_Sample[i].n_bkg << endl;
		}
}


void rootfile::massFitIO(const int index){
		cout << index << endl;
		RooRealVar mn("mn", "M_{n} (GeV/c^{2})",0.90,0.985);
		RooDataSet signal("signal", "signal", mn);
		RooDataSet data("data", "data", mn);
		RooDataSet bkg("bkg", "bkg", mn);
		RooDataSet bkgetac("bkgetac", "bkgetac", mn);
		RooDataSet bkgcharged("bkgcharged", "bkgcharged", mn);
		RooDataSet bkgsideband("bkgsideband", "bkgsideband", mn);
		for(int i = 0; i < D_Sample[index].datamassn.size(); i++){
				mn = D_Sample[index].datamassn[i];
				data.add(mn);
		}

		for(int i = 0; i < D_Sample[index].mdiymassn.size(); i++){
				mn = D_Sample[index].mdiymassn[i];
				if(i > 10000) break;
				signal.add(mn);
		}

		for(int i = 0; i < D_Sample[index].bkg1massn.size(); i++){
				mn = D_Sample[index].bkg1massn[i];
				if(i > 10000) break;
				bkg.add(mn);
		}
		cout << "kernel estimation p.d.f of background......" << endl;
		//      RooKeysPdf keysbkg("keysbkg", "keysbkg", mn, bkg, RooKeysPdf::MirrorLeft, 1);

		//#define bkg
#ifndef bkg
		cout << "kernel estimation p.d.f of signal......" << endl;
		RooKeysPdf keysshape("keysshape", "keysshape", mn, signal, RooKeysPdf::MirrorLeft, 3);




		RooRealVar mean2("mean2", "mean2", 0.001, -0.05, 0.05);
		RooRealVar sigma2("sigma2", "sigma2",0.001, 0, 0.05);
		RooGaussian ga1("ga1", "ga1", mn, mean2, sigma2);
		mn.setBins(10000, "cache");
		RooFFTConvPdf shape("shape", "shape", mn, keysshape, ga1);


		RooRealVar p0("p0", "poly 0", 0.5, -50., 50.);
		RooRealVar p1("p1", "poly 1", 0.5, -50., 50.);
		RooRealVar p2("p2", "poly 2", -1.5, -10., 10.);
		RooRealVar p3("p3", "poly 3", 0.5, -50., 50.);
		RooRealVar p4("p4", "poly 4", 0, -900000., 900000.);
		RooRealVar m01("m01", "m01", 0.981);
		RooRealVar k1("k1", "k1", -5, -50, -0.0);
		RooRealVar pp1("pp1", "pp1", 0.5, 0, 1);
		RooArgusBG argus1("argus1", "argus1", mn, m01, k1);
		RooRealVar mean3("mean3", "mean3", 0.935, 0, 2);
		RooRealVar sigma3("sigma3", "sigma3",0.2, 0, 50);
		RooGaussian ga3("ga3", "ga3", mn, mean3, sigma3);

		RooArgusPoly argus3("argus3", "argus3", mn, m01, k1, pp1, p1, p2, p3);

		RooRealVar fra1("fra1", "fra1", 1000000, 0, 20000000);
		RooRealVar fra2("fra2", "fra2", 160000, 0, 20000000);
		RooAddPdf sum("sum", "sig+bak", RooArgList(shape, argus3), RooArgList(fra1, fra2));

		RooFitResult *resbkg = argus3.fitTo(bkg, Minimizer("Minuit2", "Migrad"));
		k1.setConstant(kTRUE);
		pp1.setConstant(kTRUE);
		p0.setConstant(kTRUE);
		p1.setConstant(kTRUE);
		p2.setConstant(kTRUE);
		p3.setConstant(kTRUE);
		RooFitResult *res = sum.fitTo(data, Extended(), Save(1), Minimizer("Minuit2", "Migrad"));
		mn.setRange("bkgrange2", cut_mn1, cut_mn2);
		mn.setRange("bkgrange1", 0.9, 0.98);
		RooAbsReal *intpolyX10 = shape.createIntegral(mn,  NormSet(mn), Range("bkgrange2"));
		RooAbsReal *intpolyX11 = argus3.createIntegral(mn,  NormSet(mn), Range("bkgrange2"));
		D_Sample[index].n_signal = intpolyX10->getVal()*fra1.getVal();
		D_Sample[index].n_bkg = intpolyX11->getVal()*fra2.getVal();
		D_Sample[index].n_signalerr = fra1.getError();
		D_Sample[index].n_bkgerr = fra2.getError();
		/*

		   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

		   RooPlot *xframe = mn.frame(Title("Mass Fit p.d.f."));

		   data.plotOn(xframe);

		   sum.plotOn(xframe, LineColor(4));

		   sum.plotOn(xframe, Components("keysshape"), LineColor(2), LineStyle(2));

		   sum.plotOn(xframe, Components("argus3"), LineColor(6), LineStyle(3));
		   sum.plotOn(xframe, Components("keysbkgetac"), LineColor(7), LineStyle(3));

		   sum.plotOn(xframe, Components("keysbkgcharged"), LineColor(8), LineStyle(3));
		   sum.plotOn(xframe, Components("keysbkgsideband"), LineColor(8), LineStyle(3));
		   xframe->Draw();
		   c1->SaveAs("massfit2018.eps");
		 */

#endif
}




void rootfile::massFit(const int index){
		cout << index << endl;
		double inisignal = 0;
		double inibkg = 0;
		RooRealVar mn("mn", "M_{n} (GeV/c^{2})",0.90,0.985);
		RooDataSet signal("signal", "signal", mn);
		RooDataSet data("data", "data", mn);
		RooDataSet bkg("bkg", "bkg", mn);
		RooDataSet bkgetac("bkgetac", "bkgetac", mn);
		RooDataSet bkgcharged("bkgcharged", "bkgcharged", mn);
		RooDataSet bkgsideband("bkgsideband", "bkgsideband", mn);
		inisignal = D_Sample[index].datamassn.size()*0.825 * 1.2;
		inibkg = D_Sample[index].datamassn.size()*0.125 ;
		for(int i = 0; i < D_Sample[index].datamassn.size(); i++){
				mn = D_Sample[index].datamassn[i];
				data.add(mn);
		}

		for(int i = 0; i < D_Sample[index].mdiymassn.size(); i++){
				mn = D_Sample[index].mdiymassn[i];
				if(i == 10000) break;
				signal.add(mn);
		}

		for(int i = 0; i < D_Sample[index].bkg1massn.size(); i++){
				mn = D_Sample[index].bkg1massn[i];
				if(i == 15000) break;
				bkg.add(mn);
		}
		cout << "kernel estimation p.d.f of background......" << endl;
		//      RooKeysPdf keysbkg("keysbkg", "keysbkg", mn, bkg, RooKeysPdf::MirrorLeft, 1);

		double f_etac = 0, f_charged = 0, f_sideband = 0;
		for(int i = 0; i < D_Sample[index].bkgetac.size(); i++){
				mn = D_Sample[index].bkgetac[i];
				if(D_Sample[index].bkgetac[i] < 0.90) continue;
				f_etac += 1.0;
				bkgetac.add(mn);
		}
		f_etac = setBKGMassFit(f_etac, D_Sample[index].m_year, D_Sample[index].m_channel, "etac");
		cout << "kernel estimation p.d.f of background......" << endl;
		RooKeysPdf keysbkgetac("keysbkgetac", "keysbkgetac", mn, bkgetac, RooKeysPdf::MirrorBoth, 1);
		for(int i = 0; i < D_Sample[index].bkgcharged.size(); i++){
				mn = D_Sample[index].bkgcharged[i];
				if(D_Sample[index].bkgcharged[i] < 0.90) continue;
				f_charged += 1.0;
				bkgcharged.add(mn);
		}
		f_charged = setBKGMassFit(f_charged, D_Sample[index].m_year, D_Sample[index].m_channel, "charge");

		cout << "kernel estimation p.d.f of background......" << endl;
		RooKeysPdf keysbkgcharged("keysbkgcharged", "keysbkgcharged", mn, bkgcharged, RooKeysPdf::MirrorBoth, 1);


		for(int i = 0; i < D_Sample[index].bkgsideband.size(); i++){
				mn = D_Sample[index].bkgsideband[i];
				if(D_Sample[index].bkgsideband[i] < 0.90) continue;
				f_sideband += 1.0;
				bkgsideband.add(mn);
		}
		f_sideband = f_sideband * 0.25;  // FIXME: Just a test 
		cout << "kernel estimation p.d.f of background......" << endl;
		RooKeysPdf keysbkgsideband("keysbkgsideband", "keysbkgsideband", mn, bkgsideband, RooKeysPdf::MirrorBoth, 1);

		//#define bkg
#ifndef bkg
		cout << "kernel estimation p.d.f of signal......" << endl;
		RooKeysPdf keysshape("keysshape", "keysshape", mn, signal, RooKeysPdf::MirrorLeft, 3);




		RooRealVar mean2("mean2", "mean2", 0.001, -0.05, 0.05);
		RooRealVar sigma2("sigma2", "sigma2",0.001, 0, 0.05);
		RooGaussian ga1("ga1", "ga1", mn, mean2, sigma2);
		//		mn.setBins(10000, "cache");
		RooFFTConvPdf shape("shape", "shape", mn, keysshape, ga1);


		RooRealVar p0("p0", "poly 0", 0.5, -5000., 5000.);
		RooRealVar p1("p1", "poly 1", 0.5, -5000., 5000.);
		RooRealVar p2("p2", "poly 2", 0.5, -5000., 5000.);
		RooRealVar p3("p3", "poly 3", 0.5, -5000., 5000.);
		RooRealVar p4("p4", "poly 4", 0, -900000., 900000.);
		RooRealVar m01("m01", "m01", 0.981);
		RooRealVar k1("k1", "k1", -5, -10, -0.0);
		RooRealVar pp1("pp1", "pp1", 0.5, 0.0, 1.0);
		RooArgusBG argus1("argus1", "argus1", mn, m01, k1);
		RooRealVar mean3("mean3", "mean3", 0.935, 0, 2);
		RooRealVar sigma3("sigma3", "sigma3",0.2, 0, 50);
		RooGaussian ga3("ga3", "ga3", mn, mean3, sigma3);

		RooArgusPoly argus3("argus3", "argus3", mn, m01, k1, pp1, p1, p2, p3);

		RooRealVar fra1("fra1", "fra1", inisignal, 0, 4000000);
		RooRealVar fra2("fra2", "fra2", inibkg, 0, 400000);
		RooRealVar fra3("fra3", "fra3", f_etac);
		RooRealVar fra4("fra4", "fra4", f_charged);
		RooRealVar fra5("fra5", "fra5", f_sideband);
		RooAddPdf sum("sum", "sig+bak", RooArgList(shape, argus3, keysbkgetac, keysbkgcharged, keysbkgsideband), RooArgList(fra1, fra2, fra3, fra4, fra5));
		RooAddPdf sum2("sum2", "sig+bak", RooArgList(shape, argus3, keysbkgetac, keysbkgcharged), RooArgList(fra1, fra2, fra3, fra4));

		RooFitResult *resbkg = argus3.fitTo(bkg);
		//	resbkg->Print("v");
		cout << "Finish model background" << endl;

		k1.setConstant(kTRUE);
		pp1.setConstant(kTRUE);
		p0.setConstant(kTRUE);
		p1.setConstant(kTRUE);
		p2.setConstant(kTRUE);
		p3.setConstant(kTRUE);

		cout << "Number of background ==> " << f_etac << "  " << f_charged << "  " << f_sideband << endl;
		if(IndexSample("sideband") != -1) {
				RooFitResult *res = sum.fitTo(data,  Extended(), Save(1), Minimizer("Minuit2", "Migrad"));
				int fitstatus = res->status();
				int fit_iter = 0;  // the count of fit if status > 1
				while(fitstatus > 1){
					//	fra1.setVal(fra1.getVal());
					//	fra2.setVal(fra2.getVal());
						res = sum.fitTo(data,  Extended(), Save(1), Minimizer("Minuit2", "Migrad"));
						fitstatus = res->status();
						cout << "Status : " << fitstatus << endl;
						fit_iter++;
						if(fit_iter > 5) break;
				}
				//	res->Print("v");
				cout << "Finish mass background  " << D_Sample[index].m_year << "  " << D_Sample[index].m_channel << "  " << res->status() << endl;
		}
		else{
				RooFitResult *res = sum2.fitTo(data,  Extended(), Save(1), Minimizer("Minuit2", "Migrad"));
		}
		mn.setRange("bkgrange2", cut_mn1, cut_mn2);
		mn.setRange("bkgrange1", 0.9, 0.98);
		RooAbsReal *intpolyX10 = shape.createIntegral(mn,  NormSet(mn), Range("bkgrange2"));
		RooAbsReal *intpolyX11 = argus3.createIntegral(mn,  NormSet(mn), Range("bkgrange2"));
		D_Sample[index].n_signal = intpolyX10->getVal()*fra1.getVal();
		D_Sample[index].n_bkg = intpolyX11->getVal()*fra2.getVal();
		D_Sample[index].n_signalerr = fra1.getError();
		D_Sample[index].n_bkgerr = fra2.getError();

		/*
		   TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

		   RooPlot *xframe = mn.frame(Title("Mass Fit p.d.f."));

		   data.plotOn(xframe);

		   sum.plotOn(xframe, LineColor(4));

		   sum.plotOn(xframe, Components("keysshape"), LineColor(2), LineStyle(2));

		   sum.plotOn(xframe, Components("argus3"), LineColor(6), LineStyle(3));
		   sum.plotOn(xframe, Components("keysbkgetac"), LineColor(7), LineStyle(3));

		   sum.plotOn(xframe, Components("keysbkgcharged"), LineColor(8), LineStyle(3));
		   sum.plotOn(xframe, Components("keysbkgsideband"), LineColor(8), LineStyle(3));
		   xframe->Draw();
		   c1->SaveAs("/home/liul/workarea/XiXi/MLL/XIXIRUN/RUNDATA/v9/massfit2018.eps");
		 */


#endif
}

// ==================  GPU MEMORY =============

void rootfile::InitialMemory(){
		for(int i = 0; i < 4; i++){
				for(int ch = 0; ch < 2; ch++){
						for(int j = 0; j < 10; j ++){
								angdata[i][ch][j] = new double * [10];
								for(int l = 0; l < 10; l++){
										*(angdata[i][ch][j] + l) =  new double [NUM];
								}
						}
				}
		}
		TString s_corr = "/data/liul/workarea/XIXI/fit/boost/correctionv8.root";
		fbkgcorr= new TFile(s_corr, "read");
		fpicorr= new TFile("/data/liul/workarea/XIXI/fit/boost/picorrv8.root", "read");
		fpi0corr= new TFile("/data/liul/workarea/XIXI/fit/boost/pion0corrv8.root", "read");



		TString iniyear[4] = {"2009", "2012", "2018", "2019"};
		TString inich[2] = {"xixipm", "xixipp"};
		r1 = new TRandom();

		for(int i = 0; i < 4; i++){
				for(int j = 0; j < 2; j ++){
						angdis[i][j] = new AngDisXiXi();
						angdis[i][j]->InitialInt();
						angdis[i][j]->InitialIntmDIY();
						double pp[2][8];
						pp[0][0] = 0.611;   pp[1][0] = 0.611;
						pp[0][1] = 1.2665;   pp[1][1] = 1.2665;
						pp[0][2] = -0.3722;  pp[1][2] = -0.3722;
						pp[0][3] = -0.0154;    pp[1][3] = -0.0154;
						pp[0][4] = 0.3722;   pp[1][4] = 0.3722;
						pp[0][5] = 0.0154;    pp[1][5] = 0.0154;
						pp[0][6] = 0.6727;   pp[1][6] = 0.7703;
						pp[0][7] = -0.7703;  pp[1][7] = -0.6727;

						angdis[i][j]->SetParameter(pp[j]);
						angdis[i][j]->InitialInt();
						angdis[i][j]->InitialIntmDIY();

						hcorrbkg[i][j] = (TH2D*)fbkgcorr->Get(inich[j] + iniyear[i] + "bkg");
						hpicorr[i][j]  = (TH2D*)fpicorr->Get("picorr" + iniyear[i] + inich[j]);
						hpi0corr[i][j] = (TH2D*)fpi0corr->Get("pion0" + iniyear[i] + inich[j]);

						if(m_seed/1000 == 2){  // for lamdba 2xxx
								r1->SetSeed(m_seed);
								for(int m = 0; m < 10; m++){
										for(int n = 0; n <7; n++){
												double mean = hpicorr[i][j]->GetBinContent(m+1, n+1);
												double meanerr = hpicorr[i][j]->GetBinError(m+1, n+1);
												hpicorr[i][j]->SetBinContent(m+1, n+1, r1->Gaus(mean,meanerr));
										}
								}
						}

						if(m_seed/1000 == 3){  // for lamdba 2xxx
								r1->SetSeed(m_seed);
								for(int m = 0; m < 10; m++){
										for(int n = 0; n <7; n++){
												double mean = hpi0corr[i][j]->GetBinContent(m+1, n+1);
												double meanerr = hpi0corr[i][j]->GetBinError(m+1, n+1);
												hpi0corr[i][j]->SetBinContent(m+1, n+1, r1->Gaus(mean,meanerr));
										}
								}
						}

				}
		}	
}



double rootfile::cpufcnmll(double *pp){
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
		clock_t start,end;
		double loglike[4][2][6];
		for(int i = 0; i < 4; i++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
		}


		double norm[4][2];
		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				norm[iyear][ich] = 0;
				if(!m_norm.CompareTo("phsp")){
						norm[iyear][ich] = angdis[iyear][ich]->CalcToIntegral()/Double_t(D_Sample[i].NN[IndexSample("phsp")]);
				}
				else if(!m_norm.CompareTo("mdiy")){
						norm[iyear][ich] = angdis[iyear][ich]->CalcToIntegralmDIY()/Double_t(D_Sample[i].NN[IndexSample("phsp")]);
				}
		}



		start = clock();

		for(int i = 0; i < (N_Sample +1); i++){
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						int iyear = IndexYear(D_Sample[i].m_year);
						int ich = IndexChannel(D_Sample[i].m_channel);
						int isample = IndexSample(D_Sample[i].m_sample[j]);
						if(!D_Sample[i].m_sample[j].CompareTo("phsp")) continue;
						if(!D_Sample[i].m_sample[j].CompareTo("mdiy")) continue;
						loglike[iyear][ich][isample] = 0;
						for(int evt = 0; evt < D_Sample[i].NN[j]; evt++){
								*(out_amp[iyear][ich][isample] + evt) = angdis[iyear][ich]->Amp(*(*(angdata[iyear][ich][isample]+0)+evt),
												*(*(angdata[iyear][ich][isample]+1)+evt),
												*(*(angdata[iyear][ich][isample]+2)+evt),
												*(*(angdata[iyear][ich][isample]+3)+evt),
												*(*(angdata[iyear][ich][isample]+4)+evt),
												*(*(angdata[iyear][ich][isample]+5)+evt),
												*(*(angdata[iyear][ich][isample]+6)+evt),
												*(*(angdata[iyear][ich][isample]+7)+evt),
												*(*(angdata[iyear][ich][isample]+8)+evt));
								if(*(out_amp[iyear][ich][isample] + evt) <= 0){ cout << "data : " << *(out_amp[iyear][ich][isample] + evt) << endl;  return 0; }
								loglike[iyear][ich][isample] += TMath::Log((*(out_amp[iyear][ich][isample] + evt))/norm[iyear][ich]);
						}
				}
		}
		end = clock();
		double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
		cout << "CPU 3: running kernel " << time3 << " seconds" << endl;
	//	//	exit(1);
		double llk = 0;
		double l1 = 0;

		double scale_factor_1[4][2];
		double scale_factor_2[4][2];
		for(int i = 0; i < 4; i++){
				for(int j = 0; j < 2; j++){
						scale_factor_1[i][j] = 0;
						scale_factor_2[i][j] = 0;
				}
		}




		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								scale_factor_1[iyear][ich] = Double_t(D_Sample[i].NN[IndexSample("data")]);
								scale_factor_2[iyear][ich] = Double_t(D_Sample[i].NN[IndexSample("data")]);
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								if(0 == D_Sample[i].NN[IndexSample("bkg1")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= (D_Sample[i].n_bkg);
										scale_factor_2[iyear][ich] += (D_Sample[i].n_bkg / Double_t(D_Sample[i].NN[IndexSample("bkg1")])) * (D_Sample[i].n_bkg / Double_t(D_Sample[i].NN[IndexSample("bkg1")])) * Double_t(D_Sample[i].NN[IndexSample("bkg1")]);
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								if(0 == D_Sample[i].NN[IndexSample("sideband")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich]-= 0.25 * Double_t(D_Sample[i].NN[IndexSample("sideband")]);   //  		FIXME: It is just a mistake;  Zhipeng uses a tight region of mXi1
										scale_factor_2[iyear][ich]+= 0.25*0.25 * Double_t(D_Sample[i].NN[IndexSample("sideband")]) ;   //
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								if(0 == D_Sample[i].NN[IndexSample("etac")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= D_Sample[i].n_etac;
										scale_factor_2[iyear][ich] += (D_Sample[i].n_etac / Double_t(D_Sample[i].NN[IndexSample("etac")])) * (D_Sample[i].n_etac / Double_t(D_Sample[i].NN[IndexSample("etac")])) * Double_t(D_Sample[i].NN[IndexSample("etac")]);
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								if(0 == D_Sample[i].NN[IndexSample("charge")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= D_Sample[i].n_charge;
										scale_factor_2[iyear][ich] += (D_Sample[i].n_charge / Double_t(D_Sample[i].NN[IndexSample("charge")])) * (D_Sample[i].n_charge / Double_t(D_Sample[i].NN[IndexSample("charge")])) * Double_t(D_Sample[i].NN[IndexSample("charge")]);
								}
						}
				}
		}

		// TODO:with or without the global factor
		/*
		   for(int i = 0; i < 4; i++){
		   for(int j = 0; j < 2; j++){
		   scale_factor_1[i][j] = 1.0;
		   scale_factor_2[i][j] = 1.0;
		   }
		   }
		 */


		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * (-loglike[iyear][ich][j]) ;
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								if(0 == D_Sample[i].NN[IndexSample("bkg1")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_bkg * (loglike[iyear][ich][j]/Double_t(D_Sample[i].NN[IndexSample("bkg1")]));
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								if(0 == D_Sample[i].NN[IndexSample("sideband")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * 0.25 * ( loglike[iyear][ich][j] );		// FIXME  the scale of the sideband  Zhipeng uses a tight region of mXi1
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								if(0 == D_Sample[i].NN[IndexSample("etac")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_etac * ( loglike[iyear][ich][j] / Double_t(D_Sample[i].NN[IndexSample("etac")]));
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								if(0 == D_Sample[i].NN[IndexSample("charge")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_charge * ( loglike[iyear][ich][j] / Double_t(D_Sample[i].NN[IndexSample("charge")]));
								}
						}
				}
		}

		llk =  l1;		

		if(fit_step%100 == 0){
				cout << "likelihood : " << llk << "  ";
				for(int i =0; i < 10; i++){
						cout <<  pp[i] << "  ";
				}
				cout << endl;
		}

		fit_step++;
		return llk;

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
		clock_t start,end;
		double loglike[4][2][6];
		for(int i = 0; i < 4; i++){
				angdis[i][0]->SetParameter(pp1);
				angdis[i][1]->SetParameter(pp2);
		}


		double norm[4][2];
		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				norm[iyear][ich] = 0;
				if(!m_norm.CompareTo("phsp")){
						norm[iyear][ich] = angdis[iyear][ich]->CalcToIntegral()/Double_t(D_Sample[i].NN[IndexSample("phsp")]);
				}
				else if(!m_norm.CompareTo("mdiy")){
						norm[iyear][ich] = angdis[iyear][ich]->CalcToIntegralmDIY()/Double_t(D_Sample[i].NN[IndexSample("phsp")]);
				}
		}



		start = clock();
		clock_t start1,end1;
		clock_t start2,end2;
	//	clock_t start3,end3;

		for(int i = 0; i < (N_Sample +1); i++){
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						int iyear = IndexYear(D_Sample[i].m_year);
						int ich = IndexChannel(D_Sample[i].m_channel);
						int isample = IndexSample(D_Sample[i].m_sample[j]);
						if(!D_Sample[i].m_sample[j].CompareTo("phsp")) continue;
						if(!D_Sample[i].m_sample[j].CompareTo("mdiy")) continue;
						start1 = clock();
						gpu_Amp <<< (D_Sample[i].NN[j] * MATRIX_SIZE + MATRIX_SIZE * THREADS_PER_BLOCK ) / (MATRIX_SIZE * THREADS_PER_BLOCK), MATRIX_SIZE * THREADS_PER_BLOCK >>> ( 
										*(gpu_angdata[iyear][ich][isample] + 0), 
										*(gpu_angdata[iyear][ich][isample] + 1), 
										*(gpu_angdata[iyear][ich][isample] + 2), 
										*(gpu_angdata[iyear][ich][isample] + 3), 
										*(gpu_angdata[iyear][ich][isample] + 4), 
										*(gpu_angdata[iyear][ich][isample] + 5), 
										*(gpu_angdata[iyear][ich][isample] + 6), 
										*(gpu_angdata[iyear][ich][isample] + 7), 
										*(gpu_angdata[iyear][ich][isample] + 8),
										gpu_amp[iyear][ich][isample],
										(D_Sample[i].NN[j] + THREADS_PER_BLOCK)*80, 
										aa_para, ich, gpu_Matrix[iyear][ich][isample], norm[iyear][ich]);
						end1 = clock();
						cudaDeviceSynchronize(); // wait until prior kernel is finished
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 004!" << endl;
								exit(1);
						}
						start2 = clock();
						loglike[iyear][ich][isample] = gpu_float_sum_reduce(gpu_amp[iyear][ich][isample], D_Sample[i].NN[j]);
						end2 = clock();
	//					double time13 = ((double)(end1-start1))/CLOCKS_PER_SEC;
	//					cout << "GPU 13: running kernel " << time13 << " seconds" << endl;
	//					double time23 = ((double)(end2-start2))/CLOCKS_PER_SEC;
	//					cout << "GPU 23: running kernel " << time23 << " seconds" << endl;
/*
						int mat_size = (D_Sample[i].NN[j] + THREADS_PER_BLOCK) *sizeof(double);
						cudaMemcpy( out_amp[iyear][ich][isample], gpu_amp[iyear][ich][isample], mat_size, cudaMemcpyDeviceToHost );
						cudaStatus = cudaGetLastError();
						if (cudaStatus != cudaSuccess){
								cerr << "failure to call cuda kernel 002!" << endl;
								exit(1);
						}
						loglike[iyear][ich][isample] = 0;
						if(!D_Sample[i].m_sample[j].CompareTo("phsp")){
								for(int evt = 0; evt < D_Sample[i].NN[j]; evt++){
										if(*(out_amp[iyear][ich][isample] + evt) <= 0){ cout << "data : " << *(out_amp[iyear][ich][isample] + evt) << endl;  return 0; }
										loglike[iyear][ich][isample] += *(out_amp[iyear][ich][isample] + evt);
								}

						}
						else{
						start3 = clock();
								for(int evt = 0; evt < D_Sample[i].NN[j]; evt++){
										//	cout << "host C munu 0: " << 	*(host_eval + evt) << endl;
										if(*(out_amp[iyear][ich][isample] + evt) <= 0){ cout << "data : " << *(out_amp[iyear][ich][isample] + evt) << endl;  return 0; }
										loglike[iyear][ich][isample] += TMath::Log((*(out_amp[iyear][ich][isample] + evt))/norm[iyear][ich]);
								}

						end3 = clock();
						}
						double time13 = ((double)(end1-start1))/CLOCKS_PER_SEC;
						cout << "GPU 13: running kernel " << time13 << " seconds" << endl;
						double time23 = ((double)(end2-start2))/CLOCKS_PER_SEC;
						cout << "GPU 23: running kernel " << time23 << " seconds" << endl;
						double time33 = ((double)(end3-start3))/CLOCKS_PER_SEC;
						cout << "GPU 33: running kernel " << time33 << " seconds" << endl;
						*/
				}
		}
		end = clock();
		double time3 = ((double)(end-start))/CLOCKS_PER_SEC;
		cout << "GPU 3: running kernel " << time3 << " seconds" << endl;
		//	exit(1);
		double llk = 0;
		double l1 = 0;

		double scale_factor_1[4][2];
		double scale_factor_2[4][2];
		for(int i = 0; i < 4; i++){
				for(int j = 0; j < 2; j++){
						scale_factor_1[i][j] = 0;
						scale_factor_2[i][j] = 0;
				}
		}

		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								scale_factor_1[iyear][ich] = Double_t(D_Sample[i].NN[IndexSample("data")]);
								scale_factor_2[iyear][ich] = Double_t(D_Sample[i].NN[IndexSample("data")]);
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								if(0 == D_Sample[i].NN[IndexSample("bkg1")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= (D_Sample[i].n_bkg);
										scale_factor_2[iyear][ich] += (D_Sample[i].n_bkg / Double_t(D_Sample[i].NN[IndexSample("bkg1")])) * (D_Sample[i].n_bkg / Double_t(D_Sample[i].NN[IndexSample("bkg1")])) * Double_t(D_Sample[i].NN[IndexSample("bkg1")]);
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								if(0 == D_Sample[i].NN[IndexSample("sideband")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich]-= 0.25 * Double_t(D_Sample[i].NN[IndexSample("sideband")]);   //  		FIXME: It is just a mistake;  Zhipeng uses a tight region of mXi1
										scale_factor_2[iyear][ich]+= 0.25*0.25 * Double_t(D_Sample[i].NN[IndexSample("sideband")]) ;   //
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								if(0 == D_Sample[i].NN[IndexSample("etac")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= D_Sample[i].n_etac;
										scale_factor_2[iyear][ich] += (D_Sample[i].n_etac / Double_t(D_Sample[i].NN[IndexSample("etac")])) * (D_Sample[i].n_etac / Double_t(D_Sample[i].NN[IndexSample("etac")])) * Double_t(D_Sample[i].NN[IndexSample("etac")]);
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								if(0 == D_Sample[i].NN[IndexSample("charge")]){
										scale_factor_1[iyear][ich] -= 0;
										scale_factor_2[iyear][ich] += 0;
								}
								else{
										scale_factor_1[iyear][ich] -= D_Sample[i].n_charge;
										scale_factor_2[iyear][ich] += (D_Sample[i].n_charge / Double_t(D_Sample[i].NN[IndexSample("charge")])) * (D_Sample[i].n_charge / Double_t(D_Sample[i].NN[IndexSample("charge")])) * Double_t(D_Sample[i].NN[IndexSample("charge")]);
								}
						}
				}
		}

		// TODO:with or without the global factor
		/*
		   for(int i = 0; i < 4; i++){
		   for(int j = 0; j < 2; j++){
		   scale_factor_1[i][j] = 1.0;
		   scale_factor_2[i][j] = 1.0;
		   }
		   }
		 */


		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * (-loglike[iyear][ich][j]) ;
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								if(0 == D_Sample[i].NN[IndexSample("bkg1")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_bkg * (loglike[iyear][ich][j]/Double_t(D_Sample[i].NN[IndexSample("bkg1")]));
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								if(0 == D_Sample[i].NN[IndexSample("sideband")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * 0.25 * ( loglike[iyear][ich][j] );		// FIXME  the scale of the sideband  Zhipeng uses a tight region of mXi1
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								if(0 == D_Sample[i].NN[IndexSample("etac")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_etac * ( loglike[iyear][ich][j] / Double_t(D_Sample[i].NN[IndexSample("etac")]));
								}
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								if(0 == D_Sample[i].NN[IndexSample("charge")]){
										l1 += 0.;
								}
								else{
										l1 += ( scale_factor_1[iyear][ich] / scale_factor_2[iyear][ich]) * D_Sample[i].n_charge * ( loglike[iyear][ich][j] / Double_t(D_Sample[i].NN[IndexSample("charge")]));
								}
						}
				}
		}

		llk =  l1;		

		if(fit_step%100 == 0){
				cout << "likelihood : " << llk << "  ";
				for(int i =0; i < 10; i++){
						cout <<  pp[i] << "  ";
				}
				cout << endl;
		}

		fit_step++;
		return llk;
}


double rootfile::setBKGMassFit( double nn, TString year, TString channel, TString sample){
		double Njpis[4] = {224000000., 1087000000, 4387000000, 4387000000};
		double NMCetac[4] = {45000, 195000, 750000, 750000};
		double NMCcharge[4] = {3000000, 18000000, 55500000, 60000000};
		int iyear = IndexYear(year);
		int ich = IndexChannel(channel);
		if(!sample.CompareTo("etac")){
				return Njpis[iyear] * 0.017*0.0009*0.99524*0.99524*0.639*0.358* nn / NMCetac[iyear];
		}
		else if(!sample.CompareTo("charge")){
				return Njpis[iyear] * 0.0011*0.99524*0.99524*0.639*0.638* nn / NMCcharge[iyear];
		}
}



void rootfile::setBKG(){
		double Njpis[4] = {224000000., 1087000000, 4387000000, 4387000000};
		double NMCetac[4] = {45000, 195000, 750000, 750000};
		double NMCcharge[4] = {3000000, 18000000, 55500000, 60000000};
		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								D_Sample[i].n_etac = Njpis[iyear] * 0.017*0.0009*0.99524*0.99524*0.639*0.358* D_Sample[i].NN[j] / NMCetac[iyear];
								D_Sample[i].n_etacerr = D_Sample[i].n_etac * sqrt((0.4/1.7)*(0.4/1.7) 
												+ (2.6/9.0)*(2.6/9.0) + (0.5/63.9)*(0.5/63.9) + 
												(0.5/35.8)*(0.5/35.8));
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								D_Sample[i].n_charge = Njpis[iyear] * 0.0011*0.99524*0.99524*0.639*0.638* D_Sample[i].NN[j] / NMCcharge[iyear];
								D_Sample[i].n_chargeerr = D_Sample[i].n_charge * sqrt(
												(0.8/9.7)*(0.8/9.7) + (0.5/63.9)*(0.5/63.9) + 
												(0.5/63.9)*(0.5/63.9));

								cout << "background charge : " << D_Sample[i].n_charge << "	" << D_Sample[i].n_chargeerr << endl;
						}
				}
		}
}

Int_t rootfile::CorrFactor(TString year, TString channel, TString sample){
		Double_t factor = 0;
		int iyear = IndexYear(year);
		int ich = IndexChannel(channel);
		if(!sample.CompareTo("phsp")){
				int i_pion1_1 = floor((m_pion1_1_cos + 1) / 0.2);
				int j_pion1_1 = floor(m_pion1_1_pt / 0.05);
				int i_pion2_1 = floor((m_pion2_1_cos + 1) / 0.2);
				int j_pion2_1 = floor(m_pion2_1_pt / 0.05);
				int i_pion0 = floor((m_pion0_cos + 1) / 0.2);
				int j_pion0 = floor(m_pion0_rho / 0.05);
				if(j_pion0 <= 2) j_pion0 = 2;
				if(j_pion0 >= 4) j_pion0 = 4;
				if(j_pion2_1 >=5) j_pion2_1 =5;
				if(j_pion1_1 >=5) j_pion1_1 =5;
				int ich2 = ich ? 0 : 1;
				factor = hpi0corr[iyear][ich]->GetBinContent(i_pion0+1, j_pion0+1) * hpicorr[iyear][ich2]->GetBinContent(i_pion1_1+1, j_pion1_1+1) * hpicorr[iyear][ich]->GetBinContent(i_pion2_1+1, j_pion2_1+1);

		}
		else if(!sample.CompareTo("bkg1")){
				int i = floor((m_lmd_cos + 1) / 0.2);
				int j = floor((m_lmd_p - 0.5)/ 0.05);

				int i_pion1_1 = floor((m_pion1_1_cos + 1) / 0.2);
				int j_pion1_1 = floor(m_pion1_1_pt / 0.05);
				int i_pion2_1 = floor((m_pion2_1_cos + 1) / 0.2);
				int j_pion2_1 = floor(m_pion2_1_pt / 0.05);
				if(j_pion2_1 >=5) j_pion2_1 =5;
				if(j_pion1_1 >=5) j_pion1_1 =5;

				int ich2 = ich ? 0 : 1;

				if(j >=0 && j < 7){
						factor = (hcorrbkg[iyear][ich]->GetBinContent(i+1, j+1)) * hpicorr[iyear][ich2]->GetBinContent(i_pion1_1+1, j_pion1_1+1) * hpicorr[iyear][ich]->GetBinContent(i_pion2_1+1, j_pion2_1+1);
				}
				else{
						factor  = 1.0;
				}
		}

		double rdm = r1->Rndm();
		if(factor < 1.0 && rdm > factor){
				return 0;
		}
		else if(factor < 1.0 && rdm < factor){
				return 1;
		}
		else if(factor > 1.0 && rdm < factor - 1.0 && factor < 2.0){
				return 2;
		}
		else if(factor > 1.0 && rdm > factor - 1.0 && factor < 2.0){
				return 1;
		}
		else if(factor > 2.0 && rdm < factor - 2.0 && factor < 3.0){
				return 3;
		}
		else if(factor > 2.0 && rdm > factor - 2.0 && factor < 3.0){
				return 2;
		}
		else if(factor > 3.0){  	// almost no event has a factor larger than 3.
				return 3;
		}
		else 
				return 1.0;
}
Double_t rootfile::calculate_int(Double_t par0, Double_t par1, Double_t intlow, Double_t intup){
		double cosup = intup;
		double coslow = intlow;
		double intfuncup = par0*(cosup+par1*pow(cosup,3)/3);
		double intfunclow  = par0*(coslow+par1*pow(coslow,3)/3);
		double integral = intfuncup - intfunclow;
		return integral;
}


void rootfile::readDataIO(int index, int jndex){
		TString path = "/data/liul/workarea/XIXI/fit/boost";
		TString year = D_Sample[index].m_year;
		TString channel = D_Sample[index].m_channel;
		TString version = D_Sample[index].m_version;
		TString sample = D_Sample[index].m_sample[jndex];
		TString insample = sample;
		if(!sample.CompareTo("data")){
				insample = "mdiy";
		}
		if(!sample.CompareTo("phsp")){
				if(!m_norm.CompareTo("mdiy")){
						insample = "mdiy";
				}
				else if(!m_norm.CompareTo("phsp")){
						insample = "phsp";
				}
		}
		if(!sample.CompareTo("bkg1")){
				insample = "mdiy";
		}
		if(!sample.CompareTo("sideband")){
				insample = "data";
		}
		TString infile  = path + "/" + year + "/" + channel + "/" + version + "/" + insample + "/boost.root";
		cout << infile << " ==> " << sample << endl;
		//#define read_data
#ifndef read_data
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



		TFile *f1 = new TFile(infile, "read");
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
		t1->SetBranchAddress("runNo", &m_runNo);
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
		if(sample.CompareTo("etac")){
				t1->SetBranchAddress("lmd_p", &m_lmd_p);
				t1->SetBranchAddress("lmd_cos", &m_lmd_cos);
		}

		int nn = 0;
		int NEvt = t1->GetEntries();
		int low = 0;
		int high = NEvt;
		int LOW = -1;
		int HIGH = -1;
		r1->SetSeed(3251);
		if(!sample.CompareTo("data")){
				low = i_trial * NEvt/30;
				high = (i_trial + 1) * NEvt/30;
		}
		else if(!sample.CompareTo("phsp")){
				LOW =  i_trial * NEvt/30;
				HIGH = (i_trial + 1) * NEvt/30;
		}
		else if(!sample.CompareTo("bkg1")){
				LOW =  i_trial * NEvt/30;
				HIGH = (i_trial + 1) * NEvt/30;
		}

		for(int i = low; i <  high; i++){
				t1->GetEntry(i);
				if(i > LOW && i < HIGH) continue;
				if(!Selection(year, channel, sample, index)) continue;
				if(!sample.CompareTo("phsp")){
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
						if(!m_norm.CompareTo("phsp")){
								angdis[IndexYear(year)][IndexChannel(channel)]->AddToIntegral(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
						}
						else if(!m_norm.CompareTo("mdiy")){
								angdis[IndexYear(year)][IndexChannel(channel)]->AddToIntegralmDIY(the, Lthe, Lphi, Lbthe, Lbphi, pthe, pphi, apthe, apphi);
						}

				}
				else{
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
		}
		D_Sample[index].NN.push_back(nn);

		int iyear = IndexYear(year);
		int ich = IndexChannel(channel);
		int isample = IndexSample(sample);
		cout << "Hello " << iyear << ich << isample << endl;

		for(int i = 0; i < nn; i++){
				*(*(angdata[iyear][ich][isample]+0)+i) = gD1Xithe[i];
				*(*(angdata[iyear][ich][isample]+1)+i) = gD1Lthe[i];
				*(*(angdata[iyear][ich][isample]+2)+i) = gD1Lphi[i];
				*(*(angdata[iyear][ich][isample]+3)+i) = gD1Lbthe[i];
				*(*(angdata[iyear][ich][isample]+4)+i) = gD1Lbphi[i];
				*(*(angdata[iyear][ich][isample]+5)+i) = gD1pthe[i];
				*(*(angdata[iyear][ich][isample]+6)+i) = gD1pphi[i];
				*(*(angdata[iyear][ich][isample]+7)+i) = gD1apthe[i];
				*(*(angdata[iyear][ich][isample]+8)+i) = gD1apphi[i];
		}


		cout << nn << endl;
		f1->Close();
#endif
}

void rootfile::Print(){

		cout << "iyear	channel	signal		bkg		sideband	etac	charge" << endl;
		for(int i = 0; i < (N_Sample +1); i++){
				int iyear = IndexYear(D_Sample[i].m_year);
				int ich = IndexChannel(D_Sample[i].m_channel);
				cout << std::left << setw(4) << D_Sample[i].m_year << "	" << D_Sample[i].m_channel << "	";
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								cout << D_Sample[i].n_signal << "±" << D_Sample[i].n_signalerr << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								cout << D_Sample[i].n_bkg << "±" << D_Sample[i].n_bkgerr << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								cout << D_Sample[i].NN[IndexSample("sideband")]/4 << "±" << 
										sqrt(D_Sample[i].NN[IndexSample("sideband")]/4)<< "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								cout << D_Sample[i].n_etac << "±" << D_Sample[i].n_etacerr << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								cout << D_Sample[i].n_charge << "±" << D_Sample[i].n_chargeerr << "	";
						}
				}
				cout << endl;
				cout << std::left << setw(4) << D_Sample[i].m_year << "	" << D_Sample[i].m_channel << "	";
				for(int j = 0; j < D_Sample[i].m_sample.size(); j++){
						if(!D_Sample[i].m_sample[j].CompareTo("data")){
								cout << D_Sample[i].NN[IndexSample("data")] << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("bkg1")){
								cout << D_Sample[i].NN[IndexSample("bkg1")] << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("sideband")){
								cout << D_Sample[i].NN[IndexSample("sideband")] << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("etac")){
								cout << D_Sample[i].NN[IndexSample("etac")] << "	";
						}
						else if(!D_Sample[i].m_sample[j].CompareTo("charge")){
								cout << D_Sample[i].NN[IndexSample("charge")] << "	";
						}
				}
				cout << endl;
		}
}

void rootfile::SetBKGSysTest(vector<TString> vstr){
		for(int i = 0; i < vstr.size(); i++){
				if(!vstr[i].CompareTo("bkg1p1")){ // plus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_bkg +=  D_Sample[j].n_bkgerr;
						}
				}
				else if(!vstr[i].CompareTo("bkg1m1")){ // minus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_bkg -=  D_Sample[j].n_bkgerr;
						}
				}
				else if(!vstr[i].CompareTo("etac1p1")){ // plus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_etac +=  D_Sample[j].n_etacerr;
						}
				}
				else if(!vstr[i].CompareTo("etac1m1")){ // minus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_etac -=  D_Sample[j].n_etacerr;
						}
				}
				else if(!vstr[i].CompareTo("charge1p1")){ // plus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_charge +=  D_Sample[j].n_chargeerr;
								cout << "Backgroudn charge : " << D_Sample[j].n_charge << endl;
						}
				}
				else if(!vstr[i].CompareTo("charge1m1")){ // minus 1 sigma
						for(int j = 0; j < (N_Sample +1); j++){
								D_Sample[j].n_charge -=  D_Sample[j].n_chargeerr;
						}
				}
		}
}
