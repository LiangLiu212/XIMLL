#ifndef MN_AngDisXiXi_H_
#define MN_AngDisXiXi_H_
#include "AAProd1212.hh"
#include "AADecay12.hh"
#include "TString.h"
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;
class AngDisXiXi {
		//e+e- -> (Ss->L pi-) (S->nb pi+) (L->p pi-)
		public:


				AngDisXiXi(){
						isPHSPinte = false;
						datamass = new double [1000000];
						mdiymass = new double [1000000];
				}   

				~AngDisXiXi() {}
				 double Amp(double theta,
								double LThe,double LPhi,double LbThe,
								double LbPhi,double pThe,double pPhi, 
								double apThe, double apPhi) 
				{
						AAProd1212 prod(p_Jpsi_alpha, p_Jpsi_phi, theta, 0);
						AADecay12 Xi(p_Xi_alpha, p_Xi_phi, LThe, LPhi, 0);
						AADecay12 Xibar(p_Xibar_alpha, p_Xibar_phi, LbThe, LbPhi, 0);
						AADecay12 L(p_L_alpha, 0, pThe, pPhi, 0);
						AADecay12 Lbar(p_Lbar_alpha, 0, apThe, apPhi, 0);
						double tep=0;
						for(int mu=0; mu<4;mu++){// Xi loop
								for(int nu=0;nu<4;nu++){// Xibar loop
										for(int k=0;k<4;k++){
												for(int j=0;j<4;j++){
														tep+=prod.C(mu, nu)*
																Xi.A(mu,k)*Xibar.A(nu, j)*
																L.A(k,0)*Lbar.A(j, 0);
												}
										}
								}
						}
						return tep;
				}

				void InitialInt(){
						for (int i =0 ; i < 2048; i ++){
								 atoInt[i] = 0;
						}
				}

				void InitialIntmDIY(){
						for (int i =0 ; i < 2048; i ++){
								 atoIntmDIY[i] = 0;
						}
				}

				bool isphspintegral(){return isPHSPinte;};

				void AddToIntegralmDIY( double theta, double LThe,double LPhi,double LbThe, double LbPhi,double pThe,double pPhi, double apThe, double apPhi){
						isPHSPinte =true;

						double amp =  Amp(theta, LThe, LPhi, LbThe, LbPhi, pThe, pPhi, apThe, apPhi);

						AAProd1212 prod(1, 1, theta, 1);
						AADecay12 Xi(1, 1, LThe, LPhi, 1);
						AADecay12 Xibar(1, 1, LbThe, LbPhi, 1);
						AADecay12 L(1, 0, pThe, pPhi, 1);
						AADecay12 Lbar(1, 0, apThe, apPhi, 1);

						for(int mu=0; mu<4;mu++){// Xi loop
								for(int nu=0;nu<4;nu++){// Xibar loop
										for(int k=0;k<4;k++){
												for(int j=0;j<4;j++){
														for(int c1 = 0; c1< 2; c1++){
																for(int c2 = 0; c2 < 2; c2++){
																		for(int c3 = 0; c3 < 2; c3++){
																				Int_t idx = mu*512 + nu*128 + k*32 +j*8 + c1*4 + c2*2 + c3;
																				atoIntmDIY[idx]+=prod.C1(mu, nu, c1)*
																						Xi.Ab(mu,k, c2)*Xibar.Ab(nu, j, c3)*
																						L.Ab(k,0, 0)*Lbar.Ab(j, 0, 0)/amp;
																		}
																}
														}
												}
										}
								}
						}
				}


				void AddToIntegral( double theta, double LThe,double LPhi,double LbThe, double LbPhi,double pThe,double pPhi, double apThe, double apPhi  ){
						isPHSPinte =true;


						AAProd1212 prod(1, 1, theta, 1);
						AADecay12 Xi(1, 1, LThe, LPhi, 1);
						AADecay12 Xibar(1, 1, LbThe, LbPhi, 1);
						AADecay12 L(1, 0, pThe, pPhi, 1);
						AADecay12 Lbar(1, 0, apThe, apPhi, 1);

						for(int mu=0; mu<4;mu++){// Xi loop
								for(int nu=0;nu<4;nu++){// Xibar loop
										for(int k=0;k<4;k++){
												for(int j=0;j<4;j++){
														for(int c1 = 0; c1< 2; c1++){
																for(int c2 = 0; c2 < 2; c2++){
																		for(int c3 = 0; c3 < 2; c3++){
																				Int_t idx = mu*512 + nu*128 + k*32 +j*8 + c1*4 + c2*2 + c3;
																				atoInt[idx]+=prod.C1(mu, nu, c1)*
																						Xi.Ab(mu,k, c2)*Xibar.Ab(nu, j, c3)*
																						L.Ab(k,0, 0)*Lbar.Ab(j, 0, 0);
																		}
																}
														}
												}
										}
								}
						}
				}

				void PrintInt(){
						for (int i =0 ; i < 2048; i ++){
								if(atoInt[i] != 0)
										cout << atoInt[i] << endl;
						}
				}

				void SaveInt(const TString outfile, const int NN){
						ofstream oInt(outfile, ios::out | ios::binary);
						if(!oInt){
								cerr << "open error! => " << outfile << endl;
								exit(2);
						}
						oInt.write((char*)&NN, sizeof(int));
						oInt.write((char*)atoInt, sizeof(double) * 2048);
						oInt.close();
				}
				void SaveIntmDIY(const TString outfile, const int NN){
						ofstream oInt(outfile, ios::out | ios::binary);
						if(!oInt){
								cerr << "open error! => " << outfile << endl;
								exit(2);
						}
						oInt.write((char*)&NN, sizeof(int));
						oInt.write((char*)atoIntmDIY, sizeof(double) * 2048);
						oInt.close();
				}

				int ReadInt(const TString infile){
						ifstream oint(infile, ios::in | ios::binary);
						if(!oint){
								cerr << "open error! => " << infile << endl;
								exit(2);
						}
						int NN = 0;
						oint.read((char*)&NN, sizeof(int));
					//	cout << "NN : " << NN << endl;
						oint.read((char*)atoInt, sizeof(double) * 2048);
						oint.close();
						return NN;
				}
				int ReadIntmDIY(const TString infile){
						ifstream oint(infile, ios::in | ios::binary);
						if(!oint){
								cerr << "open error! => " << infile << endl;
								exit(2);
						}
						int NN = 0;
						oint.read((char*)&NN, sizeof(int));
					//	cout << "NN : " << NN << endl;
						oint.read((char*)atoIntmDIY, sizeof(double) * 2048);
					//	for(int i  = 0; i < 2048; i++){
					//			cout << atoIntmDIY[i] << endl;
					//	}
						oint.close();
						return NN;
				}

				double CalcToIntegral(){

						AAProd1212 prod(p_Jpsi_alpha, p_Jpsi_phi,1, 2);
						AADecay12 Xi(p_Xi_alpha, p_Xi_phi, 1, 1, 2);
						AADecay12 Xibar(p_Xibar_alpha, p_Xibar_phi, 1, 1, 2);
						AADecay12 L(p_L_alpha, 0, 1, 1, 2);
						AADecay12 Lbar(p_Lbar_alpha, 0, 1, 1, 2);
						double tmp = 0;

						for(int mu=0; mu<4;mu++){// Xi loop
								for(int nu=0;nu<4;nu++){// Xibar loop
										for(int k=0;k<4;k++){
												for(int j=0;j<4;j++){
														for(int c1 = 0; c1< 2; c1++){
																for(int c2 = 0; c2 < 2; c2++){
																		for(int c3 = 0; c3 < 2; c3++){
																				Int_t idx = mu*512 + nu*128 + k*32 +j*8 + c1*4 + c2*2 + c3;
																				tmp+=prod.C2(mu, nu, c1)*
																						Xi.Ac(mu,k, c2)*Xibar.Ac(nu, j, c3)*
																						L.Ac(k,0, 0)*Lbar.Ac(j, 0, 0)*atoInt[idx];
																		}
																}
														}
												}
										}
								}
						}
						return tmp;
				}

				double CalcToIntegralmDIY(){

						AAProd1212 prod(p_Jpsi_alpha, p_Jpsi_phi,1, 2);
						AADecay12 Xi(p_Xi_alpha, p_Xi_phi, 1, 1, 2);
						AADecay12 Xibar(p_Xibar_alpha, p_Xibar_phi, 1, 1, 2);
						AADecay12 L(p_L_alpha, 0, 1, 1, 2);
						AADecay12 Lbar(p_Lbar_alpha, 0, 1, 1, 2);
						double tmp = 0;

						for(int mu=0; mu<4;mu++){// Xi loop
								for(int nu=0;nu<4;nu++){// Xibar loop
										for(int k=0;k<4;k++){
												for(int j=0;j<4;j++){
														for(int c1 = 0; c1< 2; c1++){
																for(int c2 = 0; c2 < 2; c2++){
																		for(int c3 = 0; c3 < 2; c3++){
																				Int_t idx = mu*512 + nu*128 + k*32 +j*8 + c1*4 + c2*2 + c3;
																				tmp+=prod.C2(mu, nu, c1)*
																						Xi.Ac(mu,k, c2)*Xibar.Ac(nu, j, c3)*
																						L.Ac(k,0, 0)*Lbar.Ac(j, 0, 0)*atoIntmDIY[idx];
																		}
																}
														}
												}
										}
								}
						}
						return tmp;
				}



				void SetParameter(double *pp){
						p_Jpsi_alpha = pp[0];
						p_Jpsi_phi = pp[1];
						p_Xi_alpha= pp[2];
						p_Xi_phi =pp[3];
						p_Xibar_alpha=pp[4];
						p_Xibar_phi=pp[5];
						p_L_alpha=pp[6];
						p_Lbar_alpha=pp[7];
				}

				void setDataMass(int i, double mass){
						datamass[i] = mass;
						Ndata = i;
				}
				void setMCMass(int i, double mass){
						mdiymass[i] = mass;
						Nmdiy = i;
				}
				double DataMass(int n){
						return datamass[n];
				}
				double MCMass(int n){
						return mdiymass[n];
				}
				int getNdata(){
						return Ndata+1;
				}
				int getNmc(){
						return Nmdiy+1;
				}

		private:
				bool isPHSPinte;
				double  p_Jpsi_alpha;
				double 	p_Jpsi_phi;
				double 	p_Xi_alpha;
				double  p_Xi_phi;
				double 	p_Xibar_alpha;
				double  p_Xibar_phi;
				double  p_L_alpha;
				double  p_Lbar_alpha;

				double atoInt[2048];
				double atoIntmDIY[2048];
				double *datamass;
				double *mdiymass;
				int Ndata;
				int Nmdiy;
};
#endif // MN_AngDisSsS_H_
