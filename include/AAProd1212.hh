#ifndef AA_Prod1212_H_
#define AA_Prod1212_H_
#include "TMath.h"
class AAProd1212 {
		public:
				Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
				Double_t C1(Int_t i, Int_t j, Int_t k) const {return thC1[i][j][k];}
				Double_t C2(Int_t i, Int_t j, Int_t k) const {return thC2[i][j][k];}
				AAProd1212(double alpha_Jpsi, double delta_phi, double Jpsi_the, int flg) {
						// decay 1-- -> 1/2 + 1/2
						if(flg == 0){
								Double_t v  = sqrt(1-alpha_Jpsi*alpha_Jpsi);
								Double_t st = sin(Jpsi_the);
								Double_t ct = cos(Jpsi_the);
								Double_t cp = cos(delta_phi);
								Double_t sp = sin(delta_phi);
								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												Double_t c=0;
												if(k==0&&l==0)c=(1+alpha_Jpsi*ct*ct);
												if(k==0&&l==2)c=v*st*ct*sp;
												if(k==1&&l==1)c=st*st;
												if(k==1&&l==3)c=v*st*ct*cp;
												if(k==2&&l==0)c=-v*st*ct*sp;
												if(k==2&&l==2)c=alpha_Jpsi*st*st;
												if(k==3&&l==1)c=-v*st*ct*cp;
												if(k==3&&l==3)c=-alpha_Jpsi - ct*ct;
												thC[k][l]=c;
										}
								}
						}

						else if(flg == 1){
								// decay 1-- -> 1/2 + 1/2
								Double_t v  = 1;
								Double_t st = sin(Jpsi_the);
								Double_t ct = cos(Jpsi_the);
								Double_t cp = 1;
								Double_t sp = 1;
								Double_t alpha_Jpsi = 1;
								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												for(int i=0;i<2;i++){
														Double_t c=0;
														if(i==0&&k==0&&l==0)c=1;
														if(i==1&&k==0&&l==0)c=(alpha_Jpsi*ct*ct);
														if(i==0&&k==0&&l==2)c=v*st*ct*sp;
														if(i==1&&k==0&&l==2)c=0;
														if(i==0&&k==1&&l==1)c=st*st;
														if(i==1&&k==1&&l==1)c=0;
														if(i==0&&k==1&&l==3)c=v*st*ct*cp;
														if(i==1&&k==1&&l==3)c=0;
														if(i==0&&k==2&&l==0)c=-v*st*ct*sp;
														if(i==1&&k==2&&l==0)c=0;
														if(i==0&&k==2&&l==2)c=alpha_Jpsi*st*st;
														if(i==1&&k==2&&l==2)c=0;
														if(i==0&&k==3&&l==1)c=-v*st*ct*cp;
														if(i==1&&k==3&&l==1)c=0;
														if(i==0&&k==3&&l==3)c=-alpha_Jpsi ;
														if(i==1&&k==3&&l==3)c=-ct*ct;
														thC1[k][l][i]=c;
												}
										}
								}

						}
						else if(flg == 2){

								double v  = sqrt(1-alpha_Jpsi*alpha_Jpsi);
								double st = 1;
								double ct = 1;
								double cp = cos(delta_phi);
								double sp = sin(delta_phi);

								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												for(int i=0;i<2;i++){
														Double_t c=1;
														if(i==0&&k==0&&l==0)c=(1);
														if(i==1&&k==0&&l==0)c=(alpha_Jpsi*ct*ct);
														if(i==0&&k==0&&l==2)c=v*st*ct*sp;
														if(i==1&&k==0&&l==2)c=1;
														if(i==0&&k==1&&l==1)c=st*st;
														if(i==1&&k==1&&l==1)c=1;
														if(i==0&&k==1&&l==3)c=v*st*ct*cp;
														if(i==1&&k==1&&l==3)c=1;
														if(i==0&&k==2&&l==0)c=v*st*ct*sp;
														if(i==1&&k==2&&l==0)c=1;
														if(i==0&&k==2&&l==2)c=alpha_Jpsi*st*st;
														if(i==1&&k==2&&l==2)c=1;
														if(i==0&&k==3&&l==1)c=v*st*ct*cp;
														if(i==1&&k==3&&l==1)c=1;
														if(i==0&&k==3&&l==3)c=alpha_Jpsi ;
														if(i==1&&k==3&&l==3)c=1;
														thC2[k][l][i]=c;
												}
										}
								}
						}

				}

				~AAProd1212() {}
		private:
				double thC[4][4];
				double thC1[4][4][2];
				double thC2[4][4][2];
};
#endif /* AA_Prod1212_H_ */
