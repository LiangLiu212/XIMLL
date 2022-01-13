#ifndef AA_Prod1212_H_
#define AA_Prod1212_H_
#include "TMath.h"
class AAProd1212 {
		public:
				Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
				AAProd1212(double alpha_Jpsi, double delta_phi, double Jpsi_the) {
						// decay 1-- -> 1/2 + 1/2
						Double_t v  = TMath::Sqrt(1-alpha_Jpsi*alpha_Jpsi);
						Double_t st = TMath::Sin(Jpsi_the);
						Double_t ct = TMath::Cos(Jpsi_the);
						Double_t cp = TMath::Cos(delta_phi);
						Double_t sp = TMath::Sin(delta_phi);
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
				~AAProd1212() {}
		private:
				double thC[4][4];
};
#endif /* AA_Prod1212_H_ */
