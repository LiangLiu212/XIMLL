#ifndef AA_Decay12_H_
#define AA_Decay12_H_
#include "TMath.h"
class AADecay12 {
		public:
				AADecay12(const double alpha, const double xiphi, const double theta, const double phi, const int flg) {
						
						double cos_theta = cos(theta);
						double sin_theta = sin(theta);
						double cos_phi = cos(phi);
						double sin_phi = sin(phi);

						// decay 1/2 -> 1/2 + 0
						if(flg == 0){
								double gamma = sqrt(1-alpha*alpha)*cos(xiphi);
								double beta = sqrt(1-alpha*alpha)*sin(xiphi);
								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												Double_t c=0;

												if(k==0&& l == 0) c=1;

												if(k==0&& l == 3) c=alpha;

												if(k==1&& l == 0) c=alpha*cos_phi*sin_theta;

												if(k==1&& l == 1) c=gamma*cos_phi*cos_theta - beta*sin_phi;

												if(k==1&& l == 2) c=-(beta*cos_phi*cos_theta) - gamma*sin_phi;

												if(k==1&& l == 3) c=cos_phi*sin_theta;

												if(k==2&& l == 0) c=alpha*sin_phi*sin_theta;

												if(k==2&& l == 1) c=beta*cos_phi + gamma*cos_theta*sin_phi;

												if(k==2&& l == 2) c=gamma*cos_phi - beta*cos_theta*sin_phi;

												if(k==2&& l == 3) c=sin_phi*sin_theta;

												if(k==3&& l == 0) c=alpha*cos_theta;

												if(k==3&& l == 1) c=-(gamma*sin_theta);

												if(k==3&& l == 2) c=beta*sin_theta;

												if(k==3&& l == 3) c=cos_theta;

												tHa[k][l]=c;
										}
								}
						}
						else if(flg == 1){
								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												for(int i = 0; i < 2; i++){
														Double_t c=0;

														if(i==0&&k==0&& l == 0) c=1;
														if(i==1&&k==0&& l == 0) c=0;

														if(i==0&&k==0&& l == 3) c=1.;
														if(i==1&&k==0&& l == 3) c=0;

														if(i==0&&k==1&& l == 0) c=1.*cos_phi*sin_theta;
														if(i==1&&k==1&& l == 0) c=0;

														if(i==0&&k==1&& l == 1) c=1.*cos_phi*cos_theta;
														if(i==1&&k==1&& l == 1) c= - 1.*sin_phi;

														if(i==0&&k==1&& l == 2) c=-(1.*cos_phi*cos_theta);
														if(i==1&&k==1&& l == 2) c= - 1.*sin_phi;

														if(i==0&&k==1&& l == 3) c=cos_phi*sin_theta;
														if(i==1&&k==1&& l == 3) c=0;

														if(i==0&&k==2&& l == 0) c=1.*sin_phi*sin_theta;
														if(i==1&&k==2&& l == 0) c=0;

														if(i==0&&k==2&& l == 1) c=1.*cos_phi ;
														if(i==1&&k==2&& l == 1) c= 1.*cos_theta*sin_phi;

														if(i==0&&k==2&& l == 2) c=1.*cos_phi;
														if(i==1&&k==2&& l == 2) c=- 1.*cos_theta*sin_phi;

														if(i==0&&k==2&& l == 3) c=sin_phi*sin_theta;
														if(i==1&&k==2&& l == 3) c=0;

														if(i==0&&k==3&& l == 0) c=1.*cos_theta;
														if(i==1&&k==3&& l == 0) c=0;

														if(i==0&&k==3&& l == 1) c=-(1.*sin_theta);
														if(i==1&&k==3&& l == 1) c=0;

														if(i==0&&k==3&& l == 2) c=1.*sin_theta;
														if(i==1&&k==3&& l == 2) c=0;

														if(i==0&&k==3&& l == 3) c=cos_theta;
														if(i==1&&k==3&& l == 3) c=0;

														tHb[k][l][i]=c;
												}
										}
								}
						}

						else if(flg == 2){
								double gamma = sqrt(1-alpha*alpha)*cos(xiphi);
								double beta = sqrt(1-alpha*alpha)*sin(xiphi);

								for(int k=0;k<4;k++){
										for(int l=0;l<4;l++){
												for(int i=0;i<2;i++){
														Double_t c=1;

														if(i==0&&k==0&& l == 0) c=1;
														if(i==1&&k==0&& l == 0) c=1;

														if(i==0&&k==0&& l == 3) c=alpha;
														if(i==1&&k==0&& l == 3) c=alpha;

														if(i==0&&k==1&& l == 0) c=alpha;
														if(i==1&&k==1&& l == 0) c=alpha;

														if(i==0&&k==1&& l == 1) c=gamma;
														if(i==1&&k==1&& l == 1) c=beta;

														if(i==0&&k==1&& l == 2) c=(beta);
														if(i==1&&k==1&& l == 2) c= gamma;

														if(i==0&&k==1&& l == 3) c=1;
														if(i==1&&k==1&& l == 3) c=1;

														if(i==0&&k==2&& l == 0) c=alpha;
														if(i==1&&k==2&& l == 0) c=alpha;

														if(i==0&&k==2&& l == 1) c=beta;
														if(i==1&&k==2&& l == 1) c=gamma;

														if(i==0&&k==2&& l == 2) c=gamma;
														if(i==1&&k==2&& l == 2) c=beta;

														if(i==0&&k==2&& l == 3) c=1;
														if(i==1&&k==2&& l == 3) c=1;

														if(i==0&&k==3&& l == 0) c=alpha;
														if(i==1&&k==3&& l == 0) c=alpha;

														if(i==0&&k==3&& l == 1) c=(gamma);
														if(i==1&&k==3&& l == 1) c=(gamma);

														if(i==0&&k==3&& l == 2) c=beta;
														if(i==1&&k==3&& l == 2) c=beta;

														if(i==0&&k==3&& l == 3) c=1;
														if(i==1&&k==3&& l == 3) c=1;

														tHc[k][l][i]=c;
												}
										}
								}
						}
				}





				~AADecay12() {}
				Double_t A(Int_t i, Int_t j) const {return tHa[i][j];}
				Double_t Ab(Int_t i, Int_t j, Int_t k) const {return tHb[i][j][k];}
				Double_t Ac(Int_t i, Int_t j, Int_t k) const {return tHc[i][j][k];}
		private:
				double tHa[4][4];
				double tHb[4][4][2];
				double tHc[4][4][2];
};
#endif
