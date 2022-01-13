#ifndef AA_Prod1212_H_
#define AA_Prod1212_H_
#include "TMath.h"
class AAProd1212 {
		public:
				Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
				AAProd1212(const double eta, const double xi, const double phi2, const double phi3 double Jpsi_the) {
						// decay Jpsi -> 3/2 + 1/2
						Double_t h1 = sin(eta)*cos(xi)/sqrt(2);
						Double_t h2 = sin(eta)*sin(xi);
						Double_t h3 = cos(eta);
					//	Double_t v=TMath::Sqrt(1-alpha_Jpsi*alpha_Jpsi);
					//	Double_t st=TMath::Sin(Jpsi_the);
					//	Double_t ct=TMath::Cos(Jpsi_the);
					//	Double_t cp = TMath::Cos(delta_phi);
					//	Double_t sp = TMath::Sin(delta_phi);
						for(int k=0;k<4;k++){
								for(int l=0;l<16;l++){
										Double_t c=0;

										if(k==0&& l == 0) c=((power(h2,2) + power(h3,2))*(3 + cos(2*Jpsithe)))/2. + 2*power(h1,2)*power(sin(Jpsithe),2);

										if(k==0&& l == 1) c=-(sqrt(0.4)*h1*sin(2*Jpsithe)*(2*sqrt(3)*h2*sin(phi2) + 3*h3*sin(phi3)))/3.;

										if(k==0&& l == 6) c=(-((h2 - h3)*(h2 + h3)*(3 + cos(2*Jpsithe))) - 4*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(3));

										if(k==0&& l == 7) c=sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe);

										if(k==0&& l == 8) c=(2*h2*h3*cos(phi2 - phi3)*power(sin(Jpsithe),2))/sqrt(3);

										if(k==0&& l == 10) c=(2*h2*h3*power(sin(Jpsithe),2)*sin(phi2 - phi3))/sqrt(3);

										if(k==0&& l == 11) c=(2*h1*sin(2*Jpsithe)*(3*h2*sin(phi2) - sqrt(3)*h3*sin(phi3)))/(3.*sqrt(5));

										if(k==1&& l == 2) c=-(sqrt(2.0/15)*h1*h2*cos(phi2)*sin(2*Jpsithe));

										if(k==1&& l == 3) c=-((h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3))/sqrt(5)) - (2*(2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2))/sqrt(15);

										if(k==1&& l == 4) c=sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3);

										if(k==1&& l == 5) c=-((h2*h3*(3 + cos(2*Jpsithe))*sin(phi2 - phi3))/sqrt(3));

										if(k==1&& l == 12) c=sqrt(1.2)*h1*h2*cos(phi2)*sin(2*Jpsithe);

										if(k==1&& l == 13) c=-(sqrt(2.0/15)*h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3)) + sqrt(0.4)*(2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2);

										if(k==1&& l == 14) c=-(sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe));

										if(k==1&& l == 15) c=-(sqrt(2.0/3)*power(h3,2)*power(sin(Jpsithe),2));

										if(k==2&& l == 0) c=-2*sqrt(2)*h1*h2*cos(Jpsithe)*sin(Jpsithe)*sin(phi2);

										if(k==2&& l == 1) c=(h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3))/sqrt(5) - (2*(-2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2))/sqrt(15);

										if(k==2&& l == 6) c=sqrt(2.0/3)*h1*h2*sin(2*Jpsithe)*sin(phi2);

										if(k==2&& l == 7) c=-((h2*h3*(3 + cos(2*Jpsithe))*sin(phi2 - phi3))/sqrt(3));

										if(k==2&& l == 8) c=sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3);

										if(k==2&& l == 9) c=sqrt(2.0/3)*power(h3,2)*power(sin(Jpsithe),2);

										if(k==2&& l == 10) c=sqrt(2.0/3)*h1*h3*cos(phi3)*sin(2*Jpsithe);

										if(k==2&& l == 11) c=sqrt(2.0/15)*h2*h3*(3 + cos(2*Jpsithe))*cos(phi2 - phi3) + sqrt(0.4)*(-2*power(h1,2) + power(h2,2))*power(sin(Jpsithe),2);

										if(k==3&& l == 2) c=(-((power(h2,2) - 3*power(h3,2))*(3 + cos(2*Jpsithe))) + 4*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(15));

										if(k==3&& l == 3) c=(sqrt(0.4)*h1*(-2*sqrt(3)*h2*cos(phi2) + 3*h3*cos(phi3))*sin(2*Jpsithe))/3.;

										if(k==3&& l == 4) c=(2*h2*h3*power(sin(Jpsithe),2)*sin(phi2 - phi3))/sqrt(3);

										if(k==3&& l == 5) c=-(sqrt(2.0/3)*h1*h3*sin(2*Jpsithe)*sin(phi3));

										if(k==3&& l == 12) c=((3*power(h2,2) + power(h3,2))*(3 + cos(2*Jpsithe)) - 12*power(h1,2)*power(sin(Jpsithe),2))/(2.*sqrt(15));

										if(k==3&& l == 13) c=(2*h1*(3*h2*cos(phi2) + sqrt(3)*h3*cos(phi3))*sin(2*Jpsithe))/(3.*sqrt(5));

										if(k==3&& l == 14) c=(2*h2*h3*cos(phi2 - phi3)*power(sin(Jpsithe),2))/sqrt(3);


/*
										if(k==0&&l==0)c=1;
										if(k==1&&l==1)c=alpha_Jpsi*ct*ct;
										if(k==1&&l==2)c=st*st;
										if(k==1&&l==3)c=-alpha_Jpsi*st*st;
										if(k==2&&l==1)c=ct*ct;
										if(k==2&&l==2)c=v*sp*st*ct;
										if(k==2&&l==3)c=v*cp*st*ct;
										if(k==3&&l==1)c=alpha_Jpsi;
										if(k==3&&l==2)c=v*cp*st*ct;
										if(k==3&&l==3)c=v*sp*st*ct; */


										thC[k][l]=c;
								}
						}
				}
				~AAProd1212() {}
		private:
				double thC[4][16];
};
#endif /* AA_Prod1212_H_ */
