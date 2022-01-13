#ifndef AA_Decay32_H_
#define AA_Decay32_H_
#include "TMath.h"
class AADecay32 {
		public:
				AADecay32(const double alpha, const double beta, const double gamma, const double theta, const double phi) {
						// decay 3/2 -> 1/2 + 0
						for(int k=0;k<16;k++){
								for(int l=0;l<4;l++){
										Double_t c=0;

										if(k==0&& l == 0) c=1;

										if(k==0&& l == 3) c=alpha;

										if(k==1&& l == 0) c=sqrt(0.6)*alpha*sin(phi)*sin(theta);

										if(k==1&& l == 1) c=2*sqrt(0.6)*(beta*cos(phi) + gamma*cos(theta)*sin(phi));

										if(k==1&& l == 2) c=2*sqrt(0.6)*(gamma*cos(phi) + beta*cos(theta)*sin(phi));

										if(k==1&& l == 3) c=sqrt(0.6)*sin(phi)*sin(theta);

										if(k==6&& l == 0) c=-(sqrt(3)*(1 + 3*cos(2*theta)))/4.;

										if(k==6&& l == 3) c=-(sqrt(3)*alpha*(1 + 3*cos(2*theta)))/4.;

										if(k==7&& l == 0) c=-3*cos(phi)*cos(theta)*sin(theta);

										if(k==7&& l == 3) c=-3*alpha*cos(phi)*cos(theta)*sin(theta);

										if(k==8&& l == 0) c=(-3*cos(2*phi)*Power(sin(theta),2))/2.;

										if(k==8&& l == 3) c=(-3*alpha*cos(2*phi)*Power(sin(theta),2))/2.;

										if(k==10&& l == 0) c=-9*alpha*cos(phi)*cos(theta)*sin(phi)*Power(sin(theta),2);

										if(k==10&& l == 1) c=-3*beta*cos(2*phi)*cos(theta)*sin(theta) - (3*gamma*(1 + 3*cos(2*theta))*sin(2*phi)*sin(theta))/4.;

										if(k==10&& l == 2) c=-3*gamma*cos(2*phi)*cos(theta)*sin(theta) + (3*beta*(1 + 3*cos(2*theta))*sin(2*phi)*sin(theta))/4.;

										if(k==10&& l == 3) c=-9*cos(phi)*cos(theta)*sin(phi)*Power(sin(theta),2);

										if(k==11&& l == 0) c=(-9*alpha*(3 + 5*cos(2*theta))*sin(phi)*sin(theta))/(4.*sqrt(10));

										if(k==11&& l == 1) c=(-3*beta*cos(phi)*(3 + 5*cos(2*theta)))/(4.*sqrt(10)) - (3*gamma*(cos(theta) + 15*cos(3*theta))*sin(phi))/(8.*sqrt(10));

										if(k==11&& l == 2) c=(-3*gamma*cos(phi)*(3 + 5*cos(2*theta)))/(4.*sqrt(10)) + (3*beta*(cos(theta) + 15*cos(3*theta))*sin(phi))/(8.*sqrt(10));

										if(k==11&& l == 3) c=(-9*(3 + 5*cos(2*theta))*sin(phi)*sin(theta))/(4.*sqrt(10));

										tHa[k][l]=c;
								}
						}
				}   
				~AADecay32() {}
				Double_t A(Int_t i, Int_t j) const {return tHa[i][j];}
		private:
				double tHa[16][4];
};
#endif
