#ifndef MN_AngDisXiXi_H_
#define MN_AngDisXiXi_H_
#include "AAProd1212.hh"
#include "AADecay12.hh"
#include <cmath>
using namespace std;
class AngDisXiXi {
		//e+e- -> (Ss->L pi-) (S->nb pi+) (L->p pi-)
		public:
				AngDisXiXi(double Jpsi_alpha, double Jpsi_phi,
								double Xi_alpha, double Xi_phi, 
								double Xibar_alpha, double Xibar_phi, 
								double L_alpha,  
								double Lbar_alpha):
								p_Jpsi_alpha(Jpsi_alpha), p_Jpsi_phi(Jpsi_phi), 
								p_Xi_alpha(Xi_alpha), p_Xi_phi(Xi_phi), 
								p_Xibar_alpha(Xibar_alpha), p_Xibar_phi(Xibar_phi), 
								p_L_alpha(L_alpha), 
								p_Lbar_alpha(Lbar_alpha)
		{}   
				~AngDisXiXi() {}
				double operator()(double theta,
								double LThe,double LPhi,double LbThe,
								double LbPhi,double pThe,double pPhi, 
								double apThe, double apPhi) 
				{
						AAProd1212 prod(p_Jpsi_alpha, p_Jpsi_phi, theta);
						AADecay12 Xi(p_Xi_alpha, p_Xi_phi, LThe, LPhi);
						AADecay12 Xibar(p_Xibar_alpha, p_Xibar_phi, LbThe, LbPhi);
						AADecay12 L(p_L_alpha, 0, pThe, pPhi);
						AADecay12 Lbar(p_Lbar_alpha, 0, apThe, apPhi);
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
		private:
				double  p_Jpsi_alpha;
				double 	p_Jpsi_phi;
				double 	p_Xi_alpha;
				double  p_Xi_phi;
				double 	p_Xibar_alpha;
				double  p_Xibar_phi;
				double  p_L_alpha;
				double  p_Lbar_alpha;
};
#endif // MN_AngDisSsS_H_
