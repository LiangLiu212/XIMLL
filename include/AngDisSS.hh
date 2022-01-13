#ifndef MN_AngDisSsS_H_
#define MN_AngDisSsS_H_
#include "AAProd1212.hh"
#include "AAProd3212.hh"
#include "AADecay12.hh"
#include "AADecay32.hh"
#include <cmath>
using namespace std;
class AngDisSsS {
		//e+e- -> (Ss->L pi-) (S->nb pi+) (L->p pi-)
		public:
				AngDisSsS(double Jpsi_eta, double Jpsi_xi, double Jpsi_phi2, double Jpsi_phi3,
								double Ss_alpha, double Ss_beta, double Ss_gamma,
								double L_alpha, double L_beta, double L_gamma,
								double aS_alpha, double aS_beta, double aS_gamma) :
								p_Jpsi_eta(Jpsi_eta), p_Jpsi_xi(Jpsi_xi), p_Jpsi_phi2(Jpsi_phi2), p_Jpsi_phi3(Jpsi_phi3), 
								p_Ss_alpha(Ss_alpha), p_Ss_beta(Ss_beta), p_Ss_gamma(Ss_gamma),
								p_L_alpha(L_alpha), p_L_beta(L_beta), p_L_gamma(L_gamma),
								p_aS_alpha(aS_alpha), p_aS_beta(aS_beta), p_aS_gamma(aS_gamma),
		{}   
				~AngDisSsS() {}
				double operator()(double theta,
								double SsThe,double SsPhi,double LThe,
								double LPhi,double aSThe,double aSPhi) 
				{
						AAProd3212 prod(p_Jpsi_eta, p_Jpsi_xi, p_Jpsi_phi2, p_Jpsi_phi3, theta);
						AADecay32 Ss(p_Ss_alpha, p_Ss_beta, p_Ss_gamma, SsThe, SsPhi);
						AADecay12 L(p_L_alpha, p_L_beta, p_L_gamma, LThe, LPhi);
						AADecay12 aS(p_aS_alpha, p_aS_beta, p_aS_gamma, aSThe, aSPhi);
						double tep=0;
						for(int mu=0; mu<16;mu++){// Ss loop
								for(int nu=0;nu<4;nu++){// aS loop
										for(int k=0;k<4;k++)
												tep+=prod.C(nu, mu)*
														Ss.A(mu,k)*L.A(k,0)*
														aS.A(nu, 0);

								//		tep+=prod.C(j1,j2)*
								//				dS.A(j1,j2)*dSb.A(j1,j2);
								}
						}
						return tep;
				}
		private:
				double  p_Jpsi_eta;
				double 	p_Jpsi_xi;
				double 	p_Jpsi_phi2;
				double 	p_Jpsi_phi3;
				double 	p_Ss_alpha;
				double  p_Ss_beta;
				double  p_Ss_gamma;
				double  p_L_alpha;
				double  p_L_beta;
				double  p_L_gamma;
				double  p_aS_alpha;
				double  p_aS_beta;
				double  p_aS_gamma;
};
#endif // MN_AngDisSsS_H_
