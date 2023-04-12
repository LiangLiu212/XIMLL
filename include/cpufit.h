#ifndef CPU_FIT_H
#define CPU_FIT_H

#include "readfile.h"
#include "PolarElec/DT/AngDisPolXiXi.hh"
#include "PolarElec/XibarST/AngDisXXbar.hh"
#include "PolarElec/XiST/AngDisXX.hh"
#include "TThread.h"

using namespace std;
class cpufit : public readfile{
		public:
				cpufit(){
						m_file.clear();
						m_year.clear();
						m_type.clear();
						m_sample.clear();
						NUM = 10000000;
						cut_LmdDL = 0.0;
						cut_XiDL = 0.0;
						cut_XiCosTheta = 0.84;
						cut_mXi = 0.011;
						cut_mLmd1 = 0.011;
						cut_chi2kmf = 200.;
						cut_chi2Xi = 999.;
						cut_chi2Lmd = 999.;
						cut_mn2 = 0.955;
						cut_mn1 = 0.925;
						fit_step=0;
						flag_fast = false;
						double tmpNBKG[4][2] = {{861.385, 657.545}, {4425.2, 3421.11}, {16359.6, 12720.5}, {16589.0, 12382.8}};
						for(int i = 0; i < 4; i++)
								for(int j =0; j < 2; j++)
										NBKG[i][j] = tmpNBKG[i][j];
				}
				~cpufit() {;}
				virtual void InitialMemory();
				virtual void FreeMemory();
				virtual void IOReadData(const int index, const int MM);
				virtual double IOfcnmll(double *pp);
				virtual void IOMLL();
				virtual void *hFCN0(void *ptr);
				virtual void *hFCN1(void *ptr);
				virtual void *hFCN2(void *ptr);
				virtual void *hFCN3(void *ptr);
				virtual void *hFCN4(void *ptr);
				virtual void *hFCN5(void *ptr);
				virtual void *hFCN6(void *ptr);
				virtual void *hFCN7(void *ptr);
		private:
				virtual int IOreadData(const int n, double **para, const int  index, const int MM);
				AngDis *Pangdis[4][2][8]; 
				AngDisbar *Nangdis[4][2][8]; 
				AngDis *PangDis;
				AngDisbar *NangDis;

				Double_t polarisation     =  0.;
				Double_t alpha_jpsi       =  0.58;  // alpha_J/Psi 
				Double_t dphi_jpsi        =  1.2;    // relative phase, Dphi_J/Psi

				Double_t alpha_xi         = -0.37;  // alpha(Xi->Lambda pi-) 
				Double_t phi_xi           =  0.02;   // Xi- phi 
				Double_t alpha_lam        =  0.754;  // alpha (Lambda->p pi-)
				Double_t alpha_xibar      =  0.37;   // alpha (Xibar->Lambda pi+)
				Double_t alpha_lambar     = -0.754; // alpha (Lambda->p pi-)
				Double_t phi_xibar        = -0.02;  // Xi- phi 
				TThread *t[8];
};
#endif
