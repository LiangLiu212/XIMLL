#ifndef MN_AngDisXX_H_
#define MN_AngDisXX_H_
#include "AAProd1212.hh"
#include "AADecay12.hh"
#include <cmath>
#include "TMath.h"
/*
  AngDis  - for e+e- -> Xi Xibar exclusive distributions including
            decay chains
      
 */
using namespace std;

class AngDis {

  public:

    AngDis(double ap, double pp, double pe, double aX, double pX, double aL) :
      tHaPsi(ap),tHpPsi(pp), tHpE(pe), tHaX(aX), tHpX(pX), tHaL(aL)
  {

      double dx = 1.e-10;
    prod  = new AAProd1212(tHaPsi,tHpPsi,tHpE,0);
    dX    = new AADecay12(tHaX,tHpX,0,0);
    dL    = new AADecay12(tHaL,0,0,0);
    proda  = new AAProd1212(tHaPsi+dx,tHpPsi,tHpE,0);
    prodp  = new AAProd1212(tHaPsi,tHpPsi+dx,tHpE,0);
    prode  = new AAProd1212(tHaPsi,tHpPsi,tHpE+dx,0);
    dXa    = new AADecay12(tHaX+dx,tHpX,0,0);
    dXp    = new AADecay12(tHaX,tHpX+dx,0,0);
    dLa    = new AADecay12(tHaL+dx,0,0,0);
    iTest=0;
    iNZero=0;
  } 
 // void Set(double ap, double pp, double aL, double aLb, double pXi)
  void Set(double ap, double pp, double pe, double aX, double pX, double aL)
  {
    // Set parameters
      double dx = 1.e-10;
    tHaPsi=ap;
    tHpPsi=pp;
    tHpE=pe;
    tHaX=aX;
    tHpX=pX;
    tHaL=aL;
    prod->Table(tHaPsi,tHpPsi,tHpE,0);
    dX->Table(tHaX,tHpX,0,0);
    dL->Table(tHaL,0,0,0);
    proda->Table(tHaPsi+dx,tHpPsi,tHpE,0);
    prodp->Table(tHaPsi,tHpPsi+dx,tHpE,0);
    prode->Table(tHaPsi,tHpPsi,tHpE+dx,0);
    dXa->Table(tHaX+dx,tHpX,0,0);
    dXp->Table(tHaX,tHpX+dx,0,0);
    dLa->Table(tHaL+dx,0,0,0);
  }     
  ~AngDis() 
  {
    delete prod;
    delete dX;
    delete dL;
    delete proda;
    delete prodp;
    delete prode;
    delete dXa;
    delete dXp;
    delete dLa;
    
  }
  double  W(double th, double th1l,
	   double ph1l,double th1p,
           double ph1p)

    {
      // Calculate weight for event specified by the variables 
      prod->Table(tHaPsi,tHpPsi,tHpE,th);
      dX->Table(tHaX,tHpX,th1l,ph1l);
      dL->Table(tHaL,0,th1p,ph1p);
      double tep=0;
      for(int i1=0;i1<4;i1++){// Xi loop
          for(int j1=0;j1<4;j1++){// L loop
	      tep += prod->C(i1,0)*
	      dX->A(i1,j1)*dL->A(j1,0);
	    }
      }
      return tep;
    }
    double  DW1(double th,double th1l,
		double ph1l,double th1p,
           	double ph1p,int ik) 

    {
      // ik is variable to take derivate
      //                                 1 :aXi
      //                                 2 :phiXi
      //                                 3 :aL
      //                                 4 :apsi
      //                                 5 :DPhi
      //                                 6 :Pe
      // Calculate first derivative weight for event specified by the variables 
      double dx = 1.e-10;
      prod->Table(tHaPsi,tHpPsi,tHpE,th);
      dX->Table(tHaX,tHpX,th1l,ph1l);
      dL->Table(tHaL,0,th1p,ph1p);
      proda->Table(tHaPsi+dx,tHpPsi,tHpE,th);
      prodp->Table(tHaPsi,tHpPsi+dx,tHpE,th);
      prode->Table(tHaPsi,tHpPsi,tHpE+dx,th);
      dXa->Table(tHaX+dx,tHpX,th1l,ph1l);
      dXp->Table(tHaX,tHpX+dx,th1l,ph1l);
      dLa->Table(tHaL+dx,0,th1p,ph1p);
      double tep=0;
      for(int i1=0;i1<4;i1++){// Xi loop
      	  for(int j1=0;j1<4;j1++){// L loop
	      Double_t d3=dX->A(i1,j1);
	      if(ik==1)d3=dXa->A(i1,j1);
	      if(ik==2)d3=dXp->A(i1,j1);
	      Double_t d4=dL->A(j1,0);
	      if(ik==3)d4=dLa->A(j1,0);
	      Double_t cxx=prod->C(i1,0);
	      if(ik==4)cxx=proda->C(i1,0);
	      if(ik==5)cxx=prodp->C(i1,0);
	      if(ik==6)cxx=prode->C(i1,0);
	      tep += (cxx*d3*d4-(prod->C(i1,0)*dX->A(i1,j1)*dL->A(j1,0)))/dx;
	}
      }
      return tep;
    }
private:

double tHaPsi;
double tHpPsi;
double tHpE;
double tHaX;
double tHpX;
double tHaL;
  // integral variables:
  Double_t fInt[7776];Double_t fN;
  AAProd1212 *prod;
  AADecay12 *dX;
  AADecay12 *dL;
  AAProd1212 *proda;
  AAProd1212 *prodp;
  AAProd1212 *prode;
  AADecay12 *dXa;
  AADecay12 *dXp;
  AADecay12 *dLa;
  // Test flag to determine # of non zero terms
  Bool_t iTest;
  // Non zero terms and indicies:
  Int_t iNZero;
  Int_t gID[7776];
  Int_t gDI[7776];
  Int_t  iI[7776];
  Int_t iK1[7776];
  Int_t iK2[7776];
  Int_t iJ1[7776];
  Int_t iJ2[7776];
  Int_t iT;
  Int_t Id(UInt_t ig,UInt_t ip)
  {
    // Decode global 0<=ig<6^5 
    // 0<=id<6  <- id of term for each process  0<=ip<5
    UInt_t k=2*ip;
    int ii=0x3&(ig>>k);
    return ii;
  }
};

#endif // MN_AngDisXX_H_

