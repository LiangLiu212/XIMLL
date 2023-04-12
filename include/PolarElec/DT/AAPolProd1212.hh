#ifndef AA_PolProd1212_H_
#define AA_PolProd1212_H_
#include "TMath.h"
/*
  AAPolProd1212 : e+e- -> 1/2 1/2 polarized e exclusive 
     2018-06-09 v2.00 (AK) separate terms to speed up calculations
     2021-06-11 v3.00 (PA)

 */
class AAPolProd1212 {
  public:

  Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
  Double_t CT(Int_t k, Int_t l,Int_t it) const {return tT[it][l+4*k];}
  AAPolProd1212(double alpha_psi, double DeltaPhi_psi, double P_e, double th) {
    // decay 1/2 -> 1/2 + 0
    Double_t v=TMath::Sqrt(1-alpha_psi*alpha_psi);
    Double_t bP=v*TMath::Sin(DeltaPhi_psi);
    Double_t gP=v*TMath::Cos(DeltaPhi_psi);
    
    Table(alpha_psi,bP,gP,P_e,th);
  }
  AAPolProd1212(double alpha_psi, double bP, double gP, double P_e, double th) {
    // decay 1/2 -> 1/2 + 0
    Table(alpha_psi,bP,gP,P_e,th);
    iTest=kFALSE;
    
  }
  void Table(double alpha_psi, double DeltaPhi_psi, double P_e, double th) {
    Double_t v=TMath::Sqrt(1-alpha_psi*alpha_psi);
    Double_t gP=v*TMath::Cos(DeltaPhi_psi);
    Double_t bP=v*TMath::Sin(DeltaPhi_psi);
    Table(alpha_psi,bP,gP,P_e,th);
    iTest=kFALSE;
  }
  void SetTest(int tst=kTRUE){
    iTest=tst;
  }



  void Table(double alpha_psi, double bP, double gP,double P_e,double th) {
    // decay 1/2 -> 1/2 + 0
    Double_t sinth=TMath::Sin(th);
    Double_t costh=TMath::Cos(th);
    Double_t costh2=costh*costh;
    Double_t sinth2=1-costh2;
    if(iTest){
      sinth=1;
      costh=1;
      costh2=1;
      sinth2=1;
    }
    // for(int k=0;k<4;k++){
    //   for(int l=0;l<4;l++){
    //     Double_t c=0;
    //     Double_t p0=0;
    //     Double_t pa=0;
    //     Double_t pb=0;
    //     Double_t pg=0;
    //     Double_t pol=0;
    //     if(k==0&&l==0){p0=1;pa=costh2;}
    //     if(k==0&&l==1)pol = gP*sinth;
    //     if(k==0&&l==2)pb = sinth*costh;
    //     if(k==0&&l==3)pol = (1.+alpha_psi)*costh;
    //     if(k==1&&l==0)pol = gP*sinth;
    //     if(k==1&&l==1)p0 = sinth2;
    //     if(k==1&&l==3)pg = sinth*costh;
    //     if(k==2&&l==0)pb =-sinth*costh;
    //     if(k==2&&l==2)pa=sinth2;
    //     if(k==2&&l==3)pol=-bP*sinth;
    //     if(k==3&&l==0)pol=-(1.+alpha_psi)*costh;
    //     if(k==3&&l==1)pg =-sinth*costh;
    //     if(k==3&&l==2)pol=-bP*sinth;
    //     if(k==3&&l==3){pa=-1;p0=-costh2;}
        
    //     thC[k][l]=p0 + alpha_psi*pa + bP*pb + gP*pg + pol*P_e;
    //     Int_t ii=l+4*k;
    //     tT[0][ii]=p0;
    //     tT[1][ii]=pa;
    //     tT[2][ii]=pb;
    //     tT[3][ii]=pg;
    //     tT[4][ii]=pol;
    //   }
    // }
for(int k=0;k<4;k++){
  for(int l=0;l<4;l++){
    Double_t c=0;
    Double_t p0=0;
    Double_t pa=0;
    Double_t pb=0;
    Double_t pg=0;
    Double_t pp=0;
    Double_t pap=0;
    Double_t pbp=0;
    Double_t pgp=0;
    if(k==0&&l==0){p0=1;pa=costh2;}
    if(k==0&&l==1)pgp=sinth;
    if(k==0&&l==2)pb=sinth*costh;
    if(k==0&&l==3){pp=costh;pap=costh;}
    if(k==1&&l==0)pgp=sinth;
    if(k==1&&l==1)p0=sinth2;
    if(k==1&&l==3)pg=sinth*costh;
    if(k==2&&l==0)pb=-sinth*costh;
    if(k==2&&l==2)pa=sinth2;
    if(k==2&&l==3)pbp=-sinth;
    if(k==3&&l==0){pp=-costh;pap=-costh;}
    if(k==3&&l==1)pg=-sinth*costh;
    if(k==3&&l==2)pbp=-sinth;
    if(k==3&&l==3){pa=-1;p0=-costh2;}
    thC[k][l]=p0+alpha_psi*pa+bP*pb+gP*pg+alpha_psi*P_e*pap+bP*P_e*pbp+gP*P_e*pgp+P_e*pp;
    Int_t ii=l+4*k;
    tT[0][ii]=p0; //*1
    tT[1][ii]=pa; //*aP
    tT[2][ii]=pb; //*bP
    tT[3][ii]=pg; //*gP
    tT[4][ii]=pp; //*Pe
    tT[5][ii]=pap; //*aP*Pe
    tT[6][ii]=pbp; //*bP*Pe
    tT[7][ii]=pgp; //*gP*Pe
  }
}

  }
  ~AAPolProd1212() {}
private:

  double thC[4][4];
  // individual terms as fcn of angles:
  // double tT[5][16];
  double tT[8][16];
  Bool_t iTest;
};
#endif /* AA_PolProd1212_H_ */
