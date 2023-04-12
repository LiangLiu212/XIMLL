#ifndef AA_Decay12bar_H_
#define AA_Decay12bar_H_
#include "TMath.h"
/*
  AADecay12 : Decay 1/2 -> 1/2 + 0
     2018-06-09 v2.00 (AK) separate terms to speed up calculations

 */
class AADecay12bar {
  public:
  AADecay12bar(double aM, double pM, double th, double ph) {
    Double_t v=TMath::Sqrt(1-aM*aM);
    Double_t gM=v*TMath::Cos(pM);
    Double_t bM=v*TMath::Sin(pM);
    Table(aM, bM, gM, th,  ph);
    iTest=kFALSE;
  }
  AADecay12bar(double aM, double bM, double gM, double th, double ph) {
    Table(aM, bM, gM, th,  ph);
    iTest=kFALSE;
  }
  void SetTest(int tst=kTRUE){
    iTest=tst;
  }
  void Table(double aM, double pM, double th, double ph) {
    Double_t v=TMath::Sqrt(1-aM*aM);
    Double_t gM=v*TMath::Cos(pM);
    Double_t bM=v*TMath::Sin(pM);
    Table(aM, bM, gM, th,  ph);
  }

  void Table(double aM, double bM, double gM, double th, double ph) {
    // decay 1/2 -> 1/2 + 0
    fA=aM;
    fB=bM;
    fG=gM;
    Double_t sp=TMath::Sin(ph);
    Double_t cp=TMath::Cos(ph);
    Double_t st=TMath::Sin(th);
    Double_t ct=TMath::Cos(th);
    if(iTest){
      sp=1;
      cp=1;
      st=1;
      ct=1;
    }
    for(int k=0;k<4;k++){
      for(int l=0;l<4;l++){
        Double_t a=0;
        Double_t p0=0;
        Double_t pa=0;
        Double_t pb=0;
        Double_t pg=0;
        if(k==0&&l==0)p0= 1;
        if(k==0&&l==3)pa= 1;
        if(k==1&&l==0)pa= cp*st;
        if(k==1&&l==1){pg= ct*cp; pb=-sp;}
        if(k==1&&l==2){pb=-ct*cp;pg=-sp;}
        if(k==1&&l==3)p0= cp*st;
        if(k==2&&l==0)pa= sp*st;
        if(k==2&&l==1){pb= cp;pg=+ct*sp;}
        if(k==2&&l==2){pg= cp;pb=-ct*sp;}
        if(k==2&&l==3)p0= sp*st;
        if(k==3&&l==0)pa= ct;
        if(k==3&&l==1)pg=-st;
        if(k==3&&l==2)pb= st;
        if(k==3&&l==3)p0= ct;
        tHa[k][l]=p0+aM*pa+pb*bM+pg*gM;
        Int_t ii=l+4*k;
        tT[0][ii]=p0;
        tT[1][ii]=pa;
        tT[2][ii]=pb;
        tT[3][ii]=pg;
      }
    }
  }      
  ~AADecay12bar() {}
  Double_t A(Int_t k, Int_t l) const {return tHa[k][l];}
  Double_t D1Phi(Int_t k, Int_t l) const
  {
    // Calculate derivate d phi
       Int_t ii=l+4*k;
       Double_t dbdp=fG;
       Double_t dgdp=-fB;
       Double_t deriv=dbdp*tT[2][ii]+dgdp*tT[3][ii];
    return deriv;
  }
  Double_t D1A(Int_t k, Int_t l) const
  {
    // Calculate derivate d alpha
       Int_t ii=l+4*k;
       Double_t dbda=-fA*fB/(1-fA*fA);
       Double_t dgda=-fA*fG/(1-fA*fA);
       Double_t deriv=tT[1][ii]+dbda*tT[2][ii]+dgda*tT[3][ii];
    return deriv;
  }
  Double_t D2A(Int_t k, Int_t l) const
  {
    // Calculate second derivate d alpha
       Int_t ii=l+4*k;
       Double_t DD=(1-fA*fA)*(1-fA*fA);
       Double_t dbda=-fB/DD;
       Double_t dgda=-fG/DD;
       Double_t deriv=dbda*tT[2][ii]+dgda*tT[3][ii];
    return deriv;
  }
  Double_t AT(Int_t k, Int_t l,Int_t it) const {return tT[it][l+4*k];}
private:

  double tHa[4][4];
  // individual terms :
  double tT[4][16];
  Bool_t iTest;
  Double_t fA,fB,fG;
};
#endif
