#ifndef AA_Prod1212bar_H_
#define AA_Prod1212bar_H_
#include "TMath.h"
/*
  AAProd1212 : e+e- -> 1/2 1/2 polarized exclusive 
     2018-06-09 v2.00 (AK) separate terms to speed up calculations -> unpolarized case
     2021-02-04  (VB) 

 */
class AAProd1212bar {
  public:

  Double_t C(Int_t i, Int_t j) const {return thC[i][j];}
  Double_t DCA(Int_t k, Int_t l) const {
    // Calculate derivate d alpha
    Int_t ii=l+4*k;
    Double_t dbda=-fA*fB/(1-fA*fA);
    Double_t dgda=-fA*fG/(1-fA*fA);
    Double_t dpda=fP;
    Double_t deriv=tT[1][ii]+dbda*tT[2][ii]+dgda*tT[3][ii]+dpda*tT[5][ii];
    return deriv;
  }
  Double_t DCPhi(Int_t k, Int_t l) const
  {
    // Calculate derivate d phi
    Int_t ii=l+4*k;
    Double_t dbdp=fG;
    Double_t dgdp=-fB;
    Double_t dbpdp=fG*fP;
    Double_t dgpdp=-fB*fP;
    Double_t deriv=dbdp*tT[2][ii]+dgdp*tT[3][ii]+dbpdp*tT[6][ii]+dgpdp*tT[7][ii];
    return deriv;
  }
  Double_t DCPe(Int_t k, Int_t l) const
  {
    // Calculate derivate d Pe
    Int_t ii=l+4*k;
    Double_t dbpdp=fB;
    Double_t dgpdp=fG;
    Double_t dapdp=fA;
    Double_t deriv=tT[4][ii]+dapdp*tT[5][ii]+dbpdp*tT[6][ii]+dgpdp*tT[7][ii];
    return deriv;
  }
  Double_t CT(Int_t k, Int_t l,Int_t it) const {return tT[it][l+4*k];}
  AAProd1212bar(double aP, double pP, double Pe, double th) {
    // decay 1/2 -> 1/2 + 0
    Double_t v=TMath::Sqrt(1-aP*aP);
    Double_t gP=v*TMath::Cos(pP);
    Double_t bP=v*TMath::Sin(pP);
    Table(aP,bP,gP,Pe, th);
  }
  AAProd1212bar(double aP, double bP, double gP, double Pe, double th) {
    // decay 1/2 -> 1/2 + 0
    Table(aP,bP,gP,Pe, th);
    iTest=kFALSE;
    
  }
  void Table(double aP, double pP, double Pe, double th) {
    Double_t v=TMath::Sqrt(1-aP*aP);
    Double_t gP=v*TMath::Cos(pP);
    Double_t bP=v*TMath::Sin(pP);
    Table(aP,bP,gP,Pe, th);
    iTest=kFALSE;
  }
  void SetTest(int tst=kTRUE){
    iTest=tst;
  }
  void Table(double aP, double bP, double gP, double Pe, double th) {
    // decay 1/2 -> 1/2 + 0
    fA=aP;
    fP=Pe;
    fB=bP;
    fG=gP;
    Double_t st=TMath::Sin(th);
    Double_t ct=TMath::Cos(th);
    Double_t ct2=ct*ct;
    Double_t st2=1-ct2;
    if(iTest){
      st=1;
      ct=1;
      ct2=1;
      st2=1;
    }
    
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
        if(k==0&&l==0){p0=1;pa=ct2;}
        if(k==0&&l==1)pgp=st;
        if(k==0&&l==2)pb=st*ct;
        if(k==0&&l==3){pp=ct;pap=ct;}
        if(k==1&&l==0)pgp=st;
        if(k==1&&l==1)p0=st2;
        if(k==1&&l==3)pg=st*ct;
        if(k==2&&l==0)pb=-st*ct;
        if(k==2&&l==2)pa=st2;
        if(k==2&&l==3)pbp=-st;
        if(k==3&&l==0){pp=-ct;pap=-ct;}
        if(k==3&&l==1)pg=-st*ct;
        if(k==3&&l==2)pbp=-st;
        if(k==3&&l==3){pa=-1;p0=-ct2;}
        thC[k][l]=p0+aP*pa+bP*pb+gP*pg+aP*Pe*pap+bP*Pe*pbp+gP*Pe*pgp+Pe*pp;
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
  ~AAProd1212bar() {}
private:

  double thC[4][4];
  // individual terms as fcn of angles:
  double tT[8][16];
  Bool_t iTest;
  Double_t fA,fB,fG,fP;
};
#endif /* AA_Prod1212bar_H_ */
