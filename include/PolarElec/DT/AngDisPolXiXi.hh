#ifndef MN_AngDisPolXiXi_H_
#define MN_AngDisPolXiXi_H_
#include "AAPolProd1212.hh"
#include "AADecay12.hh"
#include <cmath>
#include "TMath.h"
/*
  AngDis  - for e+e- -> XiXi exclusive distributions including
            decay chains
  2018-06-09 : v2.00 (AK) 
    Fast calculation of normalization:
    Functions : 
      To prepare integrals:
        InitIntegral() - determine # of non zero terms
  AddToIntegral(...)  add an MC PHSP event to integral
        PrintIntegral() - print calculated integrals
  SaveIntegral(fname) - save integrals to text file
      or read them from prepared file
        ReadIntegral(fname)
    ---------
    Usage: create the class once per run and then
      -- for each change of parameters:
      Set(...) sets current values of decay parameters
      -- for each event:
        W(...) calculates event weight for the specified angles
               using current decay parameters
      CalcIntegral() calculate normalization factor using curent 
                decay parameters
      
 */
using namespace std;

class AngDis {

  public:

    AngDis(double pl, double ap, double pp, double aXi, double pXi, double aL, double aXib, double aLb, double pXib) :
  pol(pl),tHaPsi(ap),tHpPsi(pp),tHaXi(aXi), tHpXi(pXi), tHaL(aL), tHaXib(aXib), tHaLb(aLb), tHpXib(pXib)
  {
    prod  = new AAPolProd1212(tHaPsi,tHpPsi,pol,0);
    dXi   = new AADecay12(tHaXi,tHpXi,0,0);
    dXib  = new AADecay12(tHaXib,tHpXib,0,0);
    dL    = new AADecay12(tHaL,0,0,0);
    dLb   = new AADecay12(tHaLb,0,0,0);
    iTest=0;
    iNZero=0;
  } 
  void Set(double pl, double ap, double pp, double aXi, double pXi, double aL, double aXib, double aLb, double pXib)
  {
    // Set parameters
    pol=pl;
    tHaPsi=ap;
    tHpPsi=pp;
    tHaXi=aXi;
    tHpXi=pXi;
    tHaL=aL;
    tHaXib=aXib;
    tHaLb=aLb;
    tHpXib=pXib;
    prod->Table(tHaPsi,tHpPsi,pol,0);
    dXi->Table(tHaXi,tHpXi,0,0);
    dL->Table(tHaL,0,0,0);
    dXib->Table(tHaXib,tHpXib,0,0);
    dLb->Table(tHaLb,0,0,0);
  }     
  ~AngDis() 
  {
    delete prod;
    delete dXi;
    delete dL;
    delete dXib;
    delete dLb;
    
  }
  void SaveIntegral(const char *fname)
  {
    // Save calculated integrals to file fname

   fstream inTg(fname,ios::out);
    inTg<<iNZero<<endl;
    for(int i=0;i<iNZero;i++){
      inTg<<i<<" "<<gID[i]<<" "<<fInt[i]<<endl;
    }
    inTg.close();
  }
  void ReadIntegral(const char *fname)
  {
    // Read calculated integrals from file fname
    fstream inTg(fname,ios::in);
    inTg>>iNZero;
    for(int i=0;i<iNZero;i++){
      int id;
      inTg>>id>>gID[i]>>fInt[i];
    }
    fN=fInt[0];
    inTg.close();
  }
  void MergeIntegral(const char *fname)
  {
    // Read calculated integrals from file fname
    // 2018-11-22 (AK) and add to the existing in memory 
    //               To be used for thread merging
    Int_t tmiNZero,tmgID;
    Double_t tmfInt;
    fstream inTg(fname,ios::in);
    inTg>>tmiNZero;
    if(tmiNZero!=iNZero)cout<<"Wrong number of terms!!!\n"; 
    for(int i=0;i<iNZero;i++){
      int id;
      //      inTg>>id>>gID[i]>>fInt[i];
      inTg>>id>>tmgID>>tmfInt;
      if(tmiNZero!=iNZero)cout<<"Wrong id!!!\n";
      fInt[i]+=tmfInt;
    }
    fN=fInt[0];
    inTg.close();
  }
  void InitIntegral()
  {
    // Clear matrix with integrals
    for(int i=0;i<387420489;i++){fInt[i]=0;gDI[i]=0;}
    fN=0;
    // determine which terms are !=0
    // generally there could be 4^5=1024 terms
    //  but it looks as there are only 72 !=0
    //  
    prod->SetTest();
    dXi->SetTest();
    dL->SetTest();
    dXib->SetTest();
    dLb->SetTest();
    iTest=1;
    iNZero=1024;
    AddToIntegral(0,0,0,0,
                   0,0,0,0,
      0);
    iNZero=0;
    for(int i=0;i<387420489;i++){
      if(!(fInt[i]<1)){
        gID[iNZero]=i;
        gDI[i]=iNZero;
        iNZero++;
      //  cout<<" El : "<<i<<" "<<iNZero<<" "<<fInt[i]<<endl;
      }
    }
   // cout<<"Nzero terms : "<<iNZero<<" "<<endl;
    for(int i=0;i<iNZero;i++){
      fInt[i]=0;
    }
    iTest=0;
    fN=0;
    prod->SetTest(0);
    dXi->SetTest(0);
    dL->SetTest(0);
    dXib->SetTest(0);
    dLb->SetTest(0);
  }
  void PrintIndFast()
  {
    for(int ii=0;ii<iT;ii++){ //i1 0 prod
      Int_t i=iI[ii];
      Int_t id=gID[i];
      Int_t k1=iK1[ii];
      Int_t k2=iK2[ii];
      Int_t j1=iJ1[ii];
      Int_t j2=iJ2[ii];
    //  cout<<ii<<" "<<i<<" "<<k1<<" "<<k2<<" "<<j1<<" "<<j2<<endl;
    }
  }
  void AddToIntegral(double th,double th1L,double ph1L,double th1p,
                   double ph1p,double th2L,double ph2L,double th2p,
           double ph2p)
  {
    // add one PHSP event to the sums of the terms of linear expansion
    prod->Table(tHaPsi,tHpPsi,pol,th);
    dXi->Table(tHaXi,tHpXi,th1L,ph1L);
    dL->Table(tHaL,0,th1p,ph1p);
    dXib->Table(tHaXib,tHpXib,th2L,ph2L);
    dLb->Table(tHaLb,0,th2p,ph2p);
    fN++;
    if(iTest)iT=0;
    for(int i=0;i<iNZero;i++){ //i1 0 prod
                              // i2 1 Loop par dXi
                        // i3 2 Loop par dL
                              // i4 3 Loop par dXib
                        // i5 4 Loop par dLb
      //        Int_t id=i5*256+i4*64+i3*16+i2*4+i1;
      Int_t id=gID[i];
      if(iTest)id=i;
      // Decode id of each sub-process
      Int_t i1=Id(id,0);
      Int_t i2=Id(id,1);
      Int_t i3=Id(id,2);
      Int_t i4=Id(id,3);
      Int_t i5=Id(id,4);
      double tep=0;
      for(int k1=0;k1<4;k1++){// Xi loop
         for(int k2=0;k2<4;k2++){// Xi_bar loop
           for(int j1=0;j1<4;j1++){// L loop
              for(int j2=0;j2<4;j2++){// L_bar loop

                Double_t WW= prod->CT(k1,k2,i1);
                WW*=dXi->AT(k1,j1,i2)*
                dXib->AT(k2,j2,i4)*
                dL->AT(j1,0,i3)*dLb->AT(j2,0,i5);
    if(iTest){
      WW=TMath::Abs(WW);
      if(WW>0){
        //cout<<iT<<" "<<i<<" "<<k1<<" "<<k2<<" "<<j1<<" "<<j2<<endl;
        iI[iT]=i;
        iK1[iT]=k1;
        iK2[iT]=k2;
        iJ1[iT]=j1;
        iJ2[iT]=j2;
        iT++;
        
      }
    }
                tep+=WW;
              }
           }
         }
      }
      fInt[i]+=tep;// calc all terms
    }
  }
  void AddToIntegralFast(double th,double th1L,double ph1L,double th1p,
                   double ph1p,double th2L,double ph2L,double th2p,
           double ph2p)
  {
    // add one PHSP event to the sums of the terms of linear expansion
    // 2018-11-22 Fast version
    prod->Table(tHaPsi,tHpPsi,pol,th);
    dXi->Table(tHaXi,tHpXi,th1L,ph1L);
    dL->Table(tHaL,0,th1p,ph1p);
    dXib->Table(tHaXib,tHpXib,th2L,ph2L);
    dLb->Table(tHaLb,0,th2p,ph2p);
    fN++;
    for(int ii=0;ii<iT;ii++){ //i1 0 prod
                              // i2 1 Loop par dXi
                        // i3 2 Loop par dL
                              // i4 3 Loop par dXib
                        // i5 4 Loop par dLb
      //        Int_t id=i5*256+i4*64+i3*16+i2*4+i1;
      Int_t id=iI[ii];
      //      Int_t id=gID[i];
      Int_t k1=iK1[ii];
      Int_t k2=iK2[ii];
      Int_t j1=iJ1[ii];
      Int_t j2=iJ2[ii];
      // Decode id of each sub-process
      Int_t i1=Id(id,0);
      Int_t i2=Id(id,1);
      Int_t i3=Id(id,2);
      Int_t i4=Id(id,3);
      Int_t i5=Id(id,4);
      Int_t i=gDI[id];

      double tep=0;
      
      Double_t WW= prod->CT(k1,k2,i1);
      WW *= dXi->AT(k1,j1,i2) * dXib->AT(k2,j2,i4)*
            dL->AT(j1,0,i3) * dLb->AT(j2,0,i5);
      fInt[i]+=WW;// calc all terms
    }
  }
 void PrintIntegral()
  {
    // Print calculated integrals 
    for(int i=0;i<iNZero;i++){
    //  cout<<" El : "<<i<<" "<<gID[i]<<" "<<fInt[i]<<endl;
    }
    //cout<<"Nzero terms : "<<iNZero<<" "<<endl;
    
  }
 Double_t CalcIntegral()
  {
    // Calculate normalization from stored sums of the terms
    Double_t p1[5];
    Double_t d1[4];
    Double_t d2[4];
    Double_t d3[4];
    Double_t d4[4];
    Double_t vvP=TMath::Sqrt(1-tHaPsi*tHaPsi);
    Double_t vvX=TMath::Sqrt(1-tHaXi*tHaXi);
    Double_t vvXb=TMath::Sqrt(1-tHaXib*tHaXib);
    p1[0]=1;p1[1]=tHaPsi;p1[2]=vvP*TMath::Sin(tHpPsi);p1[3]=vvP*TMath::Cos(tHpPsi); p1[4]=pol;
    d1[0]=1;d1[1]=tHaXi;d1[2]=vvX*TMath::Sin(tHpXi);d1[3]=vvX*TMath::Cos(tHpXi);
    d2[0]=1;d2[1]=tHaL;d2[2]=0;d2[3]=0;
    d3[0]=1;d3[1]=tHaXib;d3[2]=vvXb*TMath::Sin(tHpXib);
    d3[3]=vvXb*TMath::Cos(tHpXib);
    d4[0]=1;d4[1]=tHaLb;d4[2]=0;d4[3]=0;
    double tep=0;
    for(int i=0;i<iNZero;i++){ // Loop over non zero integrals
      Int_t id=gID[i];
      Int_t i1=Id(id,0);
      Int_t i2=Id(id,1);
      Int_t i3=Id(id,2);
      Int_t i4=Id(id,3);
      Int_t i5=Id(id,4);
      tep+=p1[i1]*d1[i2]*d2[i3]*d3[i4]*d4[i5]*fInt[i];
    }
    return tep/fN;
  }
    double  W(double th,double th1L,double ph1L,double th1p,
           double ph1p,double th2L,double ph2L,double th2p,
           double ph2p) 

    {
      // Calculate weight for event specified by the variables 
      prod->Table(tHaPsi,tHpPsi,pol,th);
      dXi->Table(tHaXi,tHpXi,th1L,ph1L);
      dXib->Table(tHaXib,tHpXib,th2L,ph2L);
      dL->Table(tHaL,0,th1p,ph1p);
      dLb->Table(tHaLb,0,th2p,ph2p);
      double tep=0;
      for(int i1=0;i1<4;i1++){// Xi loop
        for(int i2=0;i2<4;i2++){// Xi_bar loop
          for(int j1=0;j1<4;j1++){// L loop
            for(int j2=0;j2<4;j2++){// L_bar loop
              tep += prod->C(i1,i2)*dXi->A(i1,j1)*dXib->A(i2,j2)*
              dL->A(j1,0)*dLb->A(j2,0);
            }
          }
        }
      }
      return tep;
    }



private:
  
double pol;
double tHaPsi;
double tHpPsi;
double tHaXi;
double tHpXi;
double tHaL;
double tHaXib;
double tHpXib;
double tHaLb;
  // integral variables:
  Double_t fInt[387420489]; Double_t fN;
  AAPolProd1212 *prod;
  AADecay12 *dXi;
  AADecay12 *dL;
  AADecay12 *dXib;
  AADecay12 *dLb;
  // Test flag to determine # of non zero terms
  Bool_t iTest;
  // Non zero terms and indicies:
  Int_t iNZero;
  Int_t gID[387420489];
  Int_t gDI[387420489];
  Int_t iI[387420489];
  Int_t iK1[387420489];
  Int_t iK2[387420489];
  Int_t iJ1[387420489];
  Int_t iJ2[387420489];
  Int_t iT;
  Int_t Id(UInt_t ig,UInt_t ip)
  {
    // Decode global 0<=ig<4^5 
    // 0<=id<4  <- id of term for each process  0<=ip<5
    UInt_t k=2*ip;
    int ii=0x3&(ig>>k);
    return ii;
  }
};

#endif // MN_AngDisXiXi_H_

