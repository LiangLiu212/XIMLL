/*
 * =====================================================================================
 *
 *       Filename:  RooArgusGauss.cxx
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/2022 23:53:54
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liu Liang (LL), liangzy@mail.ustc.edu.cn
 *   Organization:  USTC
 *
 * =====================================================================================
 */
#include "Riostream.h"
#include "RooArgusGauss.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
//#include "math.h"
#include "TMath.h"
RooArgusGauss::RooArgusGauss(const char *name, const char *title,
		RooAbsReal& _m,
		RooAbsReal& _m0,
		RooAbsReal& _c,
		RooAbsReal& _p,
		RooAbsReal& _mean,
		RooAbsReal& _sigma):
	RooAbsPdf(name,title),
	m("m","m",this,_m),
	m0("m0","m0",this,_m0),
	c("c","c",this,_c),
	p("p","p",this,_p),
	mean("mean","mean",this,_mean),
	sigma("sigma","sigma",this,_sigma){
}


RooArgusGauss::RooArgusGauss(const RooArgusGauss& other, const char* name) :
	RooAbsPdf(other,name),
	m("m",this,other.m),
	m0("m0",this,other.m0),
	c("c",this,other.c),
	p("p",this,other.p),
	mean("mean",this,other.mean),
	sigma("sigma",this,other.sigma){
}



Double_t RooArgusGauss::evaluate() const
{
  Double_t t= m/m0;
  if(t >= 1.0) return 0;
  Double_t u = 1.0 - t*t;
  double argusvalue = ( TMath::Power(2, -p) * TMath::Power(c, 2*(p+1)) / (TMath::Gamma(p+1) - TMath::Gamma(p+1, 0.5 * c*c)) ) * (t/m0) * TMath::Power(u, p) *exp(-0.5 * c*c * u);
  double value  =  exp( -0.5 * (m - mean)*(m - mean) / (sigma*sigma)) / (sigma*TMath::Sqrt(2 * TMath::Pi()));
 //	double value  = mean* m + sigma*m*m;
  //std::cout<< " return value : " << value << " " << m0 << " " << c << " " << p << " " << mean<< " " << sigma << std::endl;
  return value * argusvalue;
}


