/*
 * =====================================================================================
 *
 *       Filename:  RooArgusPoly.cxx
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
#include "RooArgusPoly.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
//#include "math.h"
#include "TMath.h"
RooArgusPoly::RooArgusPoly(const char *name, const char *title,
		RooAbsReal& _m,
		RooAbsReal& _m0,
		RooAbsReal& _c,
		RooAbsReal& _p,
	//	RooAbsReal& _f0,
		RooAbsReal& _f1,
		RooAbsReal& _f2,
		RooAbsReal& _f3):
	RooAbsPdf(name,title),
	m("m","m",this,_m),
	m0("m0","m0",this,_m0),
	c("c","c",this,_c),
	p("p","p",this,_p),
//	f0("f0","f0",this,_f0),
	f1("f1","f1",this,_f1),
	f2("f2","f2",this,_f2),
	f3("f3","f3",this,_f3){
}


RooArgusPoly::RooArgusPoly(const RooArgusPoly& other, const char* name) :
	RooAbsPdf(other,name),
	m("m",this,other.m),
	m0("m0",this,other.m0),
	c("c",this,other.c),
	p("p",this,other.p),
//	f0("f0",this,other.f0),
	f1("f1",this,other.f1),
	f2("f2",this,other.f2),
	f3("f3",this,other.f3){
}



Double_t RooArgusPoly::evaluate() const
{
  Double_t t= m/m0;
  if(t >= 1.0) return 0;
  Double_t u = 1.0 - t*t;
  double argusvalue = ( TMath::Power(2, -p) * TMath::Power(c, 2*(p+1)) / (TMath::Gamma(p+1) - TMath::Gamma(p+1, 0.5 * c*c)) ) * (t/m0) * TMath::Power(u, p) *exp(-0.5 * c*c * u);
//  double argusvalue = m*TMath::Power(u,p)*exp(c*u);
 // double value  =  exp( -0.5 * (m - mean)*(m - mean) / (sigma*sigma)) / (sigma*TMath::Sqrt(2 * TMath::Pi()));
 	double value  = f1*m + f2*m*m + f3*m*m*m + 1.0;
 // std::cout<< " return value : " << value * argusvalue << "  " << argusvalue << "  " <<  value << "  "  << m <<  "  " << t <<  "  " << u <<  "  " << t/m0 <<  "  " << (TMath::Gamma(p+1) - TMath::Gamma(p+1, 0.5 * c*c)) << "  " << TMath::Power(2, -p) << "  " << TMath::Power(c, 2*(p+1)) << "  " << TMath::Power(u, p) << "  " <<  exp(-0.5 * c*c * u);
 // std::cout<< " return value : " << value * argusvalue << "  " << argusvalue << "  " <<  value << "  " << c << "  " << p;
 // std::cout << std::endl;
  return value * argusvalue;
}


