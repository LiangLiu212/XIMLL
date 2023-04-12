/*
 * =====================================================================================
 *
 *       Filename:  RooArgusPoly.h
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  05/11/2022 23:48:50
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Liu Liang (LL), liangzy@mail.ustc.edu.cn
 *   Organization:  USTC
 *
 * =====================================================================================
 */

#ifndef ROOARGUSPOLY
#define ROOARGUSPOLY

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooArgusPoly : public RooAbsPdf {
	public:
		RooArgusPoly() {} ;
		RooArgusPoly(const char *name, const char *title,
				RooAbsReal& _m,
				RooAbsReal& _m0,
				RooAbsReal& _c,
				RooAbsReal& _p,
		//		RooAbsReal& _f0,
				RooAbsReal& _f1,
				RooAbsReal& _f2,
				RooAbsReal& _f3);
		RooArgusPoly(const RooArgusPoly& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooArgusPoly(*this,newname); }
		inline virtual ~RooArgusPoly() { }

	protected:

		RooRealProxy m ;
		RooRealProxy m0 ;
		RooRealProxy c ;
		RooRealProxy p ;
	//	RooRealProxy f0;
		RooRealProxy f1;
		RooRealProxy f2;
		RooRealProxy f3;

		Double_t evaluate() const ;
};

#endif
