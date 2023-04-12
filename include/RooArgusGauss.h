/*
 * =====================================================================================
 *
 *       Filename:  RooArgusGauss.h
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

#ifndef ROOARGUSGAUSS
#define ROOARGUSGAUSS

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

class RooArgusGauss : public RooAbsPdf {
	public:
		RooArgusGauss() {} ;
		RooArgusGauss(const char *name, const char *title,
				RooAbsReal& _m,
				RooAbsReal& _m0,
				RooAbsReal& _c,
				RooAbsReal& _p,
				RooAbsReal& _mean,
				RooAbsReal& _sigma);
		RooArgusGauss(const RooArgusGauss& other, const char* name=0) ;
		virtual TObject* clone(const char* newname) const { return new RooArgusGauss(*this,newname); }
		inline virtual ~RooArgusGauss() { }

	protected:

		RooRealProxy m ;
		RooRealProxy m0 ;
		RooRealProxy c ;
		RooRealProxy p ;
		RooRealProxy mean ;
		RooRealProxy sigma;

		Double_t evaluate() const ;
};

#endif
