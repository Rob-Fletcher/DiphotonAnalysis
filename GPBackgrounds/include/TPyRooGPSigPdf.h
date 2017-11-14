
// Author: Wim Lavrijsen   November 2010
// Author: Rob Roy Fletcher October 2017
//  Adapted from TPyFitFunction by Wim Lavrijsen.

#ifndef ROOT_TPyRooGPSigPdf
#define ROOT_TPyRooGPSigPdf

//////////////////////////////////////////////////////////////////////////////
//                                                                          //
// TPyFitFunction                                                           //
//                                                                          //
// Python base class to work with Math::IMultiGenFunction                   //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////


//- ROOT
//#include "Math/IFunction.h"
#include "Rtypes.h"

//- RooFit
#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

// Python
struct _object;
typedef _object PyObject;


class TPyRooGPSigPdf : public RooAbsPdf {
public:
// ctor/dtor, and assignment
   TPyRooGPSigPdf() {};
   TPyRooGPSigPdf( PyObject* self, const char *name, const char *title,
                RooAbsReal& _Myy,
                RooAbsReal& _nSig);
   TPyRooGPSigPdf(PyObject* self, const TPyRooGPSigPdf& other, const char* name); //Not sure if this is required for just getting it working.
   virtual ~TPyRooGPSigPdf();

   virtual TObject* clone(const char* newname) const
      { return new TPyRooGPSigPdf( fPySelf, *this, newname ); }

   ClassDef( TPyRooGPSigPdf, 1 );

protected:
   RooRealProxy Myy;
   RooRealProxy nSig;
   Double_t evaluate() const;

private:
// to prevent confusion when handing 'self' from python
   TPyRooGPSigPdf( const TPyRooGPSigPdf& src ) : RooAbsPdf( src ) {}
   TPyRooGPSigPdf& operator=( const TPyRooGPSigPdf& ) { return *this; }

private:
   PyObject* fPySelf;              //! actual python object
};

#endif
