
// Author: Wim Lavrijsen   November 2010
// Author: Rob Roy Fletcher October 2017
//  Adapted from TPyFitFunction by Wim Lavrijsen.

#ifndef ROOT_TPyRooGPPdf
#define ROOT_TPyRooGPPdf

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


class TPyRooGPPdf : public RooAbsPdf {
public:
// ctor/dtor, and assignment
   TPyRooGPPdf() {};
   TPyRooGPPdf( PyObject* self, const char *name, const char *title,
                RooAbsReal& _Myy);
   TPyRooGPPdf(PyObject* self, const TPyRooGPPdf& other, const char* name); //Not sure if this is required for just getting it working.
   virtual ~TPyRooGPPdf();

   virtual TObject* clone(const char* newname) const
      { return new TPyRooGPPdf( fPySelf, *this, newname ); }

   ClassDef( TPyRooGPPdf, 1 );

protected:
   RooRealProxy Myy;
   Double_t evaluate() const;

private:
// to prevent confusion when handing 'self' from python
   TPyRooGPPdf( const TPyRooGPPdf& src ) : RooAbsPdf( src ) {}
   TPyRooGPPdf& operator=( const TPyRooGPPdf& ) { return *this; }

private:
   PyObject* fPySelf;              //! actual python object
};

#endif
