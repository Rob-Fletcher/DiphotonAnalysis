//Author: Rob Roy Fletcher, Wim Lavrijsen October 2017
// Adapted from TPyFitFunction by Wim Lavrijsen.

// Bindings
#include "PyROOT.h"
#include "TPyRooGPPdf.h"
#include "ObjectProxy.h"
#include "MethodProxy.h"
#include "TPyBufferFactory.h"

//standard
#include <stdexcept>

//______________________________________________________________________________
//                       Python wrapper for Gaussian Process PDF
//                       ================================
//


//- data ---------------------------------------------------------------------
ClassImp(TPyRooGPPdf)

//- helper functions ---------------------------------------------------------
static PyObject* GetOverriddenPyMethod( PyObject* pyself, const char* method )
{
// Retrieve an overriden method on pyself
   PyObject* pymethod = 0;

   if ( pyself && pyself != Py_None ) {
      pymethod = PyObject_GetAttrString( (PyObject*)pyself, const_cast< char* >( method ) );
      if ( ! PyROOT::MethodProxy_CheckExact( pymethod ) )
         return pymethod;

      Py_XDECREF( pymethod );
      pymethod = 0;
   }

   return pymethod;
}

static PyObject* DispatchCall( PyObject* pyself, const char* method,
   PyObject* arg1 = NULL )
{
// Forward <method> to python (need to refactor this with TPySelector).
   PyObject* result = 0;

// get the named method and check for python side overload by not accepting the
// binding's methodproxy
   PyObject* pymethod = GetOverriddenPyMethod( pyself, method );

   if ( pymethod ) {
      result = PyObject_CallFunctionObjArgs( pymethod, arg1, NULL );
   } else {
   // means the method has not been overridden ... simply accept its not there
      result = 0;
      PyErr_Format( PyExc_AttributeError,
         "method %s needs implementing in derived class", const_cast< char* >( method ) );
   }

   Py_XDECREF( pymethod );

   return result;
}

//- constructors/destructor --------------------------------------------------
TPyRooGPPdf::TPyRooGPPdf( PyObject* self,
                          const char *name,
                          const char *title,
                          RooAbsReal& _Myy) :
                          RooAbsPdf(name,title),
                          Myy("Myy","Myy",this,_Myy),
                          fPySelf( 0 )
{
// Construct a TPyRooGPPdf derived with <self> as the underlying
   if ( self ) {
   // steal reference as this is us, as seen from python
      fPySelf = self;
   } else {
      Py_INCREF( Py_None );        // using None allows clearer diagnostics
      fPySelf = Py_None;
   }
}

TPyRooGPPdf::TPyRooGPPdf( PyObject* self,
                          const TPyRooGPPdf& other, const char* name) :
                          RooAbsPdf(other, name),
                          Myy("Myy",this,other.Myy),
                          fPySelf( 0 )
{
// Construct a TPyRooGPPdf derived with <self> as the underlying
   if ( self ) {
   // steal reference as this is us, as seen from python
      fPySelf = self;
   } else {
      Py_INCREF( Py_None );        // using None allows clearer diagnostics
      fPySelf = Py_None;
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Destructor. Only deref if still holding on to Py_None (circular otherwise).

TPyRooGPPdf::~TPyRooGPPdf()
{
   if ( fPySelf == Py_None ) {
      Py_DECREF( fPySelf );
   }
}

////////////////////////////////////////////////////////////////////////////////
/// Simply forward the call to python self.

Double_t TPyRooGPPdf::evaluate() const
{
   PyObject* PyMyy = PyFloat_FromDouble(Myy.arg().getVal());
   PyObject* pyresult = DispatchCall( fPySelf, "evaluate", PyMyy);

   if ( ! pyresult ) {
      PyErr_Print();
      throw std::runtime_error( "Failure in TPyRooGPPdf::evaluate()" );
   }

   double cppresult = (double)PyFloat_AsDouble( pyresult );
   Py_XDECREF( pyresult );

   return cppresult;
}
