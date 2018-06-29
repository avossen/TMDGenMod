/*
    Base class for terms in the cross section
*/

#ifndef _XSEC_TERM_H_
#define _XSEC_TERM_H_

#include "Common.d/Flavor.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {
   namespace XSec_Term {

      class XSec_Term_t {
      protected:
         bool factored, nonzero;
         std::string message;

      public:
         XSec_Term_t( bool b=0, std::string msg="" ) : factored(b), nonzero(0), message(msg) { /* */ };
         XSec_Term_t( std::string msg ) : factored(0), nonzero(0), message(msg) { /* */ };

         virtual ~XSec_Term_t() { /* */ };
         bool Factored() const { return factored; };
         bool NonZero() const { return nonzero; };

         virtual double Eval( const Var_t& var ) = 0;
         virtual std::string Message(){ return message; };

      };
   };
};

#endif
