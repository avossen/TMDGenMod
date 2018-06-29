/*
  Base class for throwing variables according to some (easy) distribution
  Keeps a pointer to a given random number generator (does not take over ownership)
*/

#ifndef _VARTHROWER_H_
#define _VARTHROWER_H_

#include "Common.d/Var.h"
#include "RNG.d/RNG.h"
#include "XSec.d/XSec.h"

namespace TMDGen {
   namespace Thrower {

      class VarThrower_t {

      protected:
         RNG_t *r;

      public:
         VarThrower_t( RNG_t *r_in ) : r(r_in) { /* */ };
         virtual ~VarThrower_t(){ /* */ };

         virtual int Initialize( Var_t& var, const XSec::XSec_t& xsec ) = 0;
         virtual void Throw( Var_t& var, double& pdf_val ) = 0;
         virtual void ThrowFlavor( Var_t& var, double& pdf_val ) = 0;
         virtual void Throw_pT( Var_t& var, double& pdf_val ) = 0;
      };

   };
};

#endif
