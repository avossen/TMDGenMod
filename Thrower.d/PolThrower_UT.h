/*
  Fully unpolarized (trivial case)
*/

#ifndef _POLTHROWER_UT_H_
#define _POLTHROWER_UT_H_

#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"

#include "RNG.d/RNG.h"

#include <gsl/gsl_rng.h>

namespace TMDGen {
   namespace Thrower {

      class PolThrower_UT_t : public PolThrower_t { 

         RNG_t *r;
         double P_T;

      public:
         PolThrower_UT_t( RNG_t* r_, double P_T_ ) : r(r_), P_T(P_T_) {
            if( P_T <= 0 || P_T > 1 )
               throw Error::Constructing("PolThrower_UT_t", "polarization must be in (0,1]" );
         };

         virtual ~PolThrower_UT_t(){ /* */ };

         virtual void Throw( Var_t& var ){
            var.P_L = 0;
            var.P_T = (r->EvalUnif() > 0.5 ) ? P_T : -P_T ;
         };
      };

   };
};

#endif
