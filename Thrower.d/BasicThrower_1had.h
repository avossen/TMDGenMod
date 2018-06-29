/*
     Throws variables uniformly per flavor
     Most similar to old FORTRAN w/ weights option
     For single hadrons
*/

#ifndef _BASICTHROWER_1HAD_H_
#define _BASICTHROWER_1HAD_H_

#include "Common.d/Var.h"
#include "RNG.d/RNG.h"
#include "XSec.d/XSec.h"

#include "Thrower.d/BasicThrower.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_vegas.h>

namespace TMDGen {
   namespace Thrower {

      class BasicThrower_1had_t : public BasicThrower_t {

         int N_integration_warmup;
         int N_integration_calls;
         int N_max_integration_calls;

         double Vinv;

         // specialization for different dependent variables
         virtual void Init_GSL_Func( gsl_monte_function *gsl_F, double* & min, double* &max );

         // specialization for different dependent variables and
         // for throwing according to different distributions
         virtual void ThrowVariables( Var_t& var, double& pdf_val );

      public:
         BasicThrower_1had_t( RNG_t *r_in, int N_integration_warmup_, int N_integration_calls,
                         int N_max_integration_calls_, const Var_t& min, const Var_t& max );
         virtual ~BasicThrower_1had_t() { /* */ };

      };

   };
};

#endif
