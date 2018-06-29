/*
     Throws variables uniformly per flavor
     Most similar to old FORTRAN w/ weights option
*/

#ifndef _BASICTHROWER_H_
#define _BASICTHROWER_H_

#include "Common.d/Flavor.h"
#include "Common.d/Var.h"
#include "RNG.d/RNG.h"
#include "XSec.d/XSec.h"
#include "Thrower.d/VarThrower.h"

#include <gsl/gsl_monte.h>

namespace TMDGen {
   namespace Thrower {

      class BasicThrower_t : public VarThrower_t {

         int N_integration_warmup;
         int N_integration_calls;
         int N_max_integration_calls;

         double integrated_xsec[ GMC_TRANS_N_FLAVORS ];
         double integrated_xsec_abserr[ GMC_TRANS_N_FLAVORS ];
         double total_xsec;

         // specialization for different dependent variables
         virtual void Init_GSL_Func( gsl_monte_function*, double* &min, double* &max ) = 0;

         // specialization for different dependent variables and
         // for throwing according to different distributions
         virtual void ThrowVariables( Var_t& var, double& pdf_val ) = 0;

      protected:
         const Var_t &min, &max;
         Var_t width;

      public:
         BasicThrower_t( RNG_t *r_in, int N_integration_warmup_, int N_integration_calls,
                         int N_max_integration_calls_, const Var_t& min_in, const Var_t& max_in );
         virtual ~BasicThrower_t();

         virtual int Initialize(  Var_t& var, const XSec::XSec_t& xsec );
         virtual void Throw( Var_t& var, double& pdf_val );
         virtual void ThrowFlavor( Var_t& var, double& pdf_val );
         virtual void Throw_pT( Var_t& var, double& pdf_val );
      };

   };
};

#endif
