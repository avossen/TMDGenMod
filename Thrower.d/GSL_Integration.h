/*
  Functions and structures used for GSL Monte Carlo integration routines
*/

#ifndef _GSL_INTEGRATION_H_
#define _GSL_INTEGRATION_H_

#include "XSec.d/XSec.h"
#include "Common.d/Flavor.h"
#include "Common.d/Var.h"

namespace TMDGen {

   struct GSL_params_t {
      Var_t &var;
      const XSec::XSec_t *xsec;

      GSL_params_t( Var_t& var_in ) : var(var_in) { /* */ };
   };

   double GSL_IntegrationFunction_1h( double *x, size_t dim, void *params );

   double GSL_IntegrationFunction_2h( double *x, size_t dim, void *params );

};

#endif
