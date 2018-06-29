/*
  Functions and structures used for GSL Monte Carlo integration routines
*/


#include "Thrower.d/GSL_Integration.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"
#include "XSec.d/XSec.h"

#include <cmath>
#include <iostream>

namespace TMDGen {

   double GSL_IntegrationFunction_1h( double *x, size_t dim, void *params ){
      if( dim != 8 )
         throw Error::SanityCheckFailure( "GSL_IntegrationFunction_1h", "dim != 8" );

      GSL_params_t *gsl_params = static_cast< GSL_params_t* >( params );

      Var_t &var( gsl_params->var );

      var.x       = x[0];
      var.y       = x[1];
      var.z       = x[2];
      var.P_hperp = x[3];
      var.phi_h   = x[4];
      var.pT      = x[5];
      var.phi_pT  = x[6];
      var.psi     = x[7];
      var.integrating = 1;

      // remove const-ness to allow it to precompute and store aspects of the state
      return (*const_cast< XSec::XSec_t*>(gsl_params->xsec) )( var );
   };

   double GSL_IntegrationFunction_2h( double *x, size_t dim, void *params ){
      if( dim != 11 )
         throw Error::SanityCheckFailure( "GSL_IntegrationFunction_1h", "dim != 11" );

      GSL_params_t *gsl_params = static_cast< GSL_params_t* >( params );

      Var_t &var( gsl_params->var );

      var.x       = x[0];
      var.y       = x[1];
      var.z       = x[2];
      var.P_hperp = x[3];
      var.phi_h   = x[4];
      var.pT      = x[5];
      var.phi_pT  = x[6];
      var.had_0.M = x[7];
      var.phi_R   = x[8];
      var.cos_vartheta = x[9];
      var.psi     = x[10];
      var.integrating = 1;

      /*
      std::cout << "M_h = " << var.had_0.M << std::endl;
      std::cout << "cos_vartheta = " << var.cos_vartheta << std::endl;
      std::cout << "phi_R = " << var.phi_R << std::endl;

      double output = (*gsl_params->xsec)( var );
      std::cout << "\t output = " << output << std::endl;
      */

      return (*const_cast< XSec::XSec_t*>(gsl_params->xsec) )( var );
   };


};
