/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

// 1/(2 pi)^3
#define ONE_OVER_EIGHT_PI_CUBED 0.004031441804110643

#include "DF.d/BCR08_.h"
#include "DF.d/h1Lperp_BCR08_.h"

#include "Common.d/Consts.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {
   namespace DF {

      h1Lperp_BCR08_t::h1Lperp_BCR08_t() : BCR08_t() {
         constr_msg = "h1Lperp: using BCR08";
      };

      double h1Lperp_BCR08_t::func_s( const Var_t& var ) const {
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double denom ( var.pT * var.pT + L_sq_s( var.x ) );
         denom *= denom;
         denom *= denom;

         double output = PROTON_MASS * M * one_minus_x;
         one_minus_x *= one_minus_x;
         output *= one_minus_x;
         output /= denom;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return -output;
      };

      double h1Lperp_BCR08_t::func_a( const Var_t& var ) const{
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double denom ( var.pT * var.pT + L_sq_a( var.x ) );
         denom *= denom;
         denom *= denom;

         double output = PROTON_MASS * M;
         one_minus_x *= one_minus_x;
         output *= one_minus_x;
         output /= denom;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return output;
      };

      double h1Lperp_BCR08_t::func_aprime( const Var_t& var ) const{
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double denom ( var.pT * var.pT + L_sq_aprime( var.x ) );
         denom *= denom;
         denom *= denom;

         double output = PROTON_MASS * M;
         one_minus_x *= one_minus_x;
         output *= one_minus_x;
         output /= denom;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return output;
      };

   };
};
