/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

// 1/(2 pi)^3
#define ONE_OVER_EIGHT_PI_CUBED 0.004031441804110643

#include "DF.d/BCR08_.h"
#include "DF.d/h1Tperp_BCR08_.h"

#include "Common.d/Consts.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {
   namespace DF {

      h1Tperp_BCR08_t::h1Tperp_BCR08_t() : BCR08_t() {
         constr_msg = "h1Tperp: using BCR08";
      };

      double h1Tperp_BCR08_t::func_s( const Var_t& var ) const {
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         M *= one_minus_x;
         M *= M;
         M *= one_minus_x;

         double denom ( var.pT * var.pT + L_sq_s( var.x ) );
         denom *= denom;
         denom *= denom;

         M /= denom;
         M *= ONE_OVER_EIGHT_PI_CUBED;

         return M;
      };

      double h1Tperp_BCR08_t::func_a( const Var_t& var ) const{
         return 0;
      };

      double h1Tperp_BCR08_t::func_aprime( const Var_t& var ) const{
         return 0;
      };

   };
};
