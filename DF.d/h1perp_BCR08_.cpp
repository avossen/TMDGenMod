/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

// 1/(2 pi)^3
#define ONE_OVER_EIGHT_PI_CUBED 0.004031441804110643

#include "DF.d/BCR08_.h"
#include "DF.d/h1perp_BCR08_.h"

#include "Common.d/Consts.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {
   namespace DF {

      h1perp_BCR08_t::h1perp_BCR08_t() : BCR08_t() {
         constr_msg = "h1perp: using BCR08";
      };

      double h1perp_BCR08_t::func_s( const Var_t& var ) const {

         double pT2 ( var.pT * var.pT );
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double L_sq = L_sq_s( var.x );

         one_minus_x /= ( pT2 + L_sq );
         double output ( one_minus_x );
         one_minus_x *= one_minus_x;
         output *= one_minus_x;           // is now (1-x)^3 / (pT^2 + L_s^2)^3
         output *= M;
         output /= L_sq;
         output *= PROTON_MASS;
         output *= C_F_alpha_S;
         output *= 0.5;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return -output;
      };

      double h1perp_BCR08_t::func_a( const Var_t& var ) const{
         double pT2 ( var.pT * var.pT );
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double L_sq = L_sq_a( var.x );
         double denom = L_sq + pT2;
         double temp = one_minus_x / denom;
         temp *= temp;

         double output( M * PROTON_MASS / L_sq / denom );
         output *= temp;
         output *= C_F_alpha_S;
         output *= 0.5;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return -output;
      };

      double h1perp_BCR08_t::func_aprime( const Var_t& var ) const{
         double pT2 ( var.pT * var.pT );
         double M ( var.x*PROTON_MASS + m );
         double one_minus_x ( 1. - var.x );

         double L_sq = L_sq_aprime( var.x );
         double denom = L_sq + pT2;
         double temp = one_minus_x / denom;
         temp *= temp;

         double output( M * PROTON_MASS / L_sq / denom );
         output *= temp;
         output *= C_F_alpha_S;
         output *= 0.5;
         output *= ONE_OVER_EIGHT_PI_CUBED;

         return -output;
      };

   };
};
