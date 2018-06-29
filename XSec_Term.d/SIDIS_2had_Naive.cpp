/*
    Naive model for dihadron terms
*/

#include "XSec_Term.d/SIDIS_2had_Naive.h"

#include "Common.d/Consts.h"
#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"
#include "Common.d/Yfunctions.h"
#include "Common.d/Enums.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
using std::cout;
using std::cerr;
using std::endl;


namespace TMDGen {
   namespace XSec_Term {

      SIDIS_2had_Naive_t::params_t::params_t() :
         A(0),
         x_slope(0),
         z_alpha(0),
         z_beta(0),
         pT_slope(0),
         kT_slope_alpha(0),
         kT_slope_beta(0),
         kT_slope_gamma(0),
         M_h_a1(0),
         M_h_a2(0)
      {
         /* */
      };

      // constructor
      SIDIS_2had_Naive_t::SIDIS_2had_Naive_t( std::string& pol_state_in, int mom_num_in, std::string& params_in,
                                              const double& pol_factor_in,
                                              const double* y_func_array,
                                              const double* cos_mod_array )
         : XSec_Term_t(1),
           pol_factor( pol_factor_in ),
           y_func_val( y_func_array ),
           y_func_val_A( &y_func_array[yFunctions::A_idx] ),
           cos_modulation(&pol_factor_in),
           mom_num(mom_num_in),
           pol_state(pol_state_in)
      {
         // copy out parameters
         {
            std::stringstream ss;
            ss << params_in;

            ss >> params.A >> params.x_slope >> params.z_alpha >> params.z_beta >> params.pT_slope;
            ss >> params.kT_slope_alpha >> params.kT_slope_beta >> params.kT_slope_gamma >> params.M_h_a1 >> params.M_h_a2;
         };

         std::string type = "INVALID";

         if( pol_state_in == "UU" ){
            is_cos = 1;
            n_S = 0;
            if( mom_num < 1 ){
               // invalid
            } else if( mom_num < 6 ){
               type = "f_1 D_1";
               y_func_val = &( y_func_array[yFunctions::A_idx] );
               n_h = 0;

               if( mom_num == 1 ){
                  l = 1;
                  m = 0;
               } else if ( mom_num == 2 ){
                  l = 1;
                  m = 1;
               } else {
                  l = 2;
                  m = mom_num - 3;
               };
            } else if( mom_num < 24 ){
               if( mom_num < 15 ){
                  mom_num -= 6;
                  type = "h_1^perp H_1^perp";
                  y_func_val = &( y_func_array[yFunctions::B_idx] );
                  n_h = 2;
               } else {
                  mom_num -= 15;
                  type = "'cos(phi)'";
                  y_func_val = &( y_func_array[yFunctions::V_idx] );
                  n_h = 1;
               };
               if( mom_num == 0 ){
                  l = 0;
                  m = 0;
               } else if (mom_num < 4){
                  l = 1;
                  m = mom_num - 2;
               } else {
                  l = 2;
                  m = mom_num - 6;
               };
            };
         } else if ( pol_state_in == "UT" ){
            is_cos = 0;

            // just leading twist so far
            if( mom_num > -1 && mom_num < 27 ){

               if( mom_num < 9 ){
                  n_h = 1;
                  n_S = -1;
                  type = "f_1T^perp D_1";
                  y_func_val = &( y_func_array[yFunctions::A_idx] );
               } else if( mom_num < 18 ){
                  n_h = 1;
                  n_S = 1;
                  type = "h_1 H_1^perp";
                  y_func_val = &( y_func_array[yFunctions::B_idx] );
               } else {
                  n_h = 3;
                  n_S = -1;
                  type = "h_1T^perp H_1^perp";
                  y_func_val = &( y_func_array[yFunctions::B_idx] );
               };

               mom_num %= 9;
               if( mom_num == 0 ){
                  l = 0;
                  m = 0;
               } else if (mom_num < 4){
                  l = 1;
                  m = mom_num - 2;
               } else {
                  l = 2;
                  m = mom_num - 6;
               };
            };
         };

         if( type == "INVALID" ){
            message = type;
            return;
         };

         // set message
         {
            std::stringstream ss;
            ss << type << " |" << l << ", " << m << "> NAIVE";
            message = ss.str();
         };

         // set cos modulation
         if( l == 0 ){
            if( m == 0 ){
               cos_modulation = &cos_mod_array[Enum::P_00];
               nonzero = 1;
            };
         } else if ( l == 1 ){
            if( m == -1 || m == 1 ){
               cos_modulation = &cos_mod_array[Enum::P_11];
               nonzero = 1;
            } else if ( m == 0 ){
               cos_modulation = &cos_mod_array[Enum::P_10];
               nonzero = 1;
            };
         } else if ( l == 2 ){
            if( m == -2 || m == 2 ){
               cos_modulation = &cos_mod_array[Enum::P_22];
               nonzero = 1;
            } else if ( m == -1 || m == 1 ){
               cos_modulation = &cos_mod_array[Enum::P_21];
               nonzero = 1;
            } else if ( m == 0 ){
               cos_modulation = &cos_mod_array[Enum::P_20];
               nonzero = 1;
            };
         };

      };

      double SIDIS_2had_Naive_t::Eval( const Var_t& var ){
         double A = params.A;

         if( params.x_slope )
            A *= ( 1 + params.x_slope*log(var.x) );

         if( params.z_alpha )
            A *= pow( var.z, params.z_alpha );

         if( params.z_beta )
            A *= pow( var.z, 1 - params.z_beta );

         if( params.pT_slope ){
            double temp = var.pT;

            temp *= temp;
            temp /= params.pT_slope;
            temp = exp( -temp );
            temp /= (params.pT_slope*PI);

            A *= temp;
         };

         if( params.kT_slope_alpha ){
            double temp = params.kT_slope_alpha*pow( var.z, params.kT_slope_beta )*pow( 1.-var.z, params.kT_slope_gamma );

            if( temp ){
               double val = var.kT;
               val *= var.z/(1.-var.z);
               val /= temp;
               val *= val;
               val = exp( -2.*val );

               A *= val;
            };
         };

         if( params.M_h_a1 ){
            double temp = 1;
            temp += params.M_h_a1 * var.had_0.M;
            temp += params.M_h_a2 * var.had_0.M * var.had_0.M;

            A *= temp;
         };

         double arg = (n_h-m)*var.phi_h + n_S*var.phi_S + m*var.phi_R;

         A *= ( is_cos ? cos(arg) : sin(arg) );
         A /= (*y_func_val_A);
         A *= PI;

         return A * (pol_factor) * (*cos_modulation) * (*y_func_val);
      };

   };
};






