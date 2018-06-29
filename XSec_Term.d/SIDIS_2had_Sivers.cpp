/*
    Partial wave expansion of SIDIS dihadron Sivers term
*/

#include "XSec_Term.d/SIDIS_2had_Sivers.h"

#include "Common.d/Consts.h"
#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"
#include "Common.d/Yfunctions.h"
#include "Common.d/Enums.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

namespace TMDGen {
   namespace XSec_Term {

      // the weights and other prefactors in the pT kT integral

      double SIDIS_2had_Sivers_t::func_00( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h - var.phi_S );
      };

      double SIDIS_2had_Sivers_t::func_11( const Var_t& var ){
         return 2.*cos( var.phi_pT - var.phi_kT )  * sin( var.phi_R - var.phi_S);
      };
      
      double SIDIS_2had_Sivers_t::func_10( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h - var.phi_S );
      };

      double SIDIS_2had_Sivers_t::func_1m1( const Var_t& var ){
         return 2.*cos( 2.*var.phi_h - var.phi_pT - var.phi_kT )  * sin( 2.*var.phi_h - var.phi_R - var.phi_S);
      };

      double SIDIS_2had_Sivers_t::func_22( const Var_t& var ){
         return 2.*cos( var.phi_h + var.phi_pT - 2.*var.phi_kT ) * sin( -var.phi_h + 2.*var.phi_R - var.phi_S );
      };

      double SIDIS_2had_Sivers_t::func_21( const Var_t& var ){
         return 2.*cos( var.phi_pT - var.phi_kT )  * sin( var.phi_R - var.phi_S);
      };

      double SIDIS_2had_Sivers_t::func_20( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h - var.phi_S );
      };

      double SIDIS_2had_Sivers_t::func_2m1( const Var_t& var ){
         return 2.*cos( 2.*var.phi_h - var.phi_pT - var.phi_kT )  * sin( 2.*var.phi_h - var.phi_R - var.phi_S);
      };

      double SIDIS_2had_Sivers_t::func_2m2( const Var_t& var ){
         return 2.*cos( 3.*var.phi_h - var.phi_pT - 2.*var.phi_kT ) * sin( 3.*var.phi_h - 2.*var.phi_R - var.phi_S );
      };

      // constructor
      SIDIS_2had_Sivers_t::SIDIS_2had_Sivers_t( int l, int m,
                                            const double& pol_factor_in, const double* y_func, const double* cos_mod_array,
                                            const FlavArrayFuncSet_t* DF_Set, 
                                            const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), dummy_value(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           cos_modulation(&dummy_value), Siver_val(&dummy_value), Re_D1_val(&dummy_value), // dummy values for now
           g_1T_val(&dummy_value), Im_D1_val(&dummy_value), // dummy values for now
           func_ptr( &SIDIS_2had_Sivers_t::func_00 ) // dummy values for now
      {
         // invalid l,m combinations taken care of by nonzero = FALSE;

         {
            std::stringstream ss;
            ss << "f_1T^perp D_1 |" << l << ", " << m << ">";
            message = ss.str();
         };

         //std::cerr << "Sivers with l,m = " << l << ' ' << m << std::endl;

         if( DF_Set->Includes("f_1T^perp") )
            Siver_val = DF_Set->Get_Val_Ptr("f_1T^perp");

         if( DF_Set->Includes("g_1T") )
            g_1T_val = DF_Set->Get_Val_Ptr("f_1T^perp");

         D1_sign = 1.;
         overall_sign = -1.;

         if( l == 0 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_00;
               cos_modulation = &cos_mod_array[Enum::P_00];
               if( FF_Set->Includes("Re_D1_00") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_00");

               D1_sign = 0;
            };
         } else if ( l == 1 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_10;
               cos_modulation = &cos_mod_array[Enum::P_10];
               if( FF_Set->Includes("Re_D1_10") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_10");

               D1_sign = 0;

            } else if( m == 1 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_11;
               cos_modulation = &cos_mod_array[Enum::P_11];
               if( FF_Set->Includes("Re_D1_11") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_11");

               if( FF_Set->Includes("Im_D1_11") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_11");

               D1_sign = 1.;

            } else if( m == -1 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_1m1;
               cos_modulation = &cos_mod_array[Enum::P_11];
               if( FF_Set->Includes("Re_D1_11") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_11");

               if( FF_Set->Includes("Im_D1_11") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_11");

               D1_sign = -1.;
            };
         } else if ( l == 2 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_20;
               cos_modulation = &cos_mod_array[Enum::P_20];
               if( FF_Set->Includes("Re_D1_20") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_20");

               D1_sign = 0.;
            } else if( m == 1 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_21;
               cos_modulation = &cos_mod_array[Enum::P_21];

               if( FF_Set->Includes("Re_D1_21") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_21");

               if( FF_Set->Includes("Im_D1_21") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_21");

               D1_sign = 1.;
            } else if( m == -1 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_2m1;
               cos_modulation = &cos_mod_array[Enum::P_21];

               if( FF_Set->Includes("Re_D1_21") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_21");

               if( FF_Set->Includes("Im_D1_21") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_21");

               D1_sign = -1.;
            } else if( m == 2 ){
               overall_sign = 1.;

               func_ptr = &SIDIS_2had_Sivers_t::func_22;
               cos_modulation = &cos_mod_array[Enum::P_22];

               if( FF_Set->Includes("Re_D1_22") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_22");

               if( FF_Set->Includes("Im_D1_22") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_22");

               D1_sign = 1.;
            } else if( m == -2 ){
               func_ptr = &SIDIS_2had_Sivers_t::func_2m2;
               cos_modulation = &cos_mod_array[Enum::P_22];

               if( FF_Set->Includes("Re_D1_22") )
                  Re_D1_val = FF_Set->Get_Val_Ptr("Re_D1_22");

               if( FF_Set->Includes("Im_D1_22") )
                  Im_D1_val = FF_Set->Get_Val_Ptr("Im_D1_22");

               D1_sign = -1.;
            };
         };

         nonzero = ( (Re_D1_val != &dummy_value) && (Siver_val != &dummy_value) );
         nonzero |= ( (Im_D1_val != &dummy_value) && (g_1T_val != &dummy_value) );

      };


      double SIDIS_2had_Sivers_t::Eval( const Var_t& var ){
         double output = (*Siver_val) * (*Re_D1_val);
         if( D1_sign )
            output += D1_sign * (*g_1T_val) * (*Im_D1_val);

         output *= 
            overall_sign * (this->*func_ptr)( var ) *
            var.pT / PROTON_MASS *
            (*pol_factor) * (*cos_modulation) * (*y_func_val);

//            if( !output && !var.integrating ){
//               std::cout << "Sivers " << (*Siver_val) << ' ' << (*Re_D1_val) << ' ' << D1_sign << ' ' << (*g_1T_val) << ' ' <<  (*Im_D1_val) << " | ";
//               std::cout << overall_sign << ' ' << (this->*func_ptr)( var ) << ' ' << var.pT / PROTON_MASS << ' ' << (*pol_factor) << ' ' << (*cos_modulation) << ' ' << (*y_func_val) << std::endl;
//            };

         return output;
      };

   };
};
