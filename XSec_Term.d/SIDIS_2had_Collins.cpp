/*
    Partial wave expansion of SIDIS dihadron Collins term
*/

#include "XSec_Term.d/SIDIS_2had_Collins.h"

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

      // the weights and other prefactors in the pT kT integral

      double SIDIS_2had_Collins_t::func_00( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h + var.phi_S );
      };

      double SIDIS_2had_Collins_t::func_11( const Var_t& var ){
         return 2.*cos( var.phi_pT - var.phi_kT )  * sin( var.phi_R + var.phi_S);
      };
      
      double SIDIS_2had_Collins_t::func_10( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h + var.phi_S );
      };

      double SIDIS_2had_Collins_t::func_1m1( const Var_t& var ){
         return 2.*cos( 2.*var.phi_h - var.phi_pT - var.phi_kT )  * sin( 2.*var.phi_h - var.phi_R + var.phi_S);
      };

      double SIDIS_2had_Collins_t::func_22( const Var_t& var ){
         return 2.*cos( var.phi_h + var.phi_pT - 2.*var.phi_kT ) * sin( -var.phi_h + 2.*var.phi_R + var.phi_S );
      };

      double SIDIS_2had_Collins_t::func_21( const Var_t& var ){
         return 2.*cos( var.phi_pT - var.phi_kT )  * sin( var.phi_R + var.phi_S);
      };

      double SIDIS_2had_Collins_t::func_20( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT ) * sin( var.phi_h + var.phi_S );
      };

      double SIDIS_2had_Collins_t::func_2m1( const Var_t& var ){
         return 2.*cos( 2.*var.phi_h - var.phi_pT - var.phi_kT )  * sin( 2.*var.phi_h - var.phi_R + var.phi_S);
      };

      double SIDIS_2had_Collins_t::func_2m2( const Var_t& var ){
         return 2.*cos( 3.*var.phi_h - var.phi_pT - 2.*var.phi_kT ) * sin( 3.*var.phi_h - 2.*var.phi_R + var.phi_S );
      };

      // constructor
      SIDIS_2had_Collins_t::SIDIS_2had_Collins_t( int l_in, int m_in,
                                            const double& pol_factor_in, const double* y_func, const double* cos_mod_array,
                                            const FlavArrayFuncSet_t* DF_Set, 
                                            const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), dummy_value(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           cos_modulation(&dummy_value), DF_val(&dummy_value), FF_val(&dummy_value), // dummy values for now
           l(l_in), m(m_in),
           func_ptr( &SIDIS_2had_Collins_t::func_00 ) // dummy values for now
      {
         // invalid l,m combinations taken care of by nonzero = FALSE;

         {
            std::stringstream ss;
            ss << "h_1 H_1^perp |" << l << ", " << m << ">";
            message = ss.str();
         };

         //std::cerr << "Collins with l,m = " << l << ' ' << m << std::endl;

         if( DF_Set->Includes("h_1") )
            DF_val = DF_Set->Get_Val_Ptr("h_1");

         overall_sign = -1.;

         if( l == 0 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Collins_t::func_00;
               cos_modulation = &cos_mod_array[Enum::P_00];
               if( FF_Set->Includes("H_1^perp_00") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_00");

            };
         } else if ( l == 1 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Collins_t::func_10;
               cos_modulation = &cos_mod_array[Enum::P_10];
               if( FF_Set->Includes("H_1^perp_10") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_10");

            } else if( m == 1 ){
               func_ptr = &SIDIS_2had_Collins_t::func_11;
               cos_modulation = &cos_mod_array[Enum::P_11];
               if( FF_Set->Includes("H_1^perp_11") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_11");

            } else if( m == -1 ){
               func_ptr = &SIDIS_2had_Collins_t::func_1m1;
               cos_modulation = &cos_mod_array[Enum::P_11];
               if( FF_Set->Includes("H_1^perp_1m1") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_1m1");
            };
         } else if ( l == 2 ){
            if( m == 0 ){
               func_ptr = &SIDIS_2had_Collins_t::func_20;
               cos_modulation = &cos_mod_array[Enum::P_20];
               if( FF_Set->Includes("H_1^perp_20") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_20");

            } else if( m == 1 ){
               func_ptr = &SIDIS_2had_Collins_t::func_21;
               cos_modulation = &cos_mod_array[Enum::P_21];

               if( FF_Set->Includes("H_1^perp_21") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_21");

            } else if( m == -1 ){
               func_ptr = &SIDIS_2had_Collins_t::func_2m1;
               cos_modulation = &cos_mod_array[Enum::P_21];

               if( FF_Set->Includes("H_1^perp_2m1") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_2m1");

            } else if( m == 2 ){
               overall_sign = 1.;

               func_ptr = &SIDIS_2had_Collins_t::func_22;
               cos_modulation = &cos_mod_array[Enum::P_22];

               if( FF_Set->Includes("H_1^perp_22") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_22");

            } else if( m == -2 ){
               func_ptr = &SIDIS_2had_Collins_t::func_2m2;
               cos_modulation = &cos_mod_array[Enum::P_22];

               if( FF_Set->Includes("H_1^perp_2m2") )
                  FF_val = FF_Set->Get_Val_Ptr("H_1^perp_2m2");

            };
         };

         nonzero = ( (FF_val != &dummy_value) && (DF_val != &dummy_value) );

      };


      double SIDIS_2had_Collins_t::Eval( const Var_t& var ){
         //std::cerr << "Collins: " << (*cos_modulation) << std::endl;

         if( 0 && !var.integrating ){
            cout << "- " << flavor_string[ var.flavor ] << ' ';
            cout << overall_sign * var.kT / var.had_0.M << ' ';
            cout << sin( var.phi_pT + var.phi_S + m*( var.phi_R - var.kT) ) << ' ' ;
            cout << (*pol_factor) * (*cos_modulation) * (*y_func_val) << ' ';
            cout << (*DF_val) << ' ' << (*FF_val) << ' ';
            cout << overall_sign * var.kT / var.had_0.M *
               sin( var.phi_pT + var.phi_S + m*( var.phi_R - var.kT) ) *
               (*pol_factor) * (*cos_modulation) * (*y_func_val) * (*DF_val) * (*FF_val) << endl;
         };

         return
            overall_sign * var.kT / var.had_0.M *
            sin( var.phi_pT + var.phi_S + m*( var.phi_R - var.kT) ) *
            // (this->*func_ptr)( var ) *
            (*pol_factor) * (*cos_modulation) * (*y_func_val) * (*DF_val) * (*FF_val);
      };

   };
};
