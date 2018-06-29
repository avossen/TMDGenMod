/*
    Partial wave expansion of SIDIS dihadron cos2phi term
*/

#include "XSec_Term.d/SIDIS_2had_cos2phi.h"

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

      double SIDIS_2had_cos2phi_t::func_00( const Var_t& var ){
         return 2.*cos( 2 * var.phi_h - var.phi_pT - var.phi_kT ) * cos( 2.*var.phi_h );
      };

      double SIDIS_2had_cos2phi_t::func_11( const Var_t& var ){
         return 2.*cos( var.phi_h - var.phi_pT )  * cos( var.phi_h + var.phi_R );
      };

      double SIDIS_2had_cos2phi_t::func_10( const Var_t& var ){
         return 2.*cos( 2 * var.phi_h - var.phi_pT - var.phi_kT ) * cos( 2.*var.phi_h );
      };

      double SIDIS_2had_cos2phi_t::func_1m1( const Var_t& var ){
         return 2.*cos( 3*var.phi_h - var.phi_pT - 2.*var.phi_kT ) * cos( 3.*var.phi_h - var.phi_R );
      };

      double SIDIS_2had_cos2phi_t::func_22( const Var_t& var ){
         return 2.*cos( var.phi_pT - var.phi_kT ) * cos( 2.*var.phi_R );
      };

      double SIDIS_2had_cos2phi_t::func_21( const Var_t& var ){
         return 2.* cos( var.phi_h - var.phi_pT )  * cos( var.phi_h + var.phi_R );
         //return cos( (var.phi_pT + var.phi_kT) + (var.phi_R - var.phi_kT) );
      };

      double SIDIS_2had_cos2phi_t::func_20( const Var_t& var ){
         return 2.*cos( 2 * var.phi_h - var.phi_pT - var.phi_kT ) * cos( 2.*var.phi_h );
      };

      double SIDIS_2had_cos2phi_t::func_2m1( const Var_t& var ){
         return 2.*cos( 3.*var.phi_h - var.phi_pT - 2.*var.phi_kT ) * cos( 3.*var.phi_h - var.phi_R );
      };

      double SIDIS_2had_cos2phi_t::func_2m2( const Var_t& var ){
         return 2.*cos( 4.*var.phi_h - var.phi_pT - 3.*var.phi_kT ) * cos( 4.*var.phi_h - 2.*var.phi_R );
      };


      // constructor
      SIDIS_2had_cos2phi_t::SIDIS_2had_cos2phi_t( int l, int m_in,
                                            const double& pol_factor_in, const double* y_func, const double* cos_mod_array,
                                            const FlavArrayFuncSet_t* DF_Set, 
                                            const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           cos_modulation(&pol_factor_in), DF_val(&pol_factor_in), FF_val(&pol_factor_in), // dummy values for now
           func_ptr( &SIDIS_2had_cos2phi_t::func_00 ),  m(m_in)
      {
         // invalid l,m combinations taken care of by nonzero = FALSE;

         //std::cerr << "cos2phi with l,m = " << l << ' ' << m << std::endl;

         {
            std::stringstream ss;
            ss << "h_1^perp H_1^perp |" << l << ", " << m << ">";
            message = ss.str();
         };

         if( DF_Set->Includes("h_1^perp") ){
            DF_val = DF_Set->Get_Val_Ptr("h_1^perp");

            if( l == 0 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_00;
                  cos_modulation = &cos_mod_array[Enum::P_00];
                  if( FF_Set->Includes("H_1^perp_00") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_00");
                     nonzero = 1;
                  };
               };
            } else if ( l == 1 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_10;
                  cos_modulation = &cos_mod_array[Enum::P_10];
                  if( FF_Set->Includes("H_1^perp_10") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_10");
                     nonzero = 1;
                  };
               } else if( m == 1 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_11;
                  cos_modulation = &cos_mod_array[Enum::P_11];
                  if( FF_Set->Includes("H_1^perp_11") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_11");
                     nonzero = 1;
                  };
               } else if( m == -1 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_1m1;
                  cos_modulation = &cos_mod_array[Enum::P_11];
                  if( FF_Set->Includes("H_1^perp_1m1") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_1m1");
                     nonzero = 1;
                  };
               };
            } else if ( l == 2 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_20;
                  cos_modulation = &cos_mod_array[Enum::P_20];
                  if( FF_Set->Includes("H_1^perp_20") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_20");
                     nonzero = 1;
                  };
               } else if( m == 1 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_21;
                  cos_modulation = &cos_mod_array[Enum::P_21];
                  if( FF_Set->Includes("H_1^perp_21") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_21");
                     nonzero = 1;
                  };
               } else if( m == -1 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_2m1;
                  cos_modulation = &cos_mod_array[Enum::P_21];
                  if( FF_Set->Includes("H_1^perp_2m1") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_2m1");
                     nonzero = 1;
                  };
               } else if( m == 2 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_22;
                  cos_modulation = &cos_mod_array[Enum::P_22];
                  if( FF_Set->Includes("H_1^perp_22") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_22");
                     nonzero = 1;
                  };
               } else if( m == -2 ){
                  func_ptr = &SIDIS_2had_cos2phi_t::func_2m2;
                  cos_modulation = &cos_mod_array[Enum::P_22];
                  if( FF_Set->Includes("H_1^perp_2m2") ){
                     FF_val = FF_Set->Get_Val_Ptr("H_1^perp_2m2");
                     nonzero = 1;
                  };
               };
            };
         };
      };

      double SIDIS_2had_cos2phi_t::Eval( const Var_t& var ){
         //std::cerr << " SIDIS_2had_cos2phi_t::Eval() " << (*cos_modulation) << std::endl;
         //          cout << "cos2phi\t" << cos( var.phi_kT + var.phi_pT + m*(var.phi_R - var.phi_kT) ) << ' ';
         //          cout << var.kT << ' ' << (-var.pT)*PROTON_MASS/var.had_0.M << ' ' ;
         //          cout << (*pol_factor) << ' ' << (*cos_modulation) << ' ' << (*y_func_val) << ' ' << (*DF_val) << ' ' << (*FF_val);
         //          cout << endl;

//          if( !var.integrating )
//             cout << "b) " << (*DF_val) << endl;

         return
            //(this->*func_ptr)( var ) *
            cos( var.phi_kT + var.phi_pT + m*(var.phi_R - var.phi_kT) ) *
            (-var.pT) *var.kT/PROTON_MASS/var.had_0.M *
            (*pol_factor) * (*cos_modulation) * (*y_func_val) * 
            (*DF_val) *
            (*FF_val);

         //return cos( var.phi_kT + var.phi_pT + m*(var.phi_R - var.phi_kT) ) * (*cos_modulation) * (*FF_val);
      };

   };
};






