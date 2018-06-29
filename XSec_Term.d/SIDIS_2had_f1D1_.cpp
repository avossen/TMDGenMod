/*
    Partial wave expansion of SIDIS dihadron f1D1 term
*/

#include "XSec_Term.d/SIDIS_2had_f1D1_.h"

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

      double SIDIS_2had_f1D1_t::func_00( const Var_t& var ){
         return 1;
      };

      double SIDIS_2had_f1D1_t::func_11( const Var_t& var ){
         return 2*cos( var.phi_h - var.phi_kT ) * cos( var.phi_h - var.phi_R );
      };

      double SIDIS_2had_f1D1_t::func_10( const Var_t& var ){
         return 1;
      };

      double SIDIS_2had_f1D1_t::func_22( const Var_t& var ){
         return 2*cos( 2.*var.phi_h - 2.*var.phi_kT ) * cos( 2.*var.phi_h - 2.*var.phi_R );
      };

      double SIDIS_2had_f1D1_t::func_21( const Var_t& var ){
         return 2*cos( var.phi_h - var.phi_kT ) * cos( var.phi_h - var.phi_R );
      };

      double SIDIS_2had_f1D1_t::func_20( const Var_t& var ){
         return 1;
      };


      // constructor
      SIDIS_2had_f1D1_t::SIDIS_2had_f1D1_t( int l_in, int m_in,
                                            const double& pol_factor_in, const double* y_func, const double* cos_mod_array,
                                            const FlavArrayFuncSet_t* DF_Set, 
                                            const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), l(l_in), m(m_in), pol_factor(&pol_factor_in), y_func_val( y_func ),
           cos_modulation(&pol_factor_in), DF_val(&pol_factor_in), FF_val(&pol_factor_in), // dummy values for now
           func_ptr( &SIDIS_2had_f1D1_t::func_00 )
      {

         {
            std::stringstream ss;
            ss << "f1 D_1 |" << l << ", " << m << ">";
            message = ss.str();
         };

         // invalid l,m combinations taken care of by nonzero = FALSE;

         if( DF_Set->Includes("f1") ){

            DF_val = DF_Set->Get_Val_Ptr("f1");

            if( l == 0 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_00;
                  cos_modulation = &cos_mod_array[Enum::P_00];

                  if( FF_Set->Includes("Re_D1_00") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_00");
                     nonzero = 1;
                  };
               };
            } else if ( l == 1 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_10;
                  cos_modulation = &cos_mod_array[Enum::P_10];
                  if( FF_Set->Includes("Re_D1_10") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_10");
                     nonzero = 1;
                  };
               } else if( m == 1 || m == -1 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_11;
                  cos_modulation = &cos_mod_array[Enum::P_11];
                  if( FF_Set->Includes("Re_D1_11") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_11");
                     nonzero = 1;
                  };
               };
            } else if ( l == 2 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_20;
                  cos_modulation = &cos_mod_array[Enum::P_20];
                  if( FF_Set->Includes("Re_D1_20") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_20");
                     nonzero = 1;
                  };
               } else if( m == 1 || m == -1 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_21;
                  cos_modulation = &cos_mod_array[Enum::P_21];
                  if( FF_Set->Includes("Re_D1_21") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_21");
                     nonzero = 1;
                  };
               } else if( m == 2 || m == -2 ){
                  func_ptr = &SIDIS_2had_f1D1_t::func_22;
                  cos_modulation = &cos_mod_array[Enum::P_22];
                  if( FF_Set->Includes("Re_D1_22") ){
                     FF_val = FF_Set->Get_Val_Ptr("Re_D1_22");
                     nonzero = 1;
                  };
               };
            };
         };

      };

      double SIDIS_2had_f1D1_t::Eval( const Var_t& var ){
         //std::cerr << "SIDIS_2had_f1D1_t::Eval(...): |" << l << ", " << m << "> " << (*cos_modulation) << std::endl;

//           if( (*cos_modulation) ){
//              std::cout << "c " << (this->*func_ptr)( var ) << ' ' << (*pol_factor) << ' ' << (*cos_modulation);
//              std::cout << ' ' << (*y_func_val) << ' ' << (*DF_val) << ' ' << (*FF_val) << std::endl;
//           };

//         std::cout << l << '.' << m << ' ' << (*cos_modulation) << ' ' << (*y_func_val) << ' ' << (*pol_factor) << ' ';
//         std::cout << (*DF_val) << ' ' << (*FF_val) << ' ' << (this->*func_ptr)( var ) << std::endl; 

//          std::cout << "f1D1\t" << cos(m*(var.phi_R-var.phi_kT)) << ' ' << (*pol_factor) << ' ' << (*cos_modulation);
//          std::cout << ' ' << (*y_func_val) << ' ' << (*DF_val) << ' ' << (*FF_val) << std::endl;

//          if( !var.integrating ){
//             cout << "+ " << flavor_string[ var.flavor ] << ' ';
//             cout << cos( m*( var.phi_R - var.kT) ) << ' ' ;
//             cout << (*pol_factor) * (*cos_modulation) * (*y_func_val) << ' ';
//             cout << (*DF_val) << ' ' << (*FF_val) << ' ';
//             cout << cos(m*(var.phi_R-var.phi_kT))
//                * (*pol_factor) * (*cos_modulation) * (*y_func_val) * (*DF_val) * (*FF_val) << endl;
//          };

//          if( !var.integrating )
//             cout << "a) " << (*DF_val) << endl;

         //return (this->*func_ptr)( var ) 
         return cos(m*(var.phi_R-var.phi_kT))
            * (*pol_factor) * (*cos_modulation) * (*y_func_val)
            * (*DF_val)
            * (*FF_val);

         //return cos(m*(var.phi_R-var.phi_kT)) * (*FF_val) * (*cos_modulation) ;
      };

   };
};






