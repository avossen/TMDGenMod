/*
    Partial wave expansion of SIDIS dihadron cosphi term
*/

#include "XSec_Term.d/SIDIS_2had_cosphi.h"

#include "Common.d/Consts.h"
#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"
#include "Common.d/Yfunctions.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"

#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
using std::cout;
using std::cerr;
using std::endl;

namespace TMDGen {
   namespace XSec_Term {

      const double SIDIS_2had_cosphi_t::ZERO = 0;

      // the weights and other prefactors in the pT kT integral

      double SIDIS_2had_cosphi_t::func_00( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( var.phi_h - var.phi_kT );
         factor2 *= cos( var.phi_h - var.phi_pT );

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_h );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_11( const Var_t& var, double factor1, double factor2 ){
         factor1 *= 0.5;
         factor2 *= cos( var.phi_pT - var.phi_kT );

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_R );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_10( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( var.phi_h - var.phi_kT );
         factor2 *= cos( var.phi_h - var.phi_pT );

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_h );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_1m1( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( 2.*var.phi_h - 2.*var.phi_kT );
         factor2 *= cos( 2.*var.phi_h - var.phi_pT - var.phi_kT );

         factor1 += factor2;
         factor1 *= 2.*cos( 2.*var.phi_h - var.phi_R );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_22( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( var.phi_h - var.phi_kT );                   // signs for these two don't match dissertation v0.6
         factor2 *= cos( var.phi_h + var.phi_pT - 2.*var.phi_kT );   // bug in the dissertation

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_h - 2.*var.phi_R );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_21( const Var_t& var, double factor1, double factor2 ){
         factor1 *= 0.5;
         factor2 *= cos( var.phi_pT - var.phi_kT );

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_R );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_20( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( var.phi_h - var.phi_kT );
         factor2 *= cos( var.phi_h - var.phi_pT );

         factor1 += factor2;
         factor1 *= 2.*cos( var.phi_h );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_2m1( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( 2.*var.phi_h - 2.*var.phi_kT );
         factor2 *= cos( 2.*var.phi_h - var.phi_pT - var.phi_kT );

         factor1 += factor2;
         factor1 *= 2.*cos( 2.*var.phi_h - var.phi_R );

         return factor1;
      };

      double SIDIS_2had_cosphi_t::func_2m2( const Var_t& var, double factor1, double factor2 ){
         factor1 *= cos( 3.*var.phi_h - 3.*var.phi_kT );
         factor2 *= cos( 3.*var.phi_h - var.phi_pT - 2.*var.phi_kT );

         factor1 += factor2;
         factor1 *= 2.*cos( 3.*var.phi_h - 2.*var.phi_R );

         return factor1;
      };


      // constructor
      SIDIS_2had_cosphi_t::SIDIS_2had_cosphi_t( int l, int m_in,
                                            const double& pol_factor_in, const double* y_func, const double* cos_mod_array,
                                            const FlavArrayFuncSet_t* DF_Set, 
                                            const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           cos_modulation(&ZERO), h_1perp_val(&ZERO), htilde_val(&ZERO), f1_val(&ZERO), ftilde_perp_val(&ZERO),
           D1_val(&ZERO), Dtilde_perp_val(&ZERO), H_1perp_val(&ZERO), Htilde_val(&ZERO),
           func_ptr( &SIDIS_2had_cosphi_t::func_00 ), m(m_in)
      {
         // invalid l,m combinations taken care of by nonzero = FALSE;

         bool includes_DF_h = 0;
         bool includes_DF_f = 0;

         {
            std::stringstream ss;
            ss << "'cos(phi)' |" << l << ", " << m << ">";
            message = ss.str();
         };

         if( DF_Set->Includes("h_1^perp") ){
            h_1perp_val = DF_Set->Get_Val_Ptr("h_1^perp");
            includes_DF_h = 1;
         };
         if( DF_Set->Includes("h~") ){
            htilde_val = DF_Set->Get_Val_Ptr("h~");
            includes_DF_h = 1;
         };
         if( DF_Set->Includes("f1") ){
            f1_val = DF_Set->Get_Val_Ptr("f1");
            includes_DF_f = 1;
         };
         if( DF_Set->Includes("f_perp~") ){
            ftilde_perp_val = DF_Set->Get_Val_Ptr("f_perp~");
            includes_DF_f = 1;
         };

         if( includes_DF_h || includes_DF_f ){
            std::string stemp, stemp2;

            if( l == 0 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_00;
                  cos_modulation = &cos_mod_array[Enum::P_00];
                  stemp = "_00";
               };
            } else if ( l == 1 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_10;
                  cos_modulation = &cos_mod_array[Enum::P_10];
                  stemp = "_10";
                  stemp2 = "_10";
               } else if( m == 1 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_11;
                  cos_modulation = &cos_mod_array[Enum::P_11];
                  stemp = "_11";
                  stemp2 = "_11";
               } else if( m == -1 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_1m1;
                  cos_modulation = &cos_mod_array[Enum::P_11];
                  stemp = "_11";
                  stemp2 = "_1m1";
               };
            } else if ( l == 2 ){
               if( m == 0 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_20;
                  cos_modulation = &cos_mod_array[Enum::P_20];
                  stemp = "_20";
                  stemp2 = "_20";
               } else if( m == 1 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_21;
                  cos_modulation = &cos_mod_array[Enum::P_21];
                  stemp = "_21";
                  stemp2 = "_21";
               } else if( m == -1 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_2m1;
                  cos_modulation = &cos_mod_array[Enum::P_21];
                  stemp = "_21";
                  stemp2 = "_2m1";
               } else if( m == 2 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_22;
                  cos_modulation = &cos_mod_array[Enum::P_22];
                  stemp = "_22";
                  stemp2 = "_22";
               } else if( m == -2 ){
                  func_ptr = &SIDIS_2had_cosphi_t::func_2m2;
                  cos_modulation = &cos_mod_array[Enum::P_22];
                  stemp = "_22";
                  stemp2 = "_2m2";
               };
            };

            nonzero = 0;
            if( FF_Set->Includes( std::string("Re_D1") + stemp ) ){
               D1_val = FF_Set->Get_Val_Ptr( std::string("Re_D1") + stemp );
               nonzero |= includes_DF_f;
            };
            if( FF_Set->Includes( std::string("D_perp~") + stemp2 ) ){
               Dtilde_perp_val = FF_Set->Get_Val_Ptr( std::string("D_perp~") + stemp2 );
               nonzero |= includes_DF_f;
            };
            if( FF_Set->Includes( std::string("H_1^perp") + stemp2 ) ){
               H_1perp_val = FF_Set->Get_Val_Ptr( std::string("H_1^perp") + stemp2 );
               nonzero |= includes_DF_h;
            };
            if( FF_Set->Includes( std::string("H~") + stemp2 ) ){
               Htilde_val = FF_Set->Get_Val_Ptr( std::string("H~") + stemp2 );
               nonzero |= includes_DF_h;
            };

         };
      };

      double SIDIS_2had_cosphi_t::Eval( const Var_t& var ){
         double factor1 = var.pT * var.pT / PROTON_MASS_SQUARED * (*h_1perp_val);
         if( (*htilde_val) )
            factor1 += var.x * (*htilde_val);
         factor1 *= (*H_1perp_val);

         if( (*Dtilde_perp_val) )
            factor1 += var.had_0.M / PROTON_MASS * (*f1_val) * (*Dtilde_perp_val) / var.z;

         if( factor1 )
            factor1 *= var.kT / var.had_0.M;

         double factor2 = (*f1_val);
         if( (*ftilde_perp_val) )
            factor2 += var.x * (*ftilde_perp_val);

         factor2 *= (*D1_val);

         if( (*Htilde_val) )
            factor2 += PROTON_MASS / var.had_0.M * (*h_1perp_val) * (*Htilde_val) / var.z;

         factor2 *= var.pT / PROTON_MASS;

         if( factor1 != factor1 )
            throw Error::SanityCheckFailure( "XSec_Term::SIDIS_2had_cosphi::Eval(...)",
                                             "factor 1 is NaN" );

         if( factor2 != factor2 )
            throw Error::SanityCheckFailure( "XSec_Term::SIDIS_2had_cosphi::Eval(...)",
                                             "factor 2 is NaN" );


//          std::cout << "cosphi: " <<
// //             (this->*func_ptr)( var, factor1, factor2 ) << ' ' << 
// //             (-2.)*PROTON_MASS/sqrt(var.Q2)  << ' ' << 
// //             (*pol_factor)  << ' ' <<  (*cos_modulation)  << ' ' << 
//             (*y_func_val) << std::endl;

         //return (this->*func_ptr)( var, factor1, factor2 ) * 

         double val = 
            ( factor1 * cos( var.phi_kT + m*(var.phi_R - var.phi_kT) ) + 
                  factor2 * cos( var.phi_pT + m*(var.phi_R - var.phi_kT) ) ) *
            (-2.)*PROTON_MASS/sqrt(var.Q2) * 
            (*pol_factor) * (*cos_modulation) * (*y_func_val);

         if( val != val )
            throw Error::SanityCheckFailure( "XSec_Term::SIDIS_2had_cosphi::Eval(...)",
                                             "returning NaN" );

            return val;
      };

   };
};
