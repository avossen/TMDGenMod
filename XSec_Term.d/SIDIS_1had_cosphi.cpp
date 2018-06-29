/*
   SIDIS pseudo scalar twist 3 cos phi term
*/

#include "XSec_Term.d/SIDIS_1had_cosphi.h"

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

      const double SIDIS_1had_cosphi_t::ZERO = 0;

      // constructor
      SIDIS_1had_cosphi_t::SIDIS_1had_cosphi_t( const double& pol_factor_in, const double* y_func,
                                                const FlavArrayFuncSet_t* DF_Set, const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           h_1perp_val(&ZERO), htilde_val(&ZERO), f1_val(&ZERO), ftilde_perp_val(&ZERO),
           D1_val(&ZERO), Dtilde_perp_val(&ZERO), H_1perp_val(&ZERO), Htilde_val(&ZERO)
      {
         bool includes_DF_h = 0;
         bool includes_DF_f = 0;

         message = "'cos(phi)' Tw. 3";

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


         nonzero = 0;
         if( FF_Set->Includes( "D1" )){
            D1_val = FF_Set->Get_Val_Ptr( "D1" );
            nonzero |= includes_DF_f;
         };
         if( FF_Set->Includes( "D_perp~" )){
            Dtilde_perp_val = FF_Set->Get_Val_Ptr( "D_perp~" );
            nonzero |= includes_DF_f;
         };
         if( FF_Set->Includes( "H_1^perp" )){
            H_1perp_val = FF_Set->Get_Val_Ptr( "H_1^perp" );
            nonzero |= includes_DF_h;
         };
         if( FF_Set->Includes( "H~" )){
            Htilde_val = FF_Set->Get_Val_Ptr( "H~" );
            nonzero |= includes_DF_h;
         };

      };

      double SIDIS_1had_cosphi_t::Eval( const Var_t& var ){
         double factor1 = var.pT * var.pT / PROTON_MASS_SQUARED * (*h_1perp_val);
         if( (*htilde_val) )
            factor1 += var.x * (*htilde_val);
         factor1 *= (*H_1perp_val);

         if( (*Dtilde_perp_val) )
            factor1 += var.had_0.M / PROTON_MASS * (*f1_val) * (*Dtilde_perp_val);

         if( factor1 )
            factor1 *= var.kT / var.had_0.M;

         double factor2 = (*f1_val);
         if( (*ftilde_perp_val) )
            factor2 += var.x * (*ftilde_perp_val);

         factor2 *= (*D1_val);

         if( (*Htilde_val) )
            factor2 += PROTON_MASS / var.had_0.M * (*h_1perp_val) * (*Htilde_val) / var.z;

         factor2 *= var.pT / PROTON_MASS;

         return ( factor1 * cos( var.phi_kT ) +
                  factor2 * cos( var.phi_pT ) ) *
            (-2.)*PROTON_MASS/sqrt(var.Q2) * 
            (*pol_factor) * (*y_func_val);
      };


   };
};
