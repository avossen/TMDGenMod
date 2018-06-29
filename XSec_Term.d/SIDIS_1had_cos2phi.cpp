/*
   SIDIS pseudo scalar twist 2 cos2phi moment
*/

#include "XSec_Term.d/SIDIS_1had_cos2phi.h"

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

      // constructor
      SIDIS_1had_cos2phi_t::SIDIS_1had_cos2phi_t( const double& pol_factor_in, const double* y_func,
                                                  const FlavArrayFuncSet_t* DF_Set, 
                                                  const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           DF_val(&pol_factor_in), FF_val(&pol_factor_in) // dummy values for now
      {
         message = "h_1^perp H_1^perp";

         if( DF_Set->Includes("h_1^perp") ){
            DF_val = DF_Set->Get_Val_Ptr("h_1^perp");
            if( FF_Set->Includes("H_1^perp") ){
               FF_val = FF_Set->Get_Val_Ptr("H_1^perp");
               nonzero = 1;
            };
         };
      };

      double SIDIS_1had_cos2phi_t::Eval( const Var_t& var ){
         return cos( var.phi_pT + var.phi_kT ) *
            (-var.pT)*var.kT/PROTON_MASS/var.had_0.M *
            (*pol_factor) * (*y_func_val) * (*DF_val) * (*FF_val);
      };

   };
};
