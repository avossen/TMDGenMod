/*
    Twist 2 SIDIS pseudoscalar f1D1 term
*/

#include "XSec_Term.d/SIDIS_1had_f1D1_.h"

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

      // constructor
      SIDIS_1had_f1D1_t::SIDIS_1had_f1D1_t(  const double& pol_factor_in, const double* y_func,
                                            const FlavArrayFuncSet_t* DF_Set, const FlavArrayFuncSet_t* FF_Set )
         : XSec_Term_t(0), pol_factor(&pol_factor_in), y_func_val( y_func ),
           DF_val(&pol_factor_in), FF_val(&pol_factor_in) // dummy values for now
      {
         message = "f1 D_1";

         if( DF_Set->Includes("f1") ){
            DF_val = DF_Set->Get_Val_Ptr("f1");
            if( FF_Set->Includes("D1") ){
               FF_val = FF_Set->Get_Val_Ptr("D1");
               nonzero = 1;
            };
         };
      };

      double SIDIS_1had_f1D1_t::Eval( const Var_t& var ){
         return (*pol_factor) * (*y_func_val) * (*DF_val) * (*FF_val);
      };

   };
};
