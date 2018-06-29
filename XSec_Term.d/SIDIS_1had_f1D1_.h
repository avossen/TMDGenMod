/*
    Twist 2 SIDIS pseudoscalar f1D1 term
*/

#ifndef _XSEC_TERM_SIDIS_1HAD_F1D1_H_
#define _XSEC_TERM_SIDIS_1HAD_F1D1_H_

#include "XSec_Term.d/XSec_Term.h"

#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"


namespace TMDGen {
   namespace XSec_Term {

      class SIDIS_1had_f1D1_t : public XSec_Term_t {

         // pointers to the precomputed values
         const double *pol_factor,       // polarization factor, such as S_T or beam polarization, etc., times e^2 and other things
            *y_func_val,       // one of A(y), B(y), ...
            *DF_val,           // distribution function value
            *FF_val;           // fragmentation function value

      public:
         SIDIS_1had_f1D1_t( const double& pol_factor_in, const double* y_func_array,
                            const FlavArrayFuncSet_t* DF_Set, const FlavArrayFuncSet_t* FF_Set );
         virtual ~SIDIS_1had_f1D1_t(){ /* */ };

         virtual double Eval( const Var_t& var );
      };
   };
};

#endif
