/*
   SIDIS pseudo scalar twist 2 cos phi term
*/

#ifndef _XSEC_TERM_SIDIS_1HAD_COSPHI_H_
#define _XSEC_TERM_SIDIS_1HAD_COSPHI_H_

#include "XSec_Term.d/XSec_Term.h"

#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"


namespace TMDGen {
   namespace XSec_Term {

      class SIDIS_1had_cosphi_t : public XSec_Term_t {

         // pointers to the precomputed values
         const double *pol_factor,       // polarization factor, such as S_T or beam polarization, etc.
            *y_func_val,                 // one of A(y), B(y), ...
            *h_1perp_val,                // DFs
            *htilde_val,
            *f1_val,
            *ftilde_perp_val,
            *D1_val,                     // FFs
            *Dtilde_perp_val,
            *H_1perp_val,
            *Htilde_val;

         static const double ZERO;
 
      public:
         SIDIS_1had_cosphi_t( const double& pol_factor_in, const double* y_func_array,
                              const FlavArrayFuncSet_t* DF_Set, 
                              const FlavArrayFuncSet_t* FF_Set );
         virtual ~SIDIS_1had_cosphi_t(){ /* */ };

         virtual double Eval( const Var_t& var );
      };
   };
};

#endif
