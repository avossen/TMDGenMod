/*
    Partial wave expansion of SIDIS dihadron cos2phi term
*/

#ifndef _XSEC_TERM_SIDIS_2HAD_COS2PHI_H_
#define _XSEC_TERM_SIDIS_2HAD_COS2PHI_H_

#include "XSec_Term.d/XSec_Term.h"

#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"


namespace TMDGen {
   namespace XSec_Term {

      class SIDIS_2had_cos2phi_t : public XSec_Term_t {

         // pointers to the precomputed values
         const double *pol_factor,       // polarization factor, such as S_T or beam polarization, etc.
            *y_func_val,       // one of A(y), B(y), ...
            *cos_modulation,   // functions of cos_vartheta
            *DF_val,           // distribution function value
            *FF_val;           // fragmentation function value

         // pointers to functions for terms varying with l, m
         double (SIDIS_2had_cos2phi_t::*func_ptr)( const Var_t& var );
         double func_00( const Var_t& var );
         double func_11( const Var_t& var );
         double func_10( const Var_t& var );
         double func_1m1( const Var_t& var );
         double func_22( const Var_t& var );
         double func_21( const Var_t& var );
         double func_20( const Var_t& var );
         double func_2m1( const Var_t& var );
         double func_2m2( const Var_t& var );

         int m;
 
      public:
         SIDIS_2had_cos2phi_t( int l, int m,
                                 const double& pol_factor_in, const double* y_func_array, const double* cos_mod_array_in,
                                 const FlavArrayFuncSet_t* DF_Set, 
                                 const FlavArrayFuncSet_t* FF_Set );
         virtual ~SIDIS_2had_cos2phi_t(){ /* */ };

         virtual double Eval( const Var_t& var );
      };
   };
};

#endif
