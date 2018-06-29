/*
    Naive model for dihadron terms
*/

#ifndef _XSEC_TERM_SIDIS_2HAD_NAIVE_H_
#define _XSEC_TERM_SIDIS_2HAD_NAIVE_H_

#include "XSec_Term.d/XSec_Term.h"

#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {
   namespace XSec_Term {

      class SIDIS_2had_Naive_t : public XSec_Term_t {

         bool is_cos;
         int l, m, n_h, n_S;

         // pointers to the precomputed values
         const double &pol_factor,       // polarization factor, such as S_T or beam polarization, etc.
            *y_func_val,       // one of A(y), B(y), ...
            *y_func_val_A,     // A(x,y)
            *cos_modulation;   // functions of cos_vartheta

         int mom_num;

         std::string pol_state;

         struct params_t {
            double A;
            double x_slope;
            double z_alpha;
            double z_beta;
            double pT_slope;
            double kT_slope_alpha;
            double kT_slope_beta;
            double kT_slope_gamma;
            double M_h_a1;
            double M_h_a2;

            params_t();
         } params;

      public:
         SIDIS_2had_Naive_t( std::string& pol_state, int mom_num_in, std::string& params,
                             const double& pol_factor_in, const double* y_func_array, const double* cos_mod_array_in );
         virtual ~SIDIS_2had_Naive_t(){ /* */ };

         virtual double Eval( const Var_t& var );
      };
   };
};

#endif
