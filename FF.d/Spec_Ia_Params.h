/*
  Struct to hold parameter sets for Spec_Ia
*/

#ifndef _SPEC_IA_PARAMS_H_
#define _SPEC_IA_PARAMS_H_

#include <string>

namespace TMDGen {
   namespace FF {

      namespace DiHad_Spec_Ia {

         enum ParamSets_t { NONE, BR06, G10 };

         struct Spec_Ia_Params_t {
            double alpha_s, beta_s, gamma_s;
            double alpha_p, beta_p, gamma_p;
            double f_s, f_rho, f_omega, f_omega_prime;
            double M_s_factor;
            double alpha_H1perp;  // extra mult. factor on H_1^perp
            double kT_alpha, kT_beta, kT_gamma;

            const static int N_Params;

            Spec_Ia_Params_t();
            Spec_Ia_Params_t(ParamSets_t set );
            Spec_Ia_Params_t( const double* X );
            Spec_Ia_Params_t( const std::string& set );

         };

      };
   };
};

#endif
