/*
  Struct to hold parameter sets for Spec_Ia
*/

#include "Spec_Ia_Params.h"

#include <string>
#include <cmath>
#include <iostream>
using std::cerr;
using std::endl;

namespace TMDGen {
   namespace FF {

      namespace DiHad_Spec_Ia {


         const int Spec_Ia_Params_t::N_Params = 15;

         Spec_Ia_Params_t::Spec_Ia_Params_t( const double* X ){
            alpha_H1perp = 1;

            alpha_s = X[0];
            beta_s = X[1];
            gamma_s = X[2];
            alpha_p = X[3];
            beta_p = X[4];
            gamma_p = X[5];
            f_s = X[6];
            f_rho = X[7];
            f_omega = X[8];
            f_omega_prime = X[9];
            M_s_factor = X[10];
            alpha_H1perp = X[11];
            kT_alpha = X[12];
            kT_beta = X[13];
            kT_gamma = X[14];

//             cerr << "f_s           = " << f_s           << endl;
//             cerr << "f_rho         = " << f_rho         << endl;
//             cerr << "f_omega       = " << f_omega       << endl;
//             cerr << "f_omega_prime = " << f_omega_prime << endl;
//            cerr << "alpha_H1perp = " << alpha_H1perp << endl;
         };


         Spec_Ia_Params_t::Spec_Ia_Params_t( const std::string& set ) {
            // not yet programed, so just clear things
            (*this) = Spec_Ia_Params_t();
         };

         Spec_Ia_Params_t::Spec_Ia_Params_t() {
            // clear things
            // not the fastest way, but minimal chance for introducing bugs

            double X[N_Params];
            for( int i=0; i<N_Params; ++i )
               X[i] = 0;

            (*this) = Spec_Ia_Params_t( X );
         };

         Spec_Ia_Params_t::Spec_Ia_Params_t(ParamSets_t set ) {

            // zero everything via the other constructor
            double X[N_Params];
            for( int i=0; i<N_Params; ++i )
               X[i] = 0;

            (*this) = Spec_Ia_Params_t( X );

            // set according to switches
            switch(set){
            case NONE:
               break;

            case BR06: 
               // from arXiv:hep-ph/0608037v1
               alpha_s = 2.60;
               beta_s = -0.751;
               gamma_s = -0.193;
               alpha_p = 7.07;
               beta_p = -0.038;
               gamma_p = -0.085;
               f_s = 1197;
               f_rho = 93.5;
               f_omega = 0.63;
               f_omega_prime = 75.2;
               M_s_factor = 2.97;

               // from arXiv:0812.0611
               alpha_H1perp = 0.32;
               break;

            case G10: 
               // best guess
               alpha_s = 0.8;
               beta_s = -0.751;
               gamma_s = -0.193;
               alpha_p = 2.0;
               beta_p = -0.038;
               gamma_p = -0.085;
               f_s = 100;
               f_rho = 100;
               f_omega = 50;
               f_omega_prime = 50;
               M_s_factor = 1; //2.97;

               // from arXiv:0812.0611
               alpha_H1perp = 0.32;
               break;

            };
         };
      };
   };
};
