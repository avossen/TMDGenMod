
#include "Thrower.d/BasicThrower.h"
#include "Thrower.d/BasicThrower_2had.h"

#include "Common.d/Consts.h"
#include "Thrower.d/GSL_Integration.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_vegas.h>

#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

namespace TMDGen {
   namespace Thrower {


      BasicThrower_2had_t::BasicThrower_2had_t( RNG_t *r_in, int N_integration_warmup_, int N_integration_calls,
                                                int N_max_integration_calls_,
                                                const Var_t& min_in, const Var_t& max_in ) :
         BasicThrower_t( r_in,  N_integration_warmup_, N_integration_calls, N_max_integration_calls_, min_in, max_in )
      {
         width.x       = (max.x       - min.x);
         width.y       = (max.y       - min.y);
         width.z       = (max.z       - min.z);
         width.P_hperp = (max.P_hperp - min.P_hperp );
         width.pT      = (max.pT      - min.pT);
         width.had_0.M = (max.had_0.M - min.had_0.M);
         width.phi_h = TWO_PI;
         width.phi_pT = TWO_PI;
         width.phi_R = TWO_PI;
         width.psi = TWO_PI;
         width.cos_vartheta = 2.;

         Vinv = 1./width.x/width.y/width.z/width.P_hperp/width.pT/width.had_0.M/width.phi_h/width.phi_pT/width.phi_R/width.psi/width.cos_vartheta;
//          cerr << "Volume is " << 1./Vinv << ' ' << width.x << ' ' << width.y << ' ' << width.z << ' ' << width.P_hperp << std::endl;
//          cerr << ' ' << width.pT << ' ' << width.had_0.M << ' ' << width.phi_h << ' ' << width.phi_pT << std::endl;
//          cerr << ' ' << width.phi_R << ' ' << width.psi << ' ' << width.vartheta << endl;

      };


      // specialization for different dependent variables
      void BasicThrower_2had_t::Init_GSL_Func( gsl_monte_function *gsl_F, double* &min_array, double* &max_array ){
         // set correct function to integrate

         // will integrate ang. int term -- no phi_S !!!

         gsl_F->f = &GSL_IntegrationFunction_2h;
         gsl_F->dim = 11;

         min_array = new double [ gsl_F->dim ];
         max_array = new double [ gsl_F->dim ];

         min_array[0] = min.x;
         min_array[1] = min.y;
         min_array[2] = min.z;
         min_array[3] = min.P_hperp;
         min_array[4] = 0.;            // phi_h
         min_array[5] = min.pT;
         min_array[6] = 0.;            // phi_{p_T}
         min_array[7] = min.had_0.M;
         min_array[8] = 0.;            // phi_{R}
         min_array[9] = -1.;           // cos_vartheta
         min_array[10] = 0;            // psi   

         max_array[0] = max.x;
         max_array[1] = max.y;
         max_array[2] = max.z;
         max_array[3] = max.P_hperp;
         max_array[4] = TWO_PI;        // phi_h
         max_array[5] = max.pT; 
         max_array[6] = TWO_PI;        // phi_{p_T}
         max_array[7] = max.had_0.M;
         max_array[8] = TWO_PI;        // phi_{R}
         max_array[9] = 1.;            // cos_vartheta
         max_array[10] = TWO_PI;       // psi   

      };

      // for throwing according to different distributions
      void BasicThrower_2had_t::ThrowVariables( Var_t& var, double& pdf_val ){
         var.x       = r->EvalUnif()*width.x       + min.x;
         var.y       = r->EvalUnif()*width.y       + min.y;
         var.z       = r->EvalUnif()*width.z       + min.z;
         var.P_hperp = r->EvalUnif()*width.P_hperp + min.P_hperp;
         var.phi_h   = r->EvalUnif()*TWO_PI;
         var.pT      = r->EvalUnif()*width.pT      + min.pT;
         var.phi_pT  = r->EvalUnif()*TWO_PI;
         var.had_0.M = r->EvalUnif()*width.had_0.M + min.had_0.M;
         var.phi_R   = r->EvalUnif()*TWO_PI;
         var.cos_vartheta = r->EvalUnif()*2. - 1.; // r->EvalUnif()*PI;  // ( 2.0 * ( r->EvalUnif() ) - 1. );
         var.psi     = r->EvalUnif()*TWO_PI;

//          std::cout << "a ";
//          std::cout << var.vartheta << std::endl;

         pdf_val = Vinv;
         //std::cout << "t " << pdf_val << endl;

      };

   };
};
