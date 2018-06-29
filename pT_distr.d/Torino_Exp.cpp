/*
  Simple exponential function, with constant mean and sigma
  Matches Torino groups function--slighly different notation and normalization
*/


#include "pT_distr.d/Torino_Exp.h"
#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"

#include <cmath>

namespace TMDGen {
   namespace pT_distr {

      // by convention, function is written as f(x) = 1/(a sqrt(pi)) * exp( -x^2/a )
      // with a = < pT^2 > or < kT^2 >

      Torino_Exp_t::Torino_Exp_t( double ave_pT_sq_, bool is_kT_) : is_kT(is_kT_), ave_pT_sq(ave_pT_sq_) {

         if( !ave_pT_sq ){
            throw Error::Constructing ("Torino_Exp_t", "ave_pT_sq = 0" );
         };

         normfactor = 1. / ave_pT_sq / PI;

         // no constraints
         min = 0;
         max = 0;
      };

      double Torino_Exp_t::Eval(const TMDGen::Var_t& var) const{
         // fixed bug in following statement (was backwards) June 4, 2010
         double temp = is_kT ? var.kT : var.pT; // pT or kT depending on is_kT 

         temp *= temp;
         temp /= ave_pT_sq;
         temp = exp( -temp );

         temp *= normfactor;

         return temp;
      };
      
   };
};

