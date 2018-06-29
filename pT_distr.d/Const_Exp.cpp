/*
  Simple exponential function, with constant mean and sigma
*/


#include "pT_distr.d/Const_Exp.h"
#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"

#include <cmath>

namespace TMDGen {
   namespace pT_distr {

      // by convention, function is written as f(x) = R/sqrt(pi) * exp( -R^2 x^2 )
      // with R^2 = pi/(4 alpha^2)

      Const_Exp_t::Const_Exp_t( double alpha, bool is_kT_) : is_kT(is_kT_) {

         if( !alpha ){
            throw Error::Constructing ("Const_Exp_t", "alpha = 0" );
         };

         minus_R_sq = - PI / 4 / alpha / alpha;
         normfactor = 0.5 / alpha;

         // no constraints
         min = 0;
         max = 0;
      };

      double Const_Exp_t::Eval(const TMDGen::Var_t& var) const{
         // fixed bug in following statement (was backwards) June 4, 2010
         double temp = is_kT ? var.kT : var.pT; // pT or kT depending on is_kT 

         if( is_kT )
            temp *= var.z;

         temp *= temp;
         temp *= minus_R_sq;
         temp = exp( temp );

         temp *= normfactor;

         if( is_kT )
            temp *= var.z;

         return temp;
      };
      
   };
};

