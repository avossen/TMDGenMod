/*
  Simple exponential function, with constant mean and sigma
*/


#include "pT_distr.d/NonConst_Exp.h"
#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"

#include <cmath>
#include <iostream>

namespace TMDGen {
   namespace pT_distr {

      // by convention, function is written as f(x) = R/sqrt(pi) * exp( -R^2 x^2 )
      // with R^2 = pi/(4 mean^2)
      // and mean = C * x^alpha * ( 1 - x )^beta

      NonConst_Exp_t::NonConst_Exp_t( int model, double C_, bool is_kT_) : is_kT(is_kT_), C(C_) {

         if( C <= 0 ){
            throw Error::Constructing ("NonConst_Exp_t", "C <= 0" );
         };

         if( is_kT ){ // i.e. doing kT
            if( model == 0 ){ 
               // Boglione & Mulders [ PRD 60 (1999) 054007 ] 
               alpha = 0.27;
               beta = 0.20;
            } else if ( model == 1 ){
               // Hashi Set I (private communication)
               alpha = 0.5;
               beta = 0.5;
            } else if ( model == 2 ){
               // Hashi Set II (private communication)
               alpha = 0.25;
               beta = 0.15;
            } else if ( model == 3 ){
               // Hashi Set III (private communication)
               alpha = 0.2681605;
               beta = 0.1830445;
            } else {
               throw Error::Constructing ("NonConst_Exp_t", "Invalid Model Selection: must be 0-3 for kT distributions" );
            };
         } else {
            if( model == 0 ){ 
               // Boglione & Mulders [ PRD 60 (1999) 054007 ] 
               alpha = 0.68;
               beta = 0.48;
            } else {
               throw Error::Constructing ("NonConst_Exp_t", "Invalid Model Selection: must be 0 for pT distributions" );
            };
         };

         // no constraints
         min = 0;
         max = 0;
      };

      double NonConst_Exp_t::Eval(const TMDGen::Var_t& var) const{

         double x = is_kT ? var.z : var.x;
         double mean = C * pow( x, alpha ) * pow( 1.-x, beta );

         double minus_R_sq = - PI / 4 / mean / mean;
         double normfactor = 0.5 / mean;

         double temp = is_kT ? var.pT : var.kT; // pT or kT depending on is_kT

         temp *= temp;
         temp *= minus_R_sq;
         temp = exp( temp );

         temp *= normfactor;

         // Note: extra factor of z on outside due to kT vs. KT nuance
         if( is_kT )
            temp /= var.z;

         //std::cout << temp << ' ' << var.z << ' ' << var.kT << ' ' << mean << ' ' << minus_R_sq << '\t' << C << ' ' << alpha << ' ' << beta << std::endl;

         return temp;
      };
      
   };
};

