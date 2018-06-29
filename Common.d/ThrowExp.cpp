/*
  Class to throw x according to normalized power law distribution
*/

#include "Common.d/ThrowExp.h"
#include "Common.d/Exceptions.h"
#include <cmath>


namespace TMDGen {

   ThrowExp_t::ThrowExp_t( double alpha_, double min_, double max_ ) : minus_alpha(-alpha_), min(min_), max(max_) {
      if( min < 0 )
         throw Error::Constructing("ThrowExp_t", "Negative minimum range");
      if( max < 0 )
         throw Error::Constructing("ThrowExp_t", "Negative maximum range");
      if( max <= min )
         throw Error::Constructing("ThrowExp_t", "Inverted range");
      if( minus_alpha >= 0 )
         throw Error::Constructing("ThrowExp_t", "non-positive alpha");

      N = 1./(1 - exp( minus_alpha * ( max - min ) ));
   };

   // returns inverse cdf at point r ( r in [0,1) ), i.e. cdf at return value = r
   double ThrowExp_t::InverseCDF( double r ){
      double x = 0;
      if( r > 0 && r < 1 ){
         x = min + log( 1 - r/N )/minus_alpha;
      };
      return x;
   };

   double ThrowExp_t::PDF( double x ){
      return ( x > min && x < max ) ? -N*minus_alpha*exp( minus_alpha*( x - min ) ) : 0;
   };
};
