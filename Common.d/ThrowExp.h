/*
  Class to throw x according to normalized power law distribution
*/

#ifndef _THROWEXP_H_
#define _THROWEXP_H_

namespace TMDGen {

   class ThrowExp_t {
      const double minus_alpha, min, max;
      double N;

   public:
      ThrowExp_t( double alpha_, double min_, double max_ );

      double InverseCDF( double r );
      double PDF( double x );

   };
};

#endif
