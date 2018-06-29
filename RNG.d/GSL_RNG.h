/*
   GSL random number generator
*/

#ifndef _GSL_RNG_H_
#define _GSL_RNG_H_

#include "RNG.d/RNG.h"

#include <gsl/gsl_rng.h>

namespace TMDGen {

   class GSL_RNG_t : public RNG_t {
      gsl_rng* r;

   public:
      static int last_seed;

      GSL_RNG_t( const gsl_rng_type *T =  gsl_rng_ranlxs2, int seed = -1 );
      virtual ~GSL_RNG_t();

      virtual double EvalUnif();
   };

};

#endif
