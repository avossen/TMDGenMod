/*
   GSL random number generator
*/

#include "RNG.d/GSL_RNG.h"

#include <gsl/gsl_rng.h>
#include <time.h>


namespace TMDGen {

   GSL_RNG_t::GSL_RNG_t( const gsl_rng_type *T, int seed ) {
      r = gsl_rng_alloc(T);

      if( seed < 1 )
         seed = time(0);

      if( seed == last_seed )
         ++seed;

      last_seed = seed;

      gsl_rng_set( r, seed );
   };

   GSL_RNG_t::~GSL_RNG_t(){
      gsl_rng_free( r );
   };

   double GSL_RNG_t::EvalUnif(){
      return gsl_rng_uniform( r );
   };

   int GSL_RNG_t::last_seed = 1;

};
