
// the A,B,C and D functions of y

#include "Common.d/Var.h"

#include <cmath>

#ifdef MAYBE_INLINE
#undef MAYBE_INLINE
#endif

#ifdef INLINE_Y_FUNCTIONS
#define MAYBE_INLINE inline
#else
#define MAYBE_INLINE  
#endif

namespace TMDGen {
   namespace yFunctions {

      MAYBE_INLINE double A( const Var_t& var){
         double output = 1. - var.y + 0.5*var.y*var.y + 0.25*var.y*var.y*var.gamma_sq;
         output /= (1. + var.gamma_sq); // cannot be negative
         return output;
      };

      MAYBE_INLINE double B( const Var_t& var){
         double output = 1. - var.y - 0.25*var.y*var.y*var.gamma_sq;
         output /= (1. + var.gamma_sq); // cannot be negative
         return output;
      };

      MAYBE_INLINE double C( const Var_t& var){
         double output = (1. - 0.5*var.y);
         output *= var.y;
         output /= (1. + var.gamma_sq); // cannot be negative
         return output;
      };

      MAYBE_INLINE double V( const Var_t& var){
         double output = 1 - var.y - 0.25*var.y*var.y*var.gamma_sq;

         output = ( output > 0 ? sqrt(output) : 0 );

         if( output ){
            output *= 2.;
            output *= (2. - var.y);
            output /= (1. + var.gamma_sq); // cannot be negative
         };

         return output;
      };

      MAYBE_INLINE double W( const Var_t& var){
         double output = 1 - var.y - 0.25*var.y*var.y*var.gamma_sq;

         output = ( output > 0 ? sqrt(output) : 0 );

         if( output ){
            output *= 2.;
            output *= var.y;
            output /= (1. + var.gamma_sq); // cannot be negative
         };

         return output;
      };

   };
};
