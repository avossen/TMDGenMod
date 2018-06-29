// pi0 decay

#ifndef _DECAY_PI0_
#define _DECAY_PI0_

#include "Common.d/FourMom.h"
#include "RNG.d/GSL_RNG.h"


namespace TMDGen {
   namespace Decay {

      class pi0 {

      public:
         pi0();
         virtual ~pi0();

         virtual int operator() ( const FourMom_t& p_pi0, FourMom_t& gamma_1, FourMom_t& gamma_2, TMDGen::GSL_RNG_t *r ) = 0;
      };
   };
};

#endif
