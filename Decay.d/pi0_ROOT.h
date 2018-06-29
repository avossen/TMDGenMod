// pi0 decay

#ifndef _DECAY_PI0_ROOT_
#define _DECAY_PI0_ROOT_

#include "Decay.d/pi0_.h"
#include "Common.d/FourMom.h"

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

namespace TMDGen {
   namespace Decay {

      class pi0_ROOT : public pi0 {
         TGenPhaseSpace genps;

      public:
         pi0_ROOT() : pi0() { /* */ };
         virtual ~pi0_ROOT() { /* */ };

         virtual int operator() ( const FourMom_t& p_pi0, FourMom_t& gamma_1, FourMom_t& gamma_2, TMDGen::GSL_RNG_t *r );
      };
   };
};

#endif
