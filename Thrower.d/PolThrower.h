/*
  Base class for throwing poliables according to some (easy) distribution
*/

#ifndef _POLTHROWER_H_
#define _POLTHROWER_H_

#include "Common.d/Var.h"

namespace TMDGen {
   namespace Thrower {

      class PolThrower_t {

      protected:
         PolThrower_t(){ /* */ };

      public:
         virtual ~PolThrower_t(){ /* */ };
         virtual void Throw( Var_t& pol ) = 0;
      };

   };
};

#endif
