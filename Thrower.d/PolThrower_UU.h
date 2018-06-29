/*
  Fully unpolarized (trivial case)
*/

#ifndef _POLTHROWER_UU_H_
#define _POLTHROWER_UU_H_

#include "Common.d/Var.h"

namespace TMDGen {
   namespace Thrower {

      class PolThrower_UU_t : public PolThrower_t {
      public:
         PolThrower_UU_t(){ /* */ };
         virtual ~PolThrower_UU_t(){ /* */ };
         virtual void Throw( Var_t& pol ){ /* */ };
      };

   };
};

#endif
