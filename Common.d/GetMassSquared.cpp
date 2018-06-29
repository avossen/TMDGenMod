/*
  Return mass squared for given Lund PID code
*/

#include "Common.d/LundPID.h"
#include "Common.d/Consts.h"

namespace TMDGen {

   double GetMassSquared( const LundPID_t& pid ){
      double mass_sq = 0;

      switch( pid ){
      case PI_PLUS:
      case PI_MINUS:
      case PI_ZERO:
         mass_sq = PION_MASS_SQUARED;
         break;
      case K_PLUS:
      case K_MINUS:
         mass_sq = KAON_MASS_SQUARED;
         break;
      default:
         mass_sq = 0;
         break;
      };

      return mass_sq;
   };

};
