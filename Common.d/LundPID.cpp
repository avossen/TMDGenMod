/*

Defines enum to encode lund particle types

*/

#include "Common.d/LundPID.h"

namespace TMDGen {

   int Convert_Lund_to_Geant( LundPID_t pid ){
      GeantPID_t output = INVALID_G;
      switch( pid ){
      case PI_PLUS:
         output = PI_PLUS_G;
         break;
      case PI_MINUS:
         output = PI_MINUS_G;
         break;
      case PI_ZERO:
         output = PI_ZERO_G;
         break;
      case K_PLUS:
         output = K_PLUS_G;
         break;
      case K_MINUS:
         output = K_MINUS_G;
         break;
      case PHOTON:
         output = PHOTON_G;
         break;
      default:
         output = INVALID_G;
         break;
      };

      return static_cast< int >( output );
   };

};
