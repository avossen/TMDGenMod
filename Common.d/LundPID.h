/*

Defines enum to encode lund particle types

*/

#ifndef _LUNDPID_H_
#define _LUNDPID_H_

namespace TMDGen {

   enum LundPID_t { PI_PLUS = 211, PI_MINUS = -211, PI_ZERO = 111, K_PLUS = 321, K_MINUS = -321, PHOTON = 22, INVALID_LUND_PID = 0 };
   enum GeantPID_t { PI_PLUS_G = 8, PI_MINUS_G = 9, PI_ZERO_G = 7, K_PLUS_G = 11, K_MINUS_G = 12, ELEC_G = 3, POSI_G = 2, PHOTON_G = 1, INVALID_G = 0 };

   int Convert_Lund_to_Geant( LundPID_t pid );

};

#endif
