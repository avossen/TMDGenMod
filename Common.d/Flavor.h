/*

Defines enum to encode flavor types & string array of names

*/

#ifndef _FLAVOR_H_
#define _FLAVOR_H_

#include <string>

#define GMC_TRANS_N_FLAVORS 11

namespace TMDGen {

   enum flavor_t { UP_FLAV, DOWN_FLAV, STR_FLAV, CHM_FLAV, BOT_FLAV, ANTI_UP_FLAV, ANTI_DOWN_FLAV, ANTI_STR_FLAV, ANTI_CHM_FLAV, ANTI_BOT_FLAV, GLUON_FLAV };

   enum GRV_flavor_t { u_idx, d_idx, ubar_idx, dbar_idx, s_idx, glue_idx };

   // string values of the flavors
   extern const std::string flavor_string [ GMC_TRANS_N_FLAVORS ];

   // Lund codes of the flavors
   extern const int Flav_2_Lund[ GMC_TRANS_N_FLAVORS ];

};

#endif
