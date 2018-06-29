/*

Defines enum to encode flavor types

*/


#include "Common.d/Flavor.h"
#include <string>

namespace TMDGen {

   // string values of the flavors
   const std::string flavor_string[ GMC_TRANS_N_FLAVORS ] = { "u", "d", "s", "c", "b", "ubar", "dbar", "sbar", "cbar", "bbar", "gluon" };

   // convert from enum to Lund codes
   const int Flav_2_Lund[ GMC_TRANS_N_FLAVORS ] = { 2, 1, 3, 4, 5, -2, -1, -3, -4, -5, 0 };

};

