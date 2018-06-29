/*
   Defines static array of quark charges, squared, intexed by flavor_t
*/

#include "Common.d/Flavor.h"
#include "Common.d/QuarkChargeSquared.h"

namespace TMDGen {

   const double QuarkChargeSquared[ GMC_TRANS_N_FLAVORS ] = { 4./9., 1./9., 1./9., 4./9., 1./9., 4./9., 1./9., 1./9., 4./9., 1./9., 0 };

};
