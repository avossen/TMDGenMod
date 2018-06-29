
/*
  Struct to keep track of parameters
*/

#ifndef _h1_TORINO_PARAMS_H_
#define _h1_TORINO_PARAMS_H_

#include "Common.d/Flavor.h"
#include <string>

namespace TMDGen {
   namespace DF {

      struct h1_Torino_Params_t {
         double N[ GMC_TRANS_N_FLAVORS ];
         double alpha, beta;
         double ave_pT2;

         h1_Torino_Params_t( std::string init_code="Ringberg" );
      };

   };
};

#endif
