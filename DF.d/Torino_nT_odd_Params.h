
/*
  Struct to keep track of parameters
*/

#ifndef _TORINO_NT_ODD_PARAMS_H_
#define _TORINO_NT_ODD_PARAMS_H_

#include "Common.d/Flavor.h"
#include <string>

namespace TMDGen {
   namespace DF {

      struct Torino_nT_odd_Params_t {
         double A[ GMC_TRANS_N_FLAVORS ];
         double alpha[ GMC_TRANS_N_FLAVORS ];
         double beta, M1_sq, ave_pT_sq;
         double lambda[ GMC_TRANS_N_FLAVORS ];

         Torino_nT_odd_Params_t( std::string init_code="" );
      };

   };
};

#endif
