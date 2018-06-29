/*
  Function to initial LHAPDF from input file
*/

#ifndef NO_LHAPDF
#ifndef _LHAPDF_INIT_H_
#define _LHAPDF_INIT_H_

#include "Common.d/ParseInput.h"

namespace TMDGen {
   namespace DF {

      int LHAPDF_init( const sgUtil::ParseInputReturn_t& parsed_input );

   };
};

#endif
#endif
