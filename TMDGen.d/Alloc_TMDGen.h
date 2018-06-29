#ifndef _ALLOC_GMC_TRANS_H_
#define _ALLOC_GMC_TRANS_H_

#include "Common.d/ParseInput.h"
#include "TMDGen.d/TMDGen.h"

namespace TMDGen {
   TMDGen_t* Alloc_TMDGen( const char* parameter_file_name );
   TMDGen_t* Alloc_TMDGen( const sgUtil::ParseInputReturn_t& parsed_input );
};

#endif
