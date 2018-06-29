/*
   Class to construct, own, and maintain all needed fragmentation functions
*/


#ifndef _FULL_FF_SET_H_
#define _FULL_FF_SET_H_

#include "Common.d/ParseInput.h"
#include "Common.d/FlavArrayFuncSet.h"
#include <map>

namespace TMDGen {
   namespace FF {

      class Full_FF_Set_t : public FlavArrayFuncSet_t {
      public:
         Full_FF_Set_t( const sgUtil::ParseInputReturn_t& parsed_input );
         ~Full_FF_Set_t();
      };

   };
};

#endif
