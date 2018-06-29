 /*
   Class to construct, own, and maintain all needed distribution functions
*/


#ifndef _FULL_DF_SET_H_
#define _FULL_DF_SET_H_

#include "Common.d/ParseInput.h"
#include "Common.d/FlavArrayFuncSet.h"
#include <map>

namespace TMDGen {
   namespace DF {

      class Full_DF_Set_t : public FlavArrayFuncSet_t {
         void Make_CTEQ( const sgUtil::ParseInputReturn_t& parsed_input, std::string key, std::string option );
         void Make_Const( const sgUtil::ParseInputReturn_t& parsed_input, std::string DF_name, double val );
         void Make_LHAPDF( const sgUtil::ParseInputReturn_t& parsed_input, std::string DF_name, std::string line );

         void SetWrapper( const sgUtil::ParseInputReturn_t &parsed_input, std::string, int& iset, std::string& path );

      public:
         Full_DF_Set_t( const sgUtil::ParseInputReturn_t& parsed_input );
         virtual ~Full_DF_Set_t();
      };

   };
};

#endif
