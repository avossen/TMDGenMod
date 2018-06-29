/*
   Child of TMDGen_t, for SIDIS production of single hadrons
*/

#ifndef _SIDIS_1HAD_H_
#define _SIDIS_1HAD_H_

#include "TMDGen.d/TMDGen.h"
#include "Common.d/ParseInput.h"
#include "Common.d/Var.h"

#include <string>

namespace TMDGen {

   class SIDIS_1had_t : public TMDGen_t {

      // for output to ROOT
      virtual int Init_Root(const std::string& filename );

      // to ensure correct process and final state for given child
      virtual void Construct_Child( const sgUtil::ParseInputReturn_t& input );

      // for specific details of reconstructing event
      virtual int ReconstructEvent_details();

   public:
      SIDIS_1had_t( const char* parameter_file_name );
      SIDIS_1had_t( const sgUtil::ParseInputReturn_t& input );
      virtual ~SIDIS_1had_t();


   };
};

#endif
