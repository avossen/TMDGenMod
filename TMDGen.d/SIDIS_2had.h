/*
   Child of TMDGen_t, for SIDIS production of dihadrons
*/

#ifndef _SIDIS_2HAD_H_
#define _SIDIS_2HAD_H_


#include "Common.d/Enums.h"
#include "Common.d/ParseInput.h"
#include "Common.d/Var.h"

#include "TMDGen.d/TMDGen.h"

#include <string>

namespace TMDGen {

   class SIDIS_2had_t : public TMDGen_t {

      // for output to ROOT
      virtual int Init_Root(const std::string& filename );
      virtual int Output_to_Root();
      double h1_P, h2_P;                                     // so can change signs if negatively charged particle
      int h1_iLType, h2_iLType, e2_iLType;
      // int output_FFs;

      // to ensure correct process and final state for given child
      virtual void Construct_Child( const sgUtil::ParseInputReturn_t& input );

      // for specific details of reconstructing event
      virtual int ReconstructEvent_details();

   public:
      SIDIS_2had_t( const char* parameter_file_name );
      SIDIS_2had_t( const sgUtil::ParseInputReturn_t& input );
      virtual ~SIDIS_2had_t();

   };
};

#endif
