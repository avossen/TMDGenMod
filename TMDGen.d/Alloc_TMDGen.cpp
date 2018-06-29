
// common allocator for children of TMDGen_t

#include "TMDGen.d/Alloc_TMDGen.h"
#include "TMDGen.d/TMDGen.h"
#include "TMDGen.d/InputCodes.h"
#include "TMDGen.d/SIDIS_1had.h"
#include "TMDGen.d/SIDIS_2had.h"

#include "Common.d/ParseInput.h"

#include <iostream>
using std::cerr;
using std::endl;

namespace TMDGen {

   TMDGen_t* Alloc_TMDGen( const char* parameter_file_name ){

      sgUtil::ParseInputReturn_t parsed_input;

      if( sgUtil::ParseInput( parameter_file_name, parsed_input ) ){
         cerr << "Alloc_TMDGen(...): Error with input file '" << parameter_file_name << "'" << endl;
         return 0;
      };

      return Alloc_TMDGen( parsed_input );
   };

   TMDGen_t* Alloc_TMDGen( const sgUtil::ParseInputReturn_t& parsed_input ){

      TMDGen_t* gmc_trans = 0;

      // check process
      sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find("Process");

      if( it == parsed_input.end() ){
         cerr << "Alloc_TMDGen(...): no 'Process' directive given." << endl;
         return 0;
      };

      if( it->second == "SIDIS" ){
         it = parsed_input.find("Final_State");
         if( it == parsed_input.end() ){
            cerr << "Alloc_TMDGen(...): no 'Final_State' directive given." << endl;
            return 0;
         };

         if( it->second == InputCode::SingleHadron ){
            gmc_trans = new SIDIS_1had_t( parsed_input );
         } else if ( it->second == InputCode::Dihadron ){
            gmc_trans = new SIDIS_2had_t( parsed_input );
         } else {
            cerr << "Alloc_TMDGen(...) ERROR: only '" << InputCode::SingleHadron << "' and '";
            cerr << InputCode::Dihadron << "' final states for SIDIS currently supported." << endl;
            // throw Error::NotYetProgrammed("Alloc_TMDGen(...): only 'Single Hadron' final states for SIDIS currently supported.");
            return 0;
         };
      } else {
         cerr << "Alloc_TMDGen(...): only SIDIS processes currently supported" << endl;
         // throw Error::NotYetProgrammed("Alloc_TMDGen(...): only SIDIS processes currently supported");
      };

      return gmc_trans;
   };
};

