/*
  Function to initial LHAPDF from input file
*/

#ifndef NO_LHAPDF

#include "Common.d/ParseInput.h"
#include "DF.d/LHAPDF_init.h"
#include "LHAPDF/LHAPDF.h"


#include <sstream>
#include <string>
#include <iostream>

using std::cerr;
using std::endl;

#ifndef LHAPDF_MAX_SETS
#define LHAPDF_MAX_SETS 3
#endif

namespace TMDGen {
   namespace DF {

      int LHAPDF_init( const sgUtil::ParseInputReturn_t& parsed_input ){
         int ierr = 0;

         int verb = 0;
         sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find("LHAPDF_Verbosity");
         if( it != parsed_input.end() ){
            verb = atoi( it->second.data() );
         };
         LHAPDF::setVerbosity( static_cast< LHAPDF::Verbosity >( verb ) );

         it = parsed_input.find("LHAPDF_Extrapolate");
         if( it != parsed_input.end() ){
            LHAPDF::extrapolate(atoi( it->second.data() ));
         };

         it = parsed_input.find("LHAPDF_Path");
         if( it != parsed_input.end() ){
            LHAPDF::setPDFPath( it->second );
         };


         for( int i=0; i<LHAPDF_MAX_SETS && !ierr; ++i ){
            std::string base;
            {
               std::stringstream ss;
               ss << "LHAPDF_" << i << '_';
               base = ss.str();
            };

            it = parsed_input.find( base + "Name" );
            if( it != parsed_input.end() ){
               std::stringstream ss( it->second );

               std::string name;
               std::string grid;
               int member = 0;
               LHAPDF::SetType type = LHAPDF::LHPDF;

               ss >> name >> grid >> member;
               if( grid == "LHPDF" ){
                  type = LHAPDF::LHPDF;
               } else if ( grid == "LHGRID" ){
                  type = LHAPDF::LHPDF;
               } else {
                  cerr << "\t\tLHAPDF_init: on nset " << i << "\n\t\t\t";
                  cerr << "Unknown grid/pdf option: '" << grid << "', must be 'LHPDF' or 'LHGRID'." << endl;
                  ierr = 1;
               };

               if( !ierr ){
                  LHAPDF::initPDFSet( i, name, type, member );
               };
            };
         };

         return ierr;
      };

   };
};

#endif
