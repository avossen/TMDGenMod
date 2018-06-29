
#include "TMDGen.d/TMDGen.h"
#include "TMDGen.d/Alloc_TMDGen.h"

#include <iostream>
using std::cerr;
using std::endl;

int main( int argc, char* argv[] ){

   if( argc != 4 && argc != 3 ){
      cerr << "ERROR: Usage " << argv[0] << " <param_file_name> <N_events> [N_trials]" << endl;
      return 127;
   };

   const char* param_file_name = argv[1];
   int N = atoi( argv[2] );
   int N_trials = 100*N;
   if( argc > 3 )
      N_trials = atoi( argv[3] );

   cerr << endl << "Creating TMDGen object" << endl;
   TMDGen::TMDGen_t *generator = TMDGen::Alloc_TMDGen( param_file_name );

   if( !generator ){
      cerr << "ERROR instanciating TMDGen class" << endl;
      return 127;
   };

   cerr << endl << "Initializing" << endl;
   if( generator->Initialize() ){
      cerr << "ERROR initializing" << endl;
      return 127;
   };

   if( !generator->OutputCodes() ){
      cerr << "ERROR: no output codes set" << endl;
      return 127;
   };


   cerr << endl << "Generating Events" << endl;

   int n_out = 0;
   int i = 0;

   for( i=0; i<N_trials && n_out < N; ++i ){
      // generate var
      if( !generator->GenerateEvent() ){

         // compute other variables
         if( !generator->ReconstructEvent() ){

            // save the data
            generator->WriteOutEvent();
            ++n_out;
         };
      };
   };

   cerr << "After " << i << " trials, generated " << n_out << " events." << endl;

   delete generator;

   return 0;
};
