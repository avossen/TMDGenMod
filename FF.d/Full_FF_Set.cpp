/*
   Class to construct, own, and maintain all needed distribution functions
*/

#include "FF.d/Full_FF_Set.h"
#include "FF.d/D1_fDSS.h"
#include "FF.d/D1_Kretzer.h"
#include "FF.d/D1_Spec_Ia.h"
//#include "FF.d/D1_Spec_Ib.h"
#include "FF.d/H1perp_Spec_Ia.h"
#include "FF.d/FF_Spec_Ib.h"
#include "FF.d/Prop_2had.h"
#include "FF.d/Const_FF.h"


#include "TMDGen.d/InputCodes.h"

#include "Common.d/ParseInput.h"
#include "Common.d/Exceptions.h"
#include "Common.d/SetPID.h"

#include <string>
#include <sstream>
#include <iostream>
#include <cstdlib>

using std::cerr;
using std::cout;
using std::endl;

namespace TMDGen {
   namespace FF {

      Full_FF_Set_t::Full_FF_Set_t( const sgUtil::ParseInputReturn_t& parsed_input ){
         // look for all fragmentation directives and make models

         sgUtil::ParseInputReturn_const_iterator_t it;

         it = parsed_input.find("Final_State");
         if( it == parsed_input.end() )
            throw Error::Constructing( "Full_FF_Set_t", "no 'Final_State' directive given" );


         try {
            if( it->second == InputCode::SingleHadron ){

               /**** SINGLE HADRON CASE ***/


               LundPID_t had_PID;
               it = parsed_input.find("Hadron_PID");
               if( SetPID( it->second, had_PID ) ){
                  throw Error::Constructing( "Full_FF_Set_t", std::string("invalid Hadron_PID: '") + it->second + "'" );
               };

               it = parsed_input.find("D1");
               if( it != parsed_input.end() ){
                  if( it->second.substr( 0, 4 ) == "fDSS" ){

                     std::stringstream ss( it->second.substr( 4 ) );
                     int order = 0;
                     std::string path = "";

                     ss >> order >> path;

                     // also need a line for the kT distribution part
                     // as fDSS does not include kT
                     it = parsed_input.find("D1_kT");
                     if( it == parsed_input.end() ){
                        throw Error::Constructing( "Full_FF_Set_t", "'D1' directive of 'fDSS' requires additional 'D1_kT' directive." );
                     };

                     if( func_map["D1"] )
                        delete func_map["D1"];

                     func_map["D1"] = new FF::D1_fDSS_t( had_PID, order, path.data(), it->second.data() );

                  } else if( it->second.substr( 0, 7 ) == "Kretzer" ){

                     std::stringstream ss( it->second.substr( 7 ) );
                     int order = 0;
                     std::string path = "";

                     ss >> order >> path;

                     // also need a line for the kT distribution part
                     // as Kretzer does not include kT
                     it = parsed_input.find("D1_kT");
                     if( it == parsed_input.end() ){
                        throw Error::Constructing( "Full_FF_Set_t", "'D1' directive of 'Kretzer' requires additional 'D1_kT' directive." );
                     };

                     if( func_map["D1"] )
                        delete func_map["D1"];

                     func_map["D1"] = new FF::D1_Kretzer_t( had_PID, order, path.data(), it->second.data() );

                  } else if ( it->second.substr( 0, 5 ) == "Const" ){

                     double val = atof( it->second.substr(5).data() );

                     // also need a line for the kT distribution part
                     // as Kretzer does not include kT
                     it = parsed_input.find("D1_kT");
                     if( it == parsed_input.end() ){
                        throw Error::Constructing( "Full_FF_Set_t", "'D1' directive of 'Const' requires additional 'D1_kT' directive." );
                     };

                     if( func_map["D1"] )
                        delete func_map["D1"];

                     func_map["D1"] = new FF::Const_FF_t( "D1", val, it->second.data() );
                     std::cerr << "\t\tD1: Const" << std::endl;

                  } else {
                     throw Error::Constructing( "Full_FF_Set_t", std::string("invalid 'D1' directive: '") + it->second + "'" );
                  };
               } else {
                  throw Error::Constructing( "Full_FF_Set_t", std::string("No 'D1' directive given") );
               };
            } else if (  it->second == InputCode::Dihadron ){

               /**** DIHADRON CASE ***/


               // get PIDs

               LundPID_t had_1_PID, had_2_PID;

               it = parsed_input.find("Hadron_1_PID");
               if( it == parsed_input.end() )
                  throw Error::Constructing( "XSec::SIDIS_2had", "no 'Hadron_1_PID' directive given" );

               if( SetPID( it->second, had_1_PID ) )
                  throw Error::Constructing( "XSec::SIDIS_2had", std::string("invalid Hadron_1_PID: '") + it->second + "'" );

               it = parsed_input.find("Hadron_2_PID");
               if( it == parsed_input.end() )
                  throw Error::Constructing( "XSec::SIDIS_2had", "no 'Hadron_2_PID' directive given" );

               if( SetPID( it->second, had_2_PID ) )
                  throw Error::Constructing( "XSec::SIDIS_2had", std::string("invalid Hadron_2_PID: '") + it->second + "'" );

	       // check for parameters
	       double Spec_Ia_Params[DiHad_Spec_Ia::Spec_Ia_Params_t::N_Params];
	       int use_Spec_Ia_Params = 0;
	       it = parsed_input.find("FF_Spec_Ia_Params");
	       if( it != parsed_input.end() ){
                  use_Spec_Ia_Params = 1;

                  int n = DiHad_Spec_Ia::Spec_Ia_Params_t::N_Params;

                  std::stringstream ss;
                  ss << it->second;
                  for( int i=0; i<n; ++i )
                     ss >> Spec_Ia_Params[i];
	       };

	       double Spec_Ib_Params[3][DiHad_Spec_Ia::Spec_Ia_Params_t::N_Params];
	       int use_Spec_Ib_Params[3] = {0,0,0};
               for( int i=0; i<3; ++i ){
                  std::stringstream ss1;
                  ss1 << "FF_Spec_Ib_Params_" << i;

                  it = parsed_input.find( ss1.str() );

                  if( it != parsed_input.end() ){
                     use_Spec_Ib_Params[i] = 1;

                     int n = DiHad_Spec_Ia::Spec_Ia_Params_t::N_Params;

                     std::stringstream ss2;
                     ss2 << it->second;
                     for( int j=0; j<n; ++j )
                        ss2 >> Spec_Ib_Params[i][j];
                  };
	       };

               
               std::string Re_D1_name[6] = { "Re_D1_00", "Re_D1_11", "Re_D1_10",
                                             "Re_D1_22", "Re_D1_21", "Re_D1_20" };



               it = parsed_input.find( Re_D1_name[0] );
               if( it == parsed_input.end() ){
                  throw Error::Constructing( "Full_FF_Set_t",
                                             std::string("No 'Re_D1_00' directive given") );
               };



               for( int i=0; i<6; ++i ){
                  std::string& name = Re_D1_name[i];

                  //cerr << "Looking for '" << Re_D1_name[i] << "'" << endl;

                  it = parsed_input.find( name );
                  if( it != parsed_input.end() ){
                     int l = atoi( name.substr(6,1).data() );
                     int m = atoi( name.substr(7,1).data() );

                     if( it->second.substr( 0, 6 ) == "Spec_I" ){
                        if( func_map[name] )
                           delete func_map[name];

                        if( it->second[6] == 'a' ){
                           FF::D1_Spec_Ia_t *D1 = new FF::D1_Spec_Ia_t( had_1_PID, had_2_PID, l, m );

                           if( use_Spec_Ia_Params )
                              D1->CopyParams(  DiHad_Spec_Ia::Spec_Ia_Params_t( Spec_Ia_Params )  );

                           func_map[name] = D1;
//                         } else if ( it->second[6] == 'B' ){
//                            FF::D1_Spec_Ib_t *D1 = new FF::D1_Spec_Ib_t( had_1_PID, had_2_PID, l, m );

//                            for( int i=0; i<3; ++i )
//                               if( use_Spec_Ib_Params[i] )
//                                  D1->CopyParams(  i, DiHad_Spec_Ia::Spec_Ia_Params_t( Spec_Ib_Params[i] )  );

//                            func_map[name] = D1;
                        } else {
                           FF::FF_Spec_Ib_t *D1 = new FF::FF_Spec_Ib_t( had_1_PID, had_2_PID, l, m, 0 );

                           for( int i=0; i<3; ++i )
                              if( use_Spec_Ib_Params[i] )
                                 D1->CopyParams(  i, DiHad_Spec_Ia::Spec_Ia_Params_t( Spec_Ib_Params[i] )  );

                           func_map[name] = D1;
                        };
                     } else if ( it->second.substr( 0, 4 ) == "Prop" ){
                        std::stringstream ss( it->second );
                        std::string prop_to_key;
                        double params[9];
                        int i;

                        ss >> prop_to_key >> prop_to_key;  // clear first word and save key

                        for( i=0; i<9 && ss; ++i ){
                           ss >> params[i];
                           ss.peek();
                        };

                        //cerr << "*** i = " << i << endl;
                        func_map[name] = new FF::Prop_2had_t( name, l, m, prop_to_key, this, i, params );

                     } else if ( it->second.substr( 0, 5 ) == "Const" ){

                        double val = atof( it->second.substr(5).data() );

                        // also need a line for the kT distribution part
                        it = parsed_input.find( name + "_kT" );
                        if( it == parsed_input.end() ){
                           throw Error::Constructing(
                                                     "Full_FF_Set_t",
                                                     std::string("'") + name + 
                                                     "' directive of 'Const' requires additional" +
                                                     "'" + name + "_kT' directive." );
                        };

                        if( func_map[name] )
                           delete func_map[name];

                        func_map[name] = new FF::Const_FF_t( name, val, it->second.data() );
                     } else {
                        throw Error::Constructing( "Full_FF_Set_t",
                                                   std::string("invalid 'D1' directive: '") +
                                                   it->second + "'" );
                     };
                  };
               };

               // Not yet programmed Im_D1

               std::string Collins_name[9] =
                  { "H_1^perp_00", 
                    "H_1^perp_11", 
                    "H_1^perp_10", 
                    "H_1^perp_1m1", 
                    "H_1^perp_22", 
                    "H_1^perp_21", 
                    "H_1^perp_20", 
                    "H_1^perp_2m1", 
                    "H_1^perp_2m2" };


               for( int i=0; i<9; ++i ){
                  std::string& name = Collins_name[i];

                  it = parsed_input.find( name );
                  if( it != parsed_input.end() ){
                     int l = atoi( name.substr(9,1).data() );
                     int m = ( (name.at(10) == 'm') ? -atoi( name.substr(11,1).data() ) : atoi( name.substr(10,1).data() ) );

                     if( it->second.substr( 0, 6 ) == "Spec_I" ){

                        if( it->second[6] == 'a' ){
                           if( func_map[name] )
                              delete func_map[name];

                           FF::H1perp_Spec_Ia_t *H1perp = new FF::H1perp_Spec_Ia_t( had_1_PID, had_2_PID, l, m );
                           if( use_Spec_Ia_Params )
                              H1perp->CopyParams(  DiHad_Spec_Ia::Spec_Ia_Params_t( Spec_Ia_Params )  );

                           func_map[name] = H1perp;
                        } else {
                           FF::FF_Spec_Ib_t *H1perp = new FF::FF_Spec_Ib_t( had_1_PID, had_2_PID, l, m, 1 );

                           for( int i=0; i<3; ++i )
                              if( use_Spec_Ib_Params[i] )
                                 H1perp->CopyParams(  i, DiHad_Spec_Ia::Spec_Ia_Params_t( Spec_Ib_Params[i] )  );

                           func_map[name] = H1perp;
                        };
                     } else if ( it->second.substr( 0, 4 ) == "Prop" ){

                        std::stringstream ss( it->second );
                        std::string prop_to_key;
                        double params[9];
                        int i;

                        ss >> prop_to_key >> prop_to_key;  // clear first word and save key

                        for( i=0; i<9 && ss; ++i ){
                           ss >> params[i];
                           ss.peek();
                        };

                        //cerr << "*** i = " << i << endl;
                        func_map[name] = new FF::Prop_2had_t( name, l, m, prop_to_key, this, i, params );

                     } else if ( it->second.substr( 0, 5 ) == "Const" ){

                        double val = atof( it->second.substr(5).data() );

                        // also need a line for the kT distribution part
                        it = parsed_input.find( name + "_kT" );
                        if( it == parsed_input.end() ){
                           throw Error::Constructing(
                                                     "Full_FF_Set_t",
                                                     std::string("'") + name + 
                                                     "' directive of 'Const' requires additional" +
                                                     "'" + name + "_kT' directive." );
                        };

                        if( func_map[name] )
                           delete func_map[name];

                        func_map[name] = new FF::Const_FF_t( name, val, it->second.data() );
                     } else {
                        throw Error::Constructing( "Full_FF_Set_t", std::string("invalid '") + name
                                                   + "' directive: '" + it->second + "'" );
                     };
                  };
               };
            } else {
               throw Error::Constructing( "TMDGen_t", std::string("invalid 'Final_State' directive: '") + it->second + "': must be '"
                                          + InputCode::SingleHadron + "' or '" + InputCode::Dihadron + "." );
            };
         }
         catch( std::exception& e ){
            Free();

            // rethrow
            throw;
         };

         Prepare_for_Precompute();
         DisplayConstrMsg();
      };

      Full_FF_Set_t::~Full_FF_Set_t(){
         /* nothing to do */
      };

   };
};
