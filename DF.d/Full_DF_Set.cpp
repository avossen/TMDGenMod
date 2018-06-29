/*
  Class to construct, own, and maintain all needed distribution functions
*/


#include "DF.d/Full_DF_Set.h"

#ifndef NO_LHAPDF
#include "DF.d/LHAPDF_init.h"
#include "DF.d/LHAPDF_DF.h"
#endif

#include "DF.d/CTEQ_DF.h"
#include "DF_Wrapper.d/CTEQ.h"

#include "DF.d/f1_GRV98_.h"
#include "DF_Wrapper.d/GRV.h"

#include "DF.d/g1_GRSV2000_.h"
#include "DF_Wrapper.d/GRSV.h"

// Torino Sivers, Boer Mulders, Transversity
#include "DF.d/Torino_nT_odd.h"
#include "DF.d/h1_Torino.h"

#include "DF.d/f1_BCR08_.h"
#include "DF.d/f1Tperp_BCR08_.h"
#include "DF.d/g1L_BCR08_.h"
#include "DF.d/g1T_BCR08_.h"
#include "DF.d/h1_BCR08_.h"
#include "DF.d/h1Lperp_BCR08_.h"
#include "DF.d/h1perp_BCR08_.h"
#include "DF.d/h1Tperp_BCR08_.h"

#include "DF.d/Const_DF.h"

#include "Common.d/ParseInput.h"
#include "Common.d/Exceptions.h"
#include "Common.d/FlavArrayFunc.h"

#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>


namespace TMDGen {
   namespace DF {

      Full_DF_Set_t::Full_DF_Set_t( const sgUtil::ParseInputReturn_t& parsed_input ){
         // look for all fragmentation directives and make models

         sgUtil::ParseInputReturn_const_iterator_t it;

         SetWrapper( parsed_input, "CTEQ", DF_Wrapper::CTEQ_t::iset, DF_Wrapper::CTEQ_t::path );
         SetWrapper( parsed_input, "GRV98", DF_Wrapper::GRV_t::iset,  DF_Wrapper::GRV_t::path );
         SetWrapper( parsed_input, "GRSV2000",  DF_Wrapper::GRSV_t::iset,  DF_Wrapper::GRSV_t::path );

#ifndef NO_LHAPDF
         if( LHAPDF_init( parsed_input ) )
            throw Error::Constructing( "Full_DF_Set_t", "error initializing LHAPDF_init" );
#endif

         // place everything in try statement, to deallocate as needed if exception thrown
         try {

            it = parsed_input.find("f1");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "f1", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "f1", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "f1", it->second.substr(6).data() );
#endif

               } else if( it->second == "BCR08" ){
                  if( func_map["f1"] )
                     delete func_map["f1"];

                  func_map["f1"] = new DF::f1_BCR08_t();

               } else if( it->second == "GRV98" ){

                  // also need a line for the pT distribution part
                  sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find( "f1_pT");

                  if( it == parsed_input.end() ){
                     throw Error::Constructing( "Full_DF_Set_t", "'f1' directive of 'GRV98' requires additional 'f1_pT' directive." );
                  };

                  if( func_map["f1"] )
                     delete func_map["f1"];

                  func_map["f1"] = new DF::f1_GRV98_t( it->second );

               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'f1' directive: '") + it->second + "'" );
               };
            };

            // Sivers
            it = parsed_input.find("f_1T^perp");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "f_1T^perp", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "f_1T^perp", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "f_1T^perp", it->second.substr(6).data() );
#endif

               } else if ( it->second.substr( 0, 6 ) == "Torino" ){
                  if( func_map["f_1T^perp"] )
                     delete func_map["f_1T^perp"];

                  func_map["f_1T^perp"] = new DF::Torino_nT_odd_t( "f_1T^perp", it->second.substr(7) );

               } else if( it->second == "BCR08"){
                  if( func_map["f_1T^perp"] )
                     delete func_map["f_1T^perp"];

                  func_map["f_1T^perp"] = new DF::f1Tperp_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'f_1T^perp' directive: '") + it->second + "'" );
               };
            };


            // unclear: g_1L, g_1T, or g_1 and g_1T^perp
            /*
            // Add in GRSV for g_1

            // LO in A_LL
            it = parsed_input.find("g_1L");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "g_1L", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "g_1L", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "g_1L", it->second.substr(6).data() );
#endif

               } else if( it->second == "BCR08"){
                  if( func_map["g_1L"] )
                     delete func_map["g_1L"];

                  func_map["g_1L"] = new DF::g1L_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'g_1L' directive: '") + it->second + "'" );
               };
            };

            // LO in A_LT (single had)
            // also A_UT with Im(D_1)
            it = parsed_input.find("g_1T");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "g_1T", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "g_1T", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "g_1T", it->second.substr(6).data() );
#endif

               } else if( it->second == "BCR08"){
                  if( func_map["g_1T"] )
                     delete func_map["g_1T"];

                  func_map["g_1T"] = new DF::g1T_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'g_1T' directive: '") + it->second + "'" );
               };
            };
            */

            // Transversity
            it = parsed_input.find("h_1");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "h_1", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "h_1", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "h_1", it->second.substr(6).data() );
#endif

               } else if ( it->second.substr( 0, 6 ) == "Torino" ){
                  if( func_map["h_1"] )
                     delete func_map["h_1"];

                  func_map["h_1"] = new DF::h1_Torino_t( it->second.substr(7) );
               } else if( it->second == "BCR08"){
                  if( func_map["h_1"] )
                     delete func_map["h_1"];

                  func_map["h_1"] = new DF::h1_BCR08_t();

               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'h_1' directive: '") + it->second + "'" );
               };
            };

            // LO in A_LU
            it = parsed_input.find("h_1L^perp");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "h_1L^perp", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "h_1L^perp", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "h_1L^perp", it->second.substr(6).data() );
#endif

               } else if( it->second == "BCR08"){
                  if( func_map["h_1L^perp"] )
                     delete func_map["h_1L^perp"];

                  func_map["h_1L^perp"] = new DF::h1Lperp_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'h_1L^perp' directive: '") + it->second + "'" );
               };
            };

            // Pretzelocity
            it = parsed_input.find("h_1T^perp");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "h_1T^perp", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "h_1T^perp", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "h_1T^perp", it->second.substr(6).data() );
#endif

               } else if( it->second == "BCR08"){
                  if( func_map["h_1T^perp"] )
                     delete func_map["h_1T^perp"];

                  func_map["h_1T^perp"] = new DF::h1Tperp_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'h_1T^perp' directive: '") + it->second + "'" );
               };
            };

            // Boer-Mulders
            it = parsed_input.find("h_1^perp");
            if( it != parsed_input.end() ){
               if( it->second.substr( 0, 4 ) == "CTEQ" ){
                  Make_CTEQ( parsed_input, "h_1^perp", it->second.substr(4) );

               } else if ( it->second.substr( 0, 5 ) == "Const" ){
                  Make_Const( parsed_input, "h_1^perp", atof( it->second.substr(5).data() ) );

#ifndef NO_LHAPDF
               } else if ( it->second.substr( 0, 6 ) == "LHAPDF" ){
                  Make_LHAPDF( parsed_input, "h_1^perp", it->second.substr(6).data() );
#endif

               } else if ( it->second.substr( 0, 6 ) == "Torino" ){
                  if( func_map["h_1^perp"] )
                     delete func_map["h_1^perp"];

                  func_map["h_1^perp"] = new DF::Torino_nT_odd_t( "h_1^perp", it->second.substr(7) );
               } else if( it->second == "BCR08"){
                  if( func_map["h_1^perp"] )
                     delete func_map["h_1^perp"];

                  func_map["h_1^perp"] = new DF::h1perp_BCR08_t();
               } else {
                  throw Error::Constructing( "Full_DF_Set_t", std::string("invalid 'h_1^perp' directive: '") + it->second + "'" );
               };
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

      Full_DF_Set_t::~Full_DF_Set_t(){
         /* nothing to do */
      };

      void Full_DF_Set_t::Make_CTEQ( const sgUtil::ParseInputReturn_t& parsed_input, std::string key, std::string option ){
         // clear white space and get first word
         std::stringstream ss( option );
         ss >> option;

         // also need a line for the pT distribution part
         // as CTEQ does not include pT
         sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find( key + "_pT");
         if( it == parsed_input.end() ){
            throw Error::Constructing( "Full_DF_Set_t",
                                       std::string("'") + key + "' directive of 'CTEQ' requires additional '" + key + "_pT' directive." );
         };

         // note: if error constructing, will throw exception
         // so no need to save old f1 in case this one fails.
         // should probably never happen that an f1 already exists
         if( func_map[key] )
            delete func_map[key];

         func_map[key] = new DF::CTEQ_DF_t ( key, option, it->second );
      };

      void Full_DF_Set_t::Make_Const( const sgUtil::ParseInputReturn_t& parsed_input, std::string DF_name, double val ){

         // also need a line for the pT distribution part
         sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find( DF_name + "_pT");

         if( it == parsed_input.end() ){
            throw Error::Constructing( "Full_DF_Set_t", std::string("'") + DF_name + "' directive of 'Const' requires additional '"
                                       + DF_name + "_pT' directive." );
         };

         if( func_map[DF_name] )
            delete func_map[DF_name];

         func_map[DF_name] = new DF::Const_DF_t( val, it->second.data() );
      };

      void Full_DF_Set_t::Make_LHAPDF( const sgUtil::ParseInputReturn_t& parsed_input, std::string DF_name, std::string line ){
#ifndef NO_LHAPDF
         std::string weight_option;
         int nset;
         int member;

         std::stringstream ss( line );
         ss >> weight_option >> nset >> member;

         // also need a line for the pT distribution part
         sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find( DF_name + "_pT");

         if( it == parsed_input.end() ){
            throw Error::Constructing( "Full_DF_Set_t", std::string("'") + DF_name + "' directive of 'LHAPDF' requires additional '"
                                       + DF_name + "_pT' directive." );
         };

         if( func_map[DF_name] )
            delete func_map[DF_name];

         func_map[DF_name] = new DF::LHAPDF_DF_t( DF_name, weight_option, nset, member, it->second );
#endif
      };

      void Full_DF_Set_t::SetWrapper( const sgUtil::ParseInputReturn_t &parsed_input, std::string key, int& iset, std::string& path ){
         sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find( key + "_iSet");
         if( it != parsed_input.end() )
            iset = atoi( it->second.data() );

         it = parsed_input.find( key + "_Path");
         if( it != parsed_input.end() )
            path = it->second.data();
      };

   };
};
