
#include <stdexcept>
#include <string>

#include "Exceptions.h"
#include <iostream>

namespace TMDGen {
   namespace Error {

      Base_t::Base_t( const std::string& prefix, const std::string& class_name, const std::string& msg ) :
         std::runtime_error( prefix + " '" + class_name + "': " + msg ) {
         std::cerr << what() << std::endl;
      };

      Constructing::Constructing( const std::string& class_name, const std::string& msg ) :
         Base_t( "Error constructing class", class_name, msg ){ /* */ };

      SanityCheckFailure::SanityCheckFailure( const std::string& class_name, const std::string& msg ) :
         Base_t( "Sanity Check failure in class", class_name, msg ){ /* */ };

      UndefinedFunction::UndefinedFunction( const std::string& class_name ) : 
         std::logic_error( std::string("Called undefined member function of base class '") + class_name + "'" ) {
         std::cerr << what() << std::endl;
      };

      NotYetProgrammed::NotYetProgrammed( const std::string& function_name ) : 
         std::logic_error( std::string("Not yet programmed: '") + function_name + "'" ) {
         std::cerr << what() << std::endl;
      };

   };
};

