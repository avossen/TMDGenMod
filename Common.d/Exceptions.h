
#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include <stdexcept>
#include <string>

namespace TMDGen {
   namespace Error {

      class Base_t : public std::runtime_error {
      public:
         Base_t( const std::string& prefix, const std::string& class_name, const std::string& msg );
      };

      class Constructing : public Base_t {
      public:
         Constructing( const std::string& class_name, const std::string& msg );
      };

      class SanityCheckFailure : public Base_t {
      public:
         SanityCheckFailure( const std::string& class_name, const std::string& msg );
      };

      class UndefinedFunction : public std::logic_error {
      public:
         UndefinedFunction( const std::string& class_name ); 
      };

      class NotYetProgrammed : public std::logic_error {
      public:
         NotYetProgrammed( const std::string& function_name ); 
      };

   };
};


#endif
