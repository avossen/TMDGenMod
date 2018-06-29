/*
  Generic Fortran wrapper for DFs
*/

#ifndef _DF_WRAPPER_H_
#define _DF_WRAPPER_H_

#include <string>

namespace TMDGen {
   namespace DF_Wrapper {

      class Wrapper_t {
      protected:
         double last_iparton;
         double last_x;
         double last_Q2;
         double last_val;

         std::string cwd;

         Wrapper_t();
         int Move_to_Dir( std::string path );
         int Return_from_Dir();
         int File_Exists( std::string path, std::string file ) const;
         int Try_Move_to_Dir( std::string path1, std::string path2, std::string file );

         virtual double Child_Eval( int iparton, double x, double Q2 ) = 0;

      public:
         double Eval( int iparton, double x, double Q2 );
         virtual ~Wrapper_t() { /* */ };
      };
   };
};

#endif
