/*
  Generic Fortran wrapper for DFs
*/

#include "DF_Wrapper.d/Wrapper.h"
#include "Common.d/Exceptions.h"
#include <cmath>
#include <string>
#include <iostream>
#include <unistd.h>
#include <sys/stat.h> 

using std::cerr;
using std::endl;

namespace TMDGen {
   namespace DF_Wrapper {

      Wrapper_t::Wrapper_t() {
         last_iparton = -999;
         last_x = -1;
         last_Q2 = -1;
         last_val = -1;
      };

      int Wrapper_t::Try_Move_to_Dir( std::string path1, std::string path2, std::string file ){
         int ierr = 0;
         int found = 0;

         if( !path1.empty() ){
            cerr << "\t\tChecking path '" << path1 << "' for file '" << file << "'" << endl;

            found = File_Exists( path1, file );

            if( found )
               ierr = Move_to_Dir( path1 );
//             else 
//                cerr << "\t\tError using given path. Instead trying '" << path2 << "'" << endl;
         } else {
            cerr << "\t\tGiven path empty, checking current dir" << endl;

            found = File_Exists( ".", file );

            // no need to move to dir if in current path
         };

         if( !found ){
            cerr << "\t\tChecking path '" << path2 << "' for file '" << file << "'" << endl;
            found = File_Exists( path2, file );

            if( found )
               ierr = Move_to_Dir( path2 );
         };

         return ierr;
      };

      int Wrapper_t::File_Exists( std::string path, std::string file ) const{
         struct stat stFileInfo;

         return !stat( (path + "/" + file).c_str(), &stFileInfo);
      };

      int Wrapper_t::Move_to_Dir( std::string path ){
         // change to directory
         if( path[0] != '\0' ){
            char cwd_temp[1000];

            if( !getcwd( cwd_temp, 1000 ) ){
               cerr << "DF_Wrapper::Wrapper_t: Cannot store current working directory" << endl;
               return 1;
            };

            // save current working directory
            cwd = cwd_temp;

            try{ 
               std::cerr << "\t\tDF_Wrapper::Wrapper_t: Attempting to change to directory '" << path << "'" << std::endl;
               chdir( path.data() );
            }
            catch( std::exception& e){
               cerr << "DF_Wrapper::Wrapper_t:  Caught error '" << e.what() << "' while changing directories" << endl;
               return 2;
            }
            catch(...){
               cerr << "DF_Wrapper::Wrapper_t: error changing to directory '" << path << "'" << endl;
               throw;
            };
         };
         return 0;
      };

      int Wrapper_t::Return_from_Dir(){
         // change back from directory
         if( !cwd.empty() ){
            try{ 
               std::cerr << "\t\tDF_Wrapper::Wrapper_t: Attempting to return to directory '" << cwd << "'" << std::endl;
               chdir( cwd.data() );
            }
            catch( std::exception& e){
               cerr << "DF_Wrapper::Wrapper_t: Caught error '" << e.what() << "' while changing directories" << endl;
               return 1;
            }
            catch(...){
               cerr << "DF_Wrapper::Wrapper_t: error changing to directory '" + cwd + "'" << endl;
               throw;
            };
         };
         cwd.clear();
         return 0;
      };

      double Wrapper_t::Eval( int iparton, double x, double Q2 ){
         double val = last_val;

         if( x != last_x || Q2 != last_Q2 || iparton != last_iparton ){
            last_x = x;
            last_Q2 = Q2;
            last_iparton = iparton;
            val = Child_Eval( iparton, x, Q2 );
            last_val = val;
         };

         return val;
      };

   };
};
