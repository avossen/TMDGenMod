/*
  Includes ctq6pdf.F interface, so can call PartonX6 directly
  Also initializes CTEQ
  Has optimization if call CTEQ multiple times at same point
*/

#include "DF.d/CTEQ_Wrapper.h"
#include "Common.d/Exceptions.h"
#include <cmath>
#include <string>
#include <iostream>
#include <unistd.h>

using std::cerr;
using std::endl;

extern "C"{

   void gmc_trans_setctq6_( const int* iset, int* ierr );
   double partonx6_( int*, double*, double* );

   extern struct ctqpar2_t {
      int Nx, Nt, NfMx;
   } ctqpar2_;

   extern struct QCDtable_t {
      double Alambda, Nfl, Iorder;
   } qcdtable_;
};



namespace TMDGen {
   namespace CTEQ {

      int& iset( Wrapper_t::iset ); 
      std::string& path( Wrapper_t::path );

      int Wrapper_t::iset = 3;
      std::string Wrapper_t::path = "";

      Wrapper_t::Wrapper_t() {
         cerr << "\tInitializing CTEQ with iset = " << iset << " and path = '" << path << "'" << endl;

         last_iparton = -999;
         last_x = -1;
         last_Q2 = -1;
         last_val = -1;

         int ierr = 0;
         std::string cwd;

         if( path[0] != '\0' ){
            char cwd_temp[1000];

            if( !getcwd( cwd_temp, 1000 ) )
               throw Error::Constructing( "CTEQ::Wrapper_t", "Cannot store current working directory");
            cwd = cwd_temp;

            try{ 
               std::cerr << "\t\tAttempting to change to directory '" << path << "'" << std::endl;
               chdir( path.data() );
            }
            catch( std::exception& e){
               throw Error::Constructing( "CTEQ::Wrapper_t", std::string("Caught error '") + e.what() +  "' while changing directories" );
            }
            catch(...){
               throw Error::Constructing( "CTEQ::Wrapper_t", std::string("error changing to directory '") + path );
            };
         };

         gmc_trans_setctq6_( &iset, &ierr );

         if( ierr )
            throw Error::Constructing( "CTEQ::Wrapper_t", "gmc_trans_setctq6 returned an error" );

         if( !cwd.empty() ){
            try{ 
               std::cerr << "\t\tAttempting to return to directory '" << cwd << "'" << std::endl;
               chdir( cwd.data() );
            }
            catch( std::exception& e){
               throw Error::Constructing( "CTEQ::Wrapper_t", std::string("Caught error '") + e.what() +  "' while changing directories" );
            }
            catch(...){
               throw Error::Constructing( "CTEQ::Wrapper_t", std::string("error changing to directory '") + path );
            };
         };
      };

      Wrapper_t& Wrapper_t::Instance(){
         static Wrapper_t singleton;
         return singleton;
      };

      double Wrapper_t::Eval( int iparton, double x, double Q2 ){
         double val = last_val;

         if( x != last_x || Q2 != last_Q2 || iparton != last_iparton ){
            last_x = x;
            last_Q2 = Q2;
            last_iparton = iparton;
            val = 0;

            //std::cout << iparton << ' ' << ctqpar2_.NfMx << ' ' << sqrt(Q2) << ' ' << qcdtable_.Alambda << endl;
            if( iparton >= -ctqpar2_.NfMx && iparton <= ctqpar2_.NfMx ){
               if( Q2 > 0 ){
                  double Q = sqrt(Q2);
                  if( x < 1 && x > 0 && Q > qcdtable_.Alambda )
                     val = partonx6_( &iparton, &x, &Q );
               };

               if( val < 0 )
                  val = 0;
            };

            last_val = val;
         };

         return val;
      };

   };
};
