/*
  Includes ctq6pdf.F interface, so can call PartonX6 directly
  Also initializes CTEQ
  Has optimization if call CTEQ multiple times at same point
*/

#include "DF_Wrapper.d/CTEQ.h"
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
   namespace DF_Wrapper {

      int CTEQ_t::iset = 3;
      std::string CTEQ_t::path = "";

      CTEQ_t::CTEQ_t() : Wrapper_t() {
         cerr << "\tInitialized CTEQ with iset = " << iset << " and path = '" << path << "'" << endl;

         if( Try_Move_to_Dir( path, CTEQ_GRID_DIR, "pdf/cteq6d.tbl" ) )
            throw Error::Constructing( "DF_Wrapper::CTEQ_t", std::string("error finding grid file '") + "pdf/cteq6d.tbl" + "'" );

         int ierr = 0;
         gmc_trans_setctq6_( &iset, &ierr );

         if( ierr )
            throw Error::Constructing( "DF_Wrapper::CTEQ_t", "error initializing CTEQ" );

         if( Return_from_Dir() )
            throw Error::Constructing( "DF_Wrapper::CTEQ_t", "error changing directories" );

      };

      CTEQ_t& CTEQ_t::Instance(){
         static CTEQ_t singleton;
         return singleton;
      };

      double CTEQ_t::Child_Eval( int iparton, double x, double Q2 ){
         double val = 0;

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

         return val;
      };

   };
};
