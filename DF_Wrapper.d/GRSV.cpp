/*
  Wrapper for grsv2000.F
  Has optimization if call GRSV multiple times at same point
  Note: is 'x g_1(x)'
*/

#include "DF_Wrapper.d/GRSV.h"
#include "Common.d/Exceptions.h"
#include <cmath>
#include <string>
#include <iostream>
#include <unistd.h>

using std::cerr;
using std::endl;

extern "C"{
   void grsv2000_( const int *idoinit, const int *iset, const double *x, const double *Q2, const int *iflav, double *val, int *ierr );
};

namespace TMDGen {
   namespace DF_Wrapper {

      int GRSV_t::iset = 1;
      std::string GRSV_t::path = "";

      GRSV_t::GRSV_t() : Wrapper_t() {
         cerr << "\tInitialized GRSV with iset = " << iset << " and path = '" << path << "'" << endl;

         if( Try_Move_to_Dir( path, GRV_GRID_DIR, "std2000_lo_g1.grid" ) )
            throw Error::Constructing( "DF_Wrapper::GRSV_t", std::string("error finding grid file '") + "std2000_lo_g1.grid" + "'" );

         int ierr = 0;
         double blah = 0;
         int iblah = 0;
         int doinit = 1;
         grsv2000_( &doinit, &iset, &blah, &blah, &iblah, &blah, &ierr );

         if( ierr )
            throw Error::Constructing( "DF_Wrapper::GRSV_t", "error initializing GRSV" );

         if( Return_from_Dir() )
            throw Error::Constructing( "DF_Wrapper::GRSV_t", "error changing directories" );
      };

      GRSV_t& GRSV_t::Instance(){
         static GRSV_t singleton;
         return singleton;
      };

      double GRSV_t::Child_Eval( int iparton, double x, double Q2 ){
         int doinit = 0;
         double val = 0;
         int ierr = 0;
         grsv2000_( &doinit, &doinit, &x, &Q2, &iparton, &val, &ierr );

         if( ierr )
            val = 0;

         return val;
      };

   };
};
