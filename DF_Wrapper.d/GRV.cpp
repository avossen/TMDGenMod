/*
  Wrapper for grv98.F
  Has optimization if call GRV multiple times at same point
  Note: is 'x f_1(x)'
*/

#include "DF_Wrapper.d/GRV.h"
#include "Common.d/Exceptions.h"
#include <cmath>
#include <string>
#include <iostream>
#include <unistd.h>

using std::cerr;
using std::endl;

extern "C"{
   void grv98_( const int *idoinit, const int *iset, const double *x, const double *Q2, const int *iflav, double *val, int *ierr );
};

namespace TMDGen {
   namespace DF_Wrapper {

      int GRV_t::iset = 1;
      std::string GRV_t::path = "";

      GRV_t::GRV_t() : Wrapper_t() {
         cerr << "\tInitializing GRV with iset = " << iset << " and path = '" << path << "'" << endl;

         if( Try_Move_to_Dir( path, GRV_GRID_DIR, "grv98nld.grid" ) )
            throw Error::Constructing( "DF_Wrapper::GRV_t", std::string("error finding grid file '") + "grv98nld.grid" + "'" );

         int ierr = 0;
         double blah = 0;
         int iblah = 0;
         int doinit = 1;
         grv98_( &doinit, &iset, &blah, &blah, &iblah, &blah, &ierr );

         if( ierr )
            throw Error::Constructing( "DF_Wrapper::GRV_t", "error initializing GRV" );

         if( Return_from_Dir() )
            throw Error::Constructing( "DF_Wrapper::GRV_t", "error changing directories" );
      };

      GRV_t& GRV_t::Instance(){
         static GRV_t singleton;
         return singleton;
      };

      double GRV_t::Child_Eval( int iparton, double x, double Q2 ){
         int doinit = 0;
         double val = 0;
         int ierr = 0;
         grv98_( &doinit, &doinit, &x, &Q2, &iparton, &val, &ierr );

         if( ierr )
            val = 0;

         return val;
      };

   };
};
