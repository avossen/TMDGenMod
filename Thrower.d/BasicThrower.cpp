/*
     Throws variables according to continuous distribution per flavor
     Most similar to old FORTRAN w/ weights option
*/

#include "Thrower.d/BasicThrower.h"
#include "Thrower.d/GSL_Integration.h"
#include "Common.d/Consts.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_monte_vegas.h>

#include <iostream>
using std::cerr;
using std::endl;

namespace TMDGen {
   namespace Thrower {

      BasicThrower_t::BasicThrower_t( RNG_t *r_in, int N_integration_warmup_,
                                      int N_integration_calls_, int N_max_integration_calls_,
                                      const Var_t& min_in, const Var_t& max_in ) :
         VarThrower_t(r_in), N_integration_warmup( N_integration_warmup_ ), N_integration_calls(N_integration_calls_),
         N_max_integration_calls(N_max_integration_calls_), min(min_in), max(max_in), width("BasicThrower_t::width")
      {
         //
      };

      BasicThrower_t::~BasicThrower_t(){
         /* nothing to do */
      };

      int BasicThrower_t::Initialize( Var_t& var, const XSec::XSec_t& xsec ){
         // set up random number generator
         gsl_rng* r_gsl = gsl_rng_alloc ( gsl_rng_ranlxs2 );
         gsl_rng_set( r_gsl, static_cast< int >( r->EvalUnif()*64000+10 ) ); // seed with some number between 10 and 64009

         // gsl monte function
         gsl_monte_function gsl_F;
         GSL_params_t gsl_params( var );
         gsl_params.xsec = &xsec;
         gsl_F.params = &gsl_params;

         int ierr = 0;

         double *min_array=0, *max_array=0;

         Init_GSL_Func( &gsl_F, min_array, max_array );

         gsl_monte_vegas_state *gsl_state = gsl_monte_vegas_alloc( gsl_F.dim );
         gsl_state->stage = 0;

         // display before integrating
         cerr << "\tIntegrating cross section per flavor" << endl;
         cerr << endl << "\t\tflavor\tvalue\t\t\t\tierr code N_evaluations\tVEGAS \\chi^2" << endl;
         cerr <<         "\t\t------\t-----\t\t\t\t--------- -------------\t-------------" << endl;

         total_xsec = 0;
         for( int flavor = 0; flavor < GMC_TRANS_N_FLAVORS; ++flavor ){
            gsl_params.var.flavor = static_cast< flavor_t >( flavor );

            if( flavor == GLUON_FLAV ){
               // cerr << "\tWARNING: Artificially setting CX for to zero for " << flavor_string[ flavor ] << "." << endl;
               // set gluon probability to one (as do not have FF for gluons)
               integrated_xsec[flavor] = 0;
               integrated_xsec_abserr[flavor] = 0;
            } else {
               //cerr << "\tIntegrating flavor " << flavor_string[flavor] << endl;

               // initialize
               gsl_monte_vegas_init( gsl_state );

               gsl_state->stage = 0;
               ierr = gsl_monte_vegas_integrate( &gsl_F, min_array, max_array, gsl_F.dim, N_integration_warmup, r_gsl, gsl_state, &integrated_xsec[flavor], &integrated_xsec_abserr[flavor] );
               int i = N_integration_warmup;

               // ensure is not identically zero
               if( integrated_xsec_abserr[flavor] ){

                  gsl_state->stage = 1;
                  do {
                     ierr = gsl_monte_vegas_integrate( &gsl_F, min_array, max_array, gsl_F.dim, N_integration_warmup, r_gsl, gsl_state, &integrated_xsec[flavor], &integrated_xsec_abserr[flavor] );
                     i += N_integration_calls;

                     gsl_state->stage = 3;

                     // break if too many calls, if chisq is small enough, or if function is still zero
                  } while ( i<N_max_integration_calls                                                  // keep going if not too many calls
                            && (fabs (gsl_state->chisq - 1.0) > 0.25)                                  // and not converged
                            && ( integrated_xsec[flavor] || integrated_xsec_abserr[flavor] ) );        // and not identically zero nor constant

                  total_xsec += integrated_xsec[flavor];
               };

               cerr << "\t\t" << flavor << ' ' << flavor_string[ flavor ] << ":\t" << integrated_xsec[flavor];
               cerr << " +/- " << integrated_xsec_abserr[flavor];
               cerr << "\t\t" << ierr << ' ' << i << '\t' << gsl_state->chisq << endl;

            };
         };

         cerr << endl << "\tTotal cross section = " << total_xsec << endl;
         cerr << endl << "\tBranching Ratios" << endl;
         for( int flavor = 0; flavor < GMC_TRANS_N_FLAVORS-1; ++flavor ){
            cerr << "\t\t" << flavor_string[ flavor ] << ":\t" << integrated_xsec[flavor]/total_xsec;
            if( integrated_xsec[flavor] )
               cerr << "\t+/- " << integrated_xsec_abserr[flavor]/total_xsec;
            cerr << endl;
         };

         // free things
         gsl_rng_free( r_gsl );
         gsl_monte_vegas_free( gsl_state );

         delete[] min_array;
         delete[] max_array;

         return ierr;
      };

      void BasicThrower_t::Throw( Var_t& var, double& pdf_val ){
         ThrowVariables( var, pdf_val );
	 ThrowFlavor( var, pdf_val );
      };

      void BasicThrower_t::ThrowFlavor( Var_t& var, double& pdf_val ){
         // determine quark flavor
         int iflavor = 0;
         {
            double r1 = r->EvalUnif()*total_xsec;
            double r2 = integrated_xsec[iflavor];

            while( r1 > r2 && iflavor < GMC_TRANS_N_FLAVORS )
               r2 += integrated_xsec[++iflavor];
         };

         // cast int to enum and save
         var.flavor = static_cast< flavor_t >( iflavor );

         // update pdf_val
         pdf_val *= integrated_xsec[ var.flavor ]/total_xsec;
      };

      void BasicThrower_t::Throw_pT( Var_t& var, double& pdf_val ){
         pdf_val /= (max.pT      - min.pT);
         pdf_val /= TWO_PI;

         var.pT      = r->EvalUnif()*width.pT      + min.pT;
         var.phi_pT  = r->EvalUnif()*TWO_PI;
      };


   };
};
