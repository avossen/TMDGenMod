// could also load stepsize from TMDGen file


#include "Common.d/ParseInput.h"
#include "FF.d/Spec_Ia_Params.h"
#include "Extra.d/GSL_fMultiMin.hpp"
#include "Opt-D1-Spec-Ia.d/Optimizer.h"

#ifdef USE_OPEN_MP
#include "Opt-D1-Spec-Ia.d/Optimizer_SM.h"
#endif


#include <iostream>
using std::cerr;
using std::endl;

#include <gsl/gsl_vector.h>
#include <sstream>

double chisq_func (const gsl_vector *v, void *params){
   double output = static_cast< TMDGen::Opt_D1_Spec_Ia::Optimizer_t* >( params ) ->
      ComputeChiSq( TMDGen::FF::DiHad_Spec_Ia::Spec_Ia_Params_t ( v->data ) );

   std::cout << "++ chisq = " << output << endl;

   return output;
};



int main( int argc, char* argv[] ){

   if( argc != 10 ){
      cerr << "ERROR: Usage: " << endl;
      cerr << '\t' << argv[0] << " <N_bins> <N_evals_per_bin_factor> <data_hists_filename> <MC_rec_hist_filename> \\" << endl;
      cerr << "\t\t" << "<MC_rec_hist_filename> <TMDGen_filename> <x_bin> <P_hperp_bin> <M_h_bin>" << endl;
      return 127;
   };

   // parse input into meaningful variables

   int N_bins = atoi(argv[1]);
   double N_evals_factor = atof(argv[2]);
   const char* data_filename = argv[3];
   const char* MC_rec_filename = 0;
   const char* MC_4pi_filename = 0;
   if( argv[4][0] != '-' )
      MC_rec_filename = argv[4];
   if( argv[5][0] != '-' )
      MC_4pi_filename = argv[5];
   const char* TMDGen_filename = argv[6];
   std::string bins[3] = { argv[7], argv[8], argv[9] };

#ifdef USE_OPEN_MP
   TMDGen::Opt_D1_Spec_Ia::Optimizer_SM_t optimizer( N_bins, N_evals_factor, data_filename, MC_rec_filename, MC_4pi_filename,
                                                   TMDGen_filename, bins );
#else
   TMDGen::Opt_D1_Spec_Ia::Optimizer_t optimizer( N_bins, N_evals_factor, data_filename, MC_rec_filename, MC_4pi_filename,
                                                   TMDGen_filename, bins );
#endif

   cerr << "Initalizing all threads" << endl;

   if( optimizer.Init() ){
      cerr << "ERROR initializing" << endl;
      return 127;
   };

   gsl_vector *stepsize = gsl_vector_alloc( 11 );
   gsl_vector_set_all( stepsize, 0.1 );
   gsl_vector_set( stepsize, 6, 10 );
   gsl_vector_set( stepsize, 7, 10 );
   gsl_vector_set( stepsize, 8, 10 );
   gsl_vector_set( stepsize, 9, 10 );

   gsl_vector *init_params = gsl_vector_alloc( 11 );

   sgUtil::ParseInputReturn_t parsed_input;
   sgUtil::ParseInput( TMDGen_filename, parsed_input );

   sgUtil::ParseInputReturn_iterator_t it = parsed_input.find("FF_Spec_Ia_Params");
   if( it == parsed_input.end() ){
     gsl_vector_set( init_params, 0, 0.8 );
     gsl_vector_set( init_params, 1, -0.751 );
     gsl_vector_set( init_params, 2, -0.193 );
     gsl_vector_set( init_params, 3, 2.0 );
     gsl_vector_set( init_params, 4, -0.038 );
     gsl_vector_set( init_params, 5, -0.085 );
     gsl_vector_set( init_params, 6, 100 );
     gsl_vector_set( init_params, 7, 100 );
     gsl_vector_set( init_params, 8, 50 );
     gsl_vector_set( init_params, 9, 50 );
     gsl_vector_set( init_params, 10, 2.5 );
   } else {
     std::stringstream ss;
     ss << it->second;
     for( int i=0; i<11; ++i )
       ss >> *gsl_vector_ptr( init_params, i );
   };

   int N_max_iters = 100;
   double threshold = 0.1;

   it = parsed_input.find("Opt_Spec_Ia_Max_Iters");
   if( it != parsed_input.end() ){
     int n = atoi( it->second.data() );

     if( n > 0 && n < 1000000 ){
       N_max_iters = n;
     } else {
       cerr << "Invalid 'Opt_Spec_Ia_Max_Iters' directive.  Using default value of " << N_max_iters << endl;
     };
   };

   it = parsed_input.find("Opt_Spec_Ia_Threshold");
   if( it != parsed_input.end() ){
     double t = atof( it->second.data() );

     if( t >= 0 && t < 1 ){
       threshold = t;
     } else {
       cerr << "Invalid 'Opt_Spec_Ia_Max_Threshold' directive.  Using default value of " << threshold << endl;
     };
   };


   cerr << "Initalizing the GSL solver" << endl;

   TMDGen::Extras::GSL::fMultiMin_t minimizer( 11 );

   gsl_multimin_function gsl_func;
   gsl_func.f = &chisq_func;
   gsl_func.n = 11;
   gsl_func.params = static_cast< void* >( &optimizer );

   minimizer.Init( &gsl_func, init_params, stepsize );


   cerr << "Iterating" << endl;
   int converged = GSL_CONTINUE;
   for ( int i=0; i<N_max_iters && converged == GSL_CONTINUE; ++i ){

      minimizer.Iterate();
      converged = minimizer.CheckConvergence( threshold );

      std::cout << "! " << i << " : ";
      for( int i=0; i<11; i++ )
         std::cout << gsl_vector_get(minimizer.x, i) << ' ';
      std::cout << " | " << minimizer.Residual() << ' ' << converged << std::endl;
   };

   std::cout << "! final ";
   for( int i=0; i<11; i++ )
      std::cout << gsl_vector_get(minimizer.x, i) << ' ';
   std::cout << " | " << minimizer.Residual() << std::endl;



};


