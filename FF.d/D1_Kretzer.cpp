/*

--> Ported from the file pkhff.F, which begins with description

***********************************************************************
*                                                                     *
*  LO and NLO FRAGMENTATION FUNCTIONS                                 *
*  for charged pions, kaons and the inclusive sum of charged hadrons  *
*                                                                     *
*  as in S. Kretzer (kretzer@pa.msu.edu):                             *
*                                                                     * 
*  `Fragmentaion Functions from Flavour-inclusive and Flavour-tagged  *
*  e^+ e^- Annihilations';                                            *
*  Phys. Rev. D 62, 054001 (2000)                                     *
*                                                                     *
*  See above reference for details!                                   *  
*								      *	
...
*                                                                     *
*  This interpolation routine returs FFs in the range:                *
*             Z                    (between  0.01   and  1.0)         *
*             Q2 = scale in GeV^2  (between  0.8    and  1.D6)        *
*             Attention: Z \lesssim 0.05   is (very!) delicate        *
*                        ( as discussed in above reference )          *
*                                                                     *
*  The quality of the interpolation (as compared to the exact         * 
*  evolution) is about:                                               *
*  | z < 0.75  | z < 0.9 | z < 1.0 |                                  *                       
*   -------------------------------                                   *
*  |  < 3%     |  < 10%  |  > 10%  |                                  *
*                                                                     *
...
*                                                                     *
***********************************************************************

--> Ported to C++ by S. Gliske, Jan, 2010.

    Actually outputs D1( z, Q2 ) times a kT distribution

*/
 

#include "FF.d/D1_Kretzer.h"

#include "Common.d/Exceptions.h"
#include "pT_distr.d/pT_distr.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using std::cerr;
using std::endl;

#include <cmath>
#include <string>
#include <cstdlib>

namespace TMDGen {
   namespace FF {

      const double D1_Kretzer_t::Q2_min = 0.8;
      const double D1_Kretzer_t::Q2_max = 1e6;
      const double D1_Kretzer_t::z_min = 0.01;
      const double D1_Kretzer_t::z_max = 1.0;

      const int D1_Kretzer_t::N_var = 2;
      const int D1_Kretzer_t::N_per[2] = { D1_KRETZER_NZ, D1_KRETZER_NQ };

      // note: in FORTRAN, Kretzer specific the grid in Z, Q^2, and then
      // took the log when combining, before passing to FINT.
      // Here, the log has been taken already
      const double D1_Kretzer_t::zq_grid[D1_KRETZER_NZ + D1_KRETZER_NQ ] = {
         -4.656463e+00, -4.605170e+00, -3.912023e+00, -3.506558e+00, -3.218876e+00, -2.995732e+00,
         -2.813411e+00, -2.659260e+00, -2.525729e+00, -2.407946e+00, -2.353878e+00, -2.302585e+00,
         -2.079442e+00, -1.897120e+00, -1.742969e+00, -1.609438e+00, -1.491655e+00, -1.386294e+00,
         -1.290984e+00, -1.203973e+00, -1.123930e+00, -1.049822e+00, -9.808293e-01, -9.162907e-01,
         -8.556661e-01, -7.985077e-01, -7.444405e-01, -6.931472e-01, -6.443570e-01, -5.978370e-01,
         -5.533852e-01, -5.108256e-01, -4.700036e-01, -4.307829e-01, -3.930426e-01, -3.566749e-01,
         -3.215836e-01, -2.876821e-01, -2.548922e-01, -2.231436e-01, -1.923719e-01, -1.625189e-01,
         -1.335314e-01, -1.053605e-01, -7.796154e-02, -5.129329e-02, -2.531781e-02, -1.005034e-02,
         0.000000e+00,

         -2.231436e-01,
         0.000000e+00, 2.623643e-01, 5.877867e-01, 9.932518e-01, 1.386294e+00, 1.856298e+00,
         2.302585e+00, 2.772589e+00, 3.218876e+00, 3.688879e+00, 4.158883e+00,
         4.605170e+00, 5.192957e+00, 5.768321e+00, 6.345636e+00,
         6.907755e+00, 7.495542e+00, 8.070906e+00, 8.648221e+00,
         9.210340e+00, 9.998798e+00, 1.073640e+01,
         1.151293e+01, 1.230138e+01, 1.303898e+01,
         1.381551e+01
      };

      // constructor
      // loads grid from file
      D1_Kretzer_t::D1_Kretzer_t( LundPID_t hadron_PID, int order, const char* path, const char* kT_instruction ) : FlavArrayFunc_t( std::string("kT ") + kT_instruction ) {

         std::string filename( path );
         if( !filename.empty() )
            filename += "/";

         // must have shift = 0 when do_ave is !0
         shift = 0;
         do_ave = 0;

         // signs of quark & antiquark factors take care of favored vs. disfavored
         // and factors of 0 for PI_ZERO results in the average of favored and disfavored
         if( order == 0 ){
            switch ( hadron_PID ){
            case PI_PLUS:
               filename += "plo.grid";
               break;
            case PI_MINUS:
               shift = 1;
               filename += "plo.grid";
               break;
            case PI_ZERO:
               do_ave = 1;
               filename += "plo.grid";
               break;
            case K_PLUS:
               filename += "klo.grid";
               break;
            case K_MINUS:
               shift = 1;
               filename += "klo.grid";
               break;
            default:
               throw Error::Constructing( "D1_Kretzer_t", "passed invalid Lund PID" );
            };
         } else if ( order == 1 ){
            switch ( hadron_PID ){
            case PI_PLUS:
               filename += "pnlo.grid";
               break;
            case PI_MINUS:
               shift = 1;
               filename += "pnlo.grid";
               break;
            case PI_ZERO:
               filename += "pnlo.grid";
               break;
            case K_PLUS:
               filename += "knlo.grid";
               break;
            case K_MINUS:
               shift = 1;
               filename += "knlo.grid";
               break;
            default:
               throw Error::Constructing( "D1_Kretzer_t", "passed invalid Lund PID" );
            };
         } else {
            throw Error::Constructing( "D1_Kretzer_t", "Order must be 0 (LO) or 1 (NLO)" );
         };

         // load data from file

         std::ifstream fin( filename.data() );
         if( !fin ){
            throw Error::Constructing( "D1_Kretzer_t", std::string("Cannot open file '") + filename + "'" );
         };

         int nlines = 0;
         std::string line;

         for( int iz = 0; iz < D1_KRETZER_NZ-1 && !fin.eof(); ++iz ){

            // note: since have log's of grid places
            // sqrt( x ) -> exp( 0.5*log(x) );
            double norm    = pow( 1.- exp(zq_grid[iz]), 4 ) * exp( 0.5*zq_grid[iz] );
            double norm_cb = pow( 1.- exp(zq_grid[iz]), 7 ) * exp( 0.3*zq_grid[iz] );
            double norm_gluon = norm_cb * ( 1.- exp(zq_grid[iz]) );

            for( int iq = 0; iq < D1_KRETZER_NQ && !fin.eof(); ++iq, ++nlines ){
               getline( fin, line );
               for( int j = 0; j < D1_KRETZER_NPART; ++j )
                  grid[j][iq*D1_KRETZER_NZ + iz] = atof( line.substr( j*10, 10 ).data() ) / (j == glue_idx ? norm_gluon : ((j==c_idx||j==b_idx) ? norm_cb : norm ));
            };
         };

         if( nlines != 1296 ){
            throw Error::Constructing( "D1_Kretzer_t", std::string("File '") + filename + "' is too short." );
         };

         // set last z to zero
         for( int j = 0; j < D1_KRETZER_NPART; ++j )
            for( int iq = 0; iq < D1_KRETZER_NQ && !fin.eof(); ++iq, ++nlines )
               grid[j][ D1_KRETZER_NZ-1 + iq*D1_KRETZER_NZ ] = 0;

         // set function pointers

         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Kretzer_t::Gluon );

         // make msg
         {
            std::stringstream ss;
            ss << "D1: using Kretzer's model with order == " << order;
            constr_msg = ss.str();
         }
      };

      D1_Kretzer_t::~D1_Kretzer_t(){
         /* nothing to do */
      };


      void D1_Kretzer_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.z = z_min;
         var_max.z = z_max;
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         pT_func->GetRange( var_min.kT, var_max.kT );
      };

      // flags with 1 each relevant variable used
      void D1_Kretzer_t::GetRelevantVar( Var_t& var ) const{
         var.z = 1;
         var.Q2 = 1;
         var.kT = 1;
      };
      

      double D1_Kretzer_t::Up( const Var_t& var ) const {
         return inner_function_1( u_idx+shift, var );
      };

      double D1_Kretzer_t::Anti_Up( const Var_t& var ) const {
         return inner_function_1( u_bar_idx-shift, var );
      };

      double D1_Kretzer_t::Down( const Var_t& var ) const {
         return inner_function_1( d_idx+shift, var );
      };

      double D1_Kretzer_t::Anti_Down( const Var_t& var ) const {
         return inner_function_1( d_bar_idx-shift, var );
      };

      double D1_Kretzer_t::Strange( const Var_t& var ) const {
         return inner_function_1( s_idx+shift, var );
      };

      double D1_Kretzer_t::Anti_Strange( const Var_t& var ) const {
         return inner_function_1( s_bar_idx-shift, var );
      };

      double D1_Kretzer_t::Charm( const Var_t& var ) const {
         return inner_function_2( c_idx, var, 7 );
      };

      double D1_Kretzer_t::Anti_Charm( const Var_t& var ) const {
         return inner_function_2( c_idx, var, 7 );
      };

      double D1_Kretzer_t::Bottom( const Var_t& var ) const {
         return inner_function_2( c_idx, var, 7 );
      };

      double D1_Kretzer_t::Anti_Bottom( const Var_t& var ) const {
         return inner_function_2( c_idx, var, 7 );
      };


      double D1_Kretzer_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double D1_Kretzer_t::Gluon( const Var_t& var ) const {
         return inner_function_2( glue_idx, var, 8 );
      };

      double D1_Kretzer_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };

      // for up, down, and strage and anti-partners
      inline double D1_Kretzer_t::inner_function_1( const int& iflav, const Var_t& var ) const {

         double val = 0;
         if( var.z > z_min && var.z < z_max && var.Q2 > Q2_min && var.Q2 < Q2_max ){

            flavor_identifier_t flav( static_cast< flavor_identifier_t >( iflav ) );   // warning--no sanity check here!

            double log_input[2] = { log( var.z ), log( var.Q2 ) };
            val = fint_( &N_var, &log_input, &N_per, &zq_grid, &grid[ flav ] );

            if( do_ave ){
               val += fint_( &N_var, &log_input, &N_per, &zq_grid, &grid[ flav+1 ] );   // assuming shift is 0 when do_ave is !0
               val *= 0.5;
            };

            double norm = ( 1. - var.z );
            norm *= norm;
            norm *= norm;
            norm *= sqrt( var.z );

            val *= norm;
            val /= var.z;
         };

         return val;
      };


      // for charm and bottom anti-partners
      inline double D1_Kretzer_t::inner_function_2( const flavor_identifier_t& flav, const Var_t& var, int a ) const {

         double val = 0;
         if( var.z > z_min && var.z < z_max && var.Q2 > Q2_min && var.Q2 < Q2_max ){
            double log_input[2] = { log( var.z ), log( var.Q2 ) };
            val = fint_( &N_var, &log_input, &N_per, &zq_grid, &grid[ flav ] );

            val *= pow( 1.-var.z,   a );
            val *= pow(    var.z, 0.3 );
            val /= var.z;
         };

         return val;
      };

   };
};
