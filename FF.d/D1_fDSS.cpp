/*

--> From the file fDSS.F, which begins with the description

********************************************************************
*                                                                  *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
*  D.de Florian, R.Sassot, M.Stratmann   Phys.Rev.D75 114010 2007  *
*                                 *and*  Phys.Rev.D76 074033 2007  *
*                                                                  *
...
*                                                                  *
*            X                    (between  0.05   and  1.0)       *
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
*                                                                  *
...
********************************************************************

--> Ported to C++ by S. Gliske, Jan, 2010.

    Note: some variable names follow the fortran code,
    e.g. sometimes z is denoted x

    fDSS originally returned ``z D1( Q^2, z )''
    Now modified to return ``D1(Q^2, z)'', i.e. no z prefactor

    Note also returns product of D1 with a kT distribution

*/
 

#include "FF.d/D1_fDSS.h"

#include "Common.d/Exceptions.h"
#include "pT_distr.d/pT_distr.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using std::cerr;
using std::endl;

#include <cmath>
#include <cstdlib>

namespace TMDGen {
   namespace FF {

      const double D1_fDSS_t::Q2_min = 1.0;
      const double D1_fDSS_t::Q2_max = 1e5;
      const double D1_fDSS_t::z_min = 0.05;
      const double D1_fDSS_t::z_max = 1.0;

      const int D1_fDSS_t::N_var = 2;
      const int D1_fDSS_t::N_per[2] = { D1_DFSS_NX, D1_DFSS_NQ };

      const double D1_fDSS_t::xq_grid[D1_DFSS_NX + D1_DFSS_NQ ] = { 0.01,  0.02,  0.03,  0.04,  0.05,  0.06,  0.07,
                                                                    0.08,  0.09,  0.095, 0.1,   0.125, 0.15,  0.175,
                                                                    0.2,   0.225, 0.25,  0.275, 0.3,   0.325, 0.35,
                                                                    0.375, 0.4,   0.45,  0.5,   0.55,  0.6,   0.65,
                                                                    0.7,   0.75,  0.8,   0.85,  0.9,   0.93,  1.0, 
                                                                    1.e0, 1.25e0, 1.5e0, 2.5e0, 4.0e0, 6.4e0,
                                                                    1.0e1, 1.5e1, 2.5e1, 4.0e1, 6.4e1, 1.0e2,
                                                                    1.8e2, 3.2e2, 5.8e2, 1.0e3, 1.8e3, 3.2e3,
                                                                    5.8e3, 1.0e4, 1.8e4, 3.2e4, 5.8e4, 1.0e5   }; 

      const double D1_fDSS_t::log_xq_grid[D1_DFSS_NX + D1_DFSS_NQ ] = { -4.60517019, -3.91202301, -3.50655790, -3.21887582, -2.99573227, -2.81341072, -2.65926004,
                                                                        -2.52572864, -2.40794561, -2.35387839, -2.30258509, -2.07944154, -1.89711998, -1.74296931,
                                                                        -1.60943791, -1.49165488, -1.38629436, -1.29098418, -1.20397280, -1.12393010, -1.04982212,
                                                                        -0.98082925, -0.91629073, -0.79850770, -0.69314718, -0.59783700, -0.51082562, -0.43078292,
                                                                        -0.35667494, -0.28768207, -0.22314355, -0.16251893, -0.10536052, -0.07257069,  0.00000000,
                                                                        0.00000000,  0.22314355,  0.40546511,  0.91629073,  1.38629436,  1.85629799,
                                                                        2.30258509,  2.70805020,  3.21887582,  3.68887945,  4.15888308,  4.60517019,
                                                                        5.19295685,  5.76832100,  6.36302810,  6.90775528,  7.49554194,  8.07090609,
                                                                        8.66561320,  9.21034037,  9.79812704, 10.37349118, 10.96819829, 11.51292546 };

      // constructor
      // loads grid from file
      D1_fDSS_t::D1_fDSS_t( LundPID_t hadron_PID, int order, const char* path, const char* kT_instruction ) : FlavArrayFunc_t( std::string("kT ") + kT_instruction ) {

         //for( int i=0; i<D1_DFSS_NX + D1_DFSS_NQ; i++ )
         //   std::cout << i+1 << ' ' << log_xq_grid[i] << endl;

         std::string filename( path );
         if( !filename.empty() )
            filename += "/";

         // signs of quark & antiquark factors take care of favored vs. disfavored
         // and factors of 0 for PI_ZERO results in the average of favored and disfavored
         if( order == 0 ){
            switch ( hadron_PID ){
            case PI_PLUS:
               quark_factor = 1;
               antiquark_factor = -1;
               filename += "PILO.GRID";
               break;
            case PI_MINUS:
               quark_factor = -1;
               antiquark_factor = 1;
               filename += "PILO.GRID";
               break;
            case PI_ZERO:
               quark_factor = 0;
               antiquark_factor = 0;
               filename += "PILO.GRID";
               break;
            case K_PLUS:
               quark_factor = 1;
               antiquark_factor = -1;
               filename += "KALO.GRID";
               break;
            case K_MINUS:
               quark_factor = -1;
               antiquark_factor = 1;
               filename += "KALO.GRID";
               break;
            default:
               throw Error::Constructing( "D1_fDSS_t", "passed invalid Lund PID" );
            };
         } else if ( order == 1 ){
            switch ( hadron_PID ){
            case PI_PLUS:
               quark_factor = 1;
               antiquark_factor = -1;
               filename += "PINLO.GRID";
               break;
            case PI_MINUS:
               quark_factor = -1;
               antiquark_factor = 1;
               filename += "PINLO.GRID";
               break;
            case PI_ZERO:
               quark_factor = 0;
               antiquark_factor = 0;
               filename += "PINLO.GRID";
               break;
            case K_PLUS:
               quark_factor = 1;
               antiquark_factor = -1;
               filename += "KANLO.GRID";
               break;
            case K_MINUS:
               quark_factor = -1;
               antiquark_factor = 1;
               filename += "KANLO.GRID";
               break;
            default:
               throw Error::Constructing( "D1_fDSS_t", "passed invalid Lund PID" );
            };
         } else {
            throw Error::Constructing( "D1_fDSS_t", "Order must be 0 (LO) or 1 (NLO)" );
         };

         // load data from file

         std::ifstream fin( filename.data() );
         if( !fin ){
            throw Error::Constructing( "D1_fDSS_t", std::string("Cannot open file '") + filename + "'" );
         };

         int nlines = 0;
         std::string line;

         for( int ix = 0; ix < D1_DFSS_NX-1 && !fin.eof(); ++ix ){
            double norm = pow( 1.- xq_grid[ix], 4 ) * sqrt( xq_grid[ix] );
            double norm_cb = pow( 1.- xq_grid[ix], 7 ) * pow( xq_grid[ix], 0.3 );

            for( int iq = 0; iq < D1_DFSS_NQ && !fin.eof(); ++iq, ++nlines ){
               getline( fin, line );
               for( int j = 0; j < D1_DFSS_NPART; ++j ){
                  grid[j][ iq*D1_DFSS_NX + ix ] = atof( line.substr( j*10, 10 ).data() ) / ( j==3||j==4 ? norm_cb : norm );
               };
            };
         };

         if( nlines != 816 ){
            throw Error::Constructing( "D1_fDSS_t", std::string("File '") + filename + "' is too short." );
         };

         // set last x to zero
         for( int j = 0; j < D1_DFSS_NPART; ++j )
            for( int iq = 0; iq < D1_DFSS_NQ && !fin.eof(); ++iq, ++nlines )
               grid[j][ D1_DFSS_NX-1 + iq*D1_DFSS_NX ] = 0;

         // set function pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_fDSS_t::Gluon );


         // make msg
         {
            std::stringstream ss;

            ss << "D1: using fDSS with order == " << order;

            constr_msg = ss.str();
         }
      };

      D1_fDSS_t::~D1_fDSS_t(){
         /* nothing to do */
      };


      void D1_fDSS_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.z = z_min;
         var_max.z = z_max;
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         pT_func->GetRange( var_min.kT, var_max.kT );
      };

      // flags with 1 each relevant variable used
      void D1_fDSS_t::GetRelevantVar( Var_t& var ) const{
         var.z = 1;
         var.Q2 = 1;
         var.kT = 1;
      };
      

      double D1_fDSS_t::Up( const Var_t& var ) const {
         //q std::cout << "u    ";
         return inner_function_1( var, quark_factor, u_idx, u_val_idx );
      };

      double D1_fDSS_t::Anti_Up( const Var_t& var ) const {
         //q std::cout << "ubar ";
         return inner_function_1( var, antiquark_factor, u_idx, u_val_idx );
      };

      double D1_fDSS_t::Down( const Var_t& var ) const {
         //q std::cout << "d    ";
         return inner_function_1( var, quark_factor, d_idx, d_val_idx );
      };

      double D1_fDSS_t::Anti_Down( const Var_t& var ) const {
         //q std::cout << "dbar ";
         return inner_function_1( var, antiquark_factor, d_idx, d_val_idx );
      };

      double D1_fDSS_t::Strange( const Var_t& var ) const {
         //q std::cout << "s    ";
         return inner_function_1( var, quark_factor, s_idx, s_val_idx );
      };

      double D1_fDSS_t::Anti_Strange( const Var_t& var ) const {
         //q std::cout << "sbar ";
         return inner_function_1( var, antiquark_factor, s_idx, s_val_idx );
      };

      double D1_fDSS_t::Charm( const Var_t& var ) const {
         return inner_function_2( var, c_idx );
      };

      double D1_fDSS_t::Anti_Charm( const Var_t& var ) const {
         return inner_function_2( var, c_idx );
      };

      double D1_fDSS_t::Bottom( const Var_t& var ) const {
         return inner_function_2( var, b_idx );
      };

      double D1_fDSS_t::Anti_Bottom( const Var_t& var ) const {
         return inner_function_2( var, b_idx );
      };

      double D1_fDSS_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double D1_fDSS_t::Gluon( const Var_t& var ) const {
         return Undefined( var );
      };

      double D1_fDSS_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };


      inline double D1_fDSS_t::inner_function_1( const Var_t& var, double factor, flavor_identifier_t flav1, flavor_identifier_t flav2 ) const {
         double output = 0;
         if( var.z > z_min && var.z < z_max && var.Q2 > Q2_min && var.Q2 < Q2_max ){
            double log_input[2] = { log(var.z), log(var.Q2) };

            // call the fortran interpolation program
            output = fint_( &N_var, &log_input, &N_per, &log_xq_grid, &grid[ flav1 ] );
            if( factor ){
               output += factor*fint_( &N_var, &log_input, &N_per, &log_xq_grid, &grid[ flav2 ] );
               output *= 0.5;
            };

            double norm = ( 1 - var.z );
            norm *= norm;
            norm *= norm;
            norm *= sqrt( var.z );

            output *= norm;

            //            for( int i=110; i<120; ++i )
            //   std::cout << ".c " << grid[flav1][i] << endl;
            /*
            double a = norm * fint_( &N_var, &log_input, &N_per, &log_xq_grid, &grid[ flav1 ] );
            double b = norm * fint_( &N_var, &log_input, &N_per, &log_xq_grid, &grid[ flav2 ] );

            std::cout << " c " << var.z << ' ' << a << ' ' << b << endl;

            output = (a+b)/2;
            */
         };


         return output;
      };


      inline double D1_fDSS_t::inner_function_2( const Var_t& var, flavor_identifier_t flav1 ) const {
         double output = 0;
         if( var.z > z_min && var.z < z_max && var.Q2 > Q2_min && var.Q2 < Q2_max ){
            double log_input[2] = { log(var.z), log(var.Q2) };

            // call the fortran interpolation program
            output = fint_( &N_var, &log_input, &N_per, &xq_grid, &grid[ flav1 ] );
            output *= 0.5;

            double norm = ( 1 - var.z );
            norm *= norm;
            norm *= norm;
            norm *= sqrt( var.z );

            output *= norm;
         };

         return output / var.z;
      };
   };
};
