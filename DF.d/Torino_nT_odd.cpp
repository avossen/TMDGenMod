/*
   Torino Group's transversity function
   arXiv:hep-ph/0701006v3
   arXiv:0812.4366v1
*/

#include "DF.d/Torino_nT_odd.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"
#include "Common.d/Flavor.h"
#include "Common.d/Consts.h"

#include <cmath>
#include <iostream>
#include <string>

using std::cout;
using std::cerr;
using std::endl;

namespace TMDGen {
   namespace DF {

      const double Torino_nT_odd_t::x_min = 1e-4;
      const double Torino_nT_odd_t::x_max = 1;
      const double Torino_nT_odd_t::Q2_min = 0.8;
      const double Torino_nT_odd_t::Q2_max = 1e6;
      const double Torino_nT_odd_t::pT_min = 0;
      const double Torino_nT_odd_t::pT_max = 1e6; // completely arbitrary


      Torino_nT_odd_t::Torino_nT_odd_t( std::string DF_name, std::string param_code ) :
         FlavArrayFunc_t( std::string("pT NULL")),
         params( DF_name + "_" + param_code),
         f1_GRV(DF_Wrapper::GRV_t::Instance())
      {
         // make msg
         constr_msg = DF_name + std::string(": Torino '") + param_code + "'";

         if( DF_name != "h_1^perp" && DF_name != "f_1T^perp" )
            throw Error::Constructing("Torino_nT_odd_t", std::string("Invalid DF name: '") + DF_name + "'" );

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Torino_nT_odd_t::Gluon );

         slope = -(1./params.ave_pT_sq + 1./params.M1_sq);
         //cerr << "slope = " << slope << ", " << params.ave_pT_sq << ' ' << params.M1_sq << endl; 

         prefactor = sqrt(2*exp(1)) * PROTON_MASS/sqrt(params.M1_sq) / PI / params.ave_pT_sq;
      };

      Torino_nT_odd_t::~Torino_nT_odd_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void Torino_nT_odd_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         var_min.pT = pT_min;
         var_max.pT = pT_max;
      };

      // flags with 1 each relevant variable used
      void Torino_nT_odd_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double Torino_nT_odd_t::Up( const Var_t& var ) const{
         return inner_function( UP_FLAV, u_idx, var );
      };

      double Torino_nT_odd_t::Down( const Var_t& var ) const{
         return inner_function( DOWN_FLAV, d_idx, var );
      };

      double Torino_nT_odd_t::Strange( const Var_t& var ) const{
         return inner_function( STR_FLAV, s_idx, var );
      };

      double Torino_nT_odd_t::Charm( const Var_t& var ) const{
         return 0;
      };

      double Torino_nT_odd_t::Bottom( const Var_t& var ) const{
         return 0;
      };

      double Torino_nT_odd_t::Anti_Up( const Var_t& var ) const{
         return inner_function( ANTI_UP_FLAV, ubar_idx, var );
      };

      double Torino_nT_odd_t::Anti_Down( const Var_t& var ) const{
         return inner_function( ANTI_DOWN_FLAV, dbar_idx, var );
      };

      double Torino_nT_odd_t::Anti_Strange( const Var_t& var ) const{
         return inner_function( ANTI_STR_FLAV, s_idx, var );
      };

      double Torino_nT_odd_t::Anti_Charm( const Var_t& var ) const{
         return 0;
      };

      double Torino_nT_odd_t::Anti_Bottom( const Var_t& var ) const{
         return 0;
      };

      double Torino_nT_odd_t::Gluon( const Var_t& var ) const{
         return 0;
      };

      double Torino_nT_odd_t::Quark( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

      double Torino_nT_odd_t::Sea( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

      double Torino_nT_odd_t::inner_function( flavor_t flav1, GRV_flavor_t flav2, const Var_t& var ) const {
         double A = params.A[ flav1 ];
         double val = 0;

         if( A ){
            val = f1_GRV.Eval( flav2, var.x, var.Q2);
            if( val ){
               val *= A;
               val *= pow( var.x, params.alpha[ flav1 ] );
               val *= pow( 1.-var.x, params.beta );
               val *= pow( params.alpha[ flav1 ] + params.beta, params.alpha[ flav1 ] + params.beta );
               val /= pow( params.alpha[ flav1 ], params.alpha[ flav1 ] );
               val /= pow( params.beta, params.beta );
               val *= exp( var.pT*var.pT*slope );
               val *= prefactor;
            };

            //             cout << constr_msg << "\t";
//             cout << f1_GRV.Eval( flav2, var.x, var.Q2) << ' ';
//             cout << ' ' << pow( var.x, params.alpha[ flav1 ] );
//             cout << ' ' << pow( 1.-var.x, params.beta );
//             cout << ' ' << pow( params.alpha[ flav1 ] + params.beta, params.alpha[ flav1 ] + params.beta );
//             cout << ' ' << pow( params.alpha[ flav1 ], params.alpha[ flav1 ] );
//             cout << ' ' << pow( params.beta, params.beta );
//             cout << ' ' << exp( var.pT*var.pT*slope );
//             cout << ' ' << prefactor;
//             cout << endl;

//             val = ( ( params.lambda[ flav1 ] < -999 ) ? -fabs(val) : params.lambda[ flav1 ]*val );
//             cout << val << endl;
         };

         return val;
      };


   };
};
