/*
   Torino Group's transversity function
   arXiv:hep-ph/0701006v3
   arXiv:0812.4366v1
*/

#include "DF.d/h1_Torino.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"
#include "Common.d/Flavor.h"
#include "Common.d/Consts.h"

#include <cmath>
#include <iostream>
#include <string>

using std::cerr;
using std::endl;

namespace TMDGen {
   namespace DF {

      const double h1_Torino_t::x_min = 1e-4;
      const double h1_Torino_t::x_max = 1;
      const double h1_Torino_t::Q2_min = 0.8;
      const double h1_Torino_t::Q2_max = 1e6;
      const double h1_Torino_t::pT_min = 0;
      const double h1_Torino_t::pT_max = 1e6; // completely arbitrary


      h1_Torino_t::h1_Torino_t( std::string param_code ) :
         FlavArrayFunc_t( std::string("pT NULL")),
         params(param_code),
         f1_GRV(DF_Wrapper::GRV_t::Instance()),
         g1_GRSV(DF_Wrapper::GRSV_t::Instance())
      {
         // make msg
         constr_msg = std::string("h_1: Torino '") + param_code + "'";

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &h1_Torino_t::Gluon );
      };

      h1_Torino_t::~h1_Torino_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void h1_Torino_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         var_min.pT = pT_min;
         var_max.pT = pT_max;
      };

      // flags with 1 each relevant variable used
      void h1_Torino_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double h1_Torino_t::Up( const Var_t& var ) const{
         return inner_function( UP_FLAV, u_idx, var );
      };

      double h1_Torino_t::Down( const Var_t& var ) const{
         return inner_function( DOWN_FLAV, d_idx, var );
      };

      double h1_Torino_t::Strange( const Var_t& var ) const{
         return inner_function( STR_FLAV, s_idx, var );
      };

      double h1_Torino_t::Charm( const Var_t& var ) const{
         return 0;
      };

      double h1_Torino_t::Bottom( const Var_t& var ) const{
         return 0;
      };

      double h1_Torino_t::Anti_Up( const Var_t& var ) const{
         return inner_function( ANTI_UP_FLAV, ubar_idx, var );
      };

      double h1_Torino_t::Anti_Down( const Var_t& var ) const{
         return inner_function( ANTI_DOWN_FLAV, dbar_idx, var );
      };

      double h1_Torino_t::Anti_Strange( const Var_t& var ) const{
         return inner_function( ANTI_STR_FLAV, s_idx, var );
      };

      double h1_Torino_t::Anti_Charm( const Var_t& var ) const{
         return 0;
      };

      double h1_Torino_t::Anti_Bottom( const Var_t& var ) const{
         return 0;
      };

      double h1_Torino_t::Gluon( const Var_t& var ) const{
         return 0;
      };

      double h1_Torino_t::Quark( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

      double h1_Torino_t::Sea( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

      double h1_Torino_t::inner_function( flavor_t flav1, GRV_flavor_t flav2, const Var_t& var ) const {
         double N = params.N[ flav1 ];
         double val = 0;

         if( N ){
            val = f1_GRV.Eval( flav2, var.x, var.Q2);
            val += g1_GRSV.Eval( flav2, var.x, var.Q2);
            if( val ){
               val *= N;
               val *= pow( var.x, params.alpha );
               val *= pow( 1.-var.x, params.beta );
               val *= pow( params.alpha + params.beta, params.alpha + params.beta );
               val /= pow( params.alpha, params.alpha );
               val /= pow( params.beta, params.beta );
               val *= 0.5;
               val *= exp( -var.pT*var.pT / params.ave_pT2 );
               val /= PI;
               val /= params.ave_pT2;
            };
         };

         return val;
      };


   };
};
