/*
  DF from GRSV2000 wrapper times a p_T distribution
*/

#include "DF.d/g1_GRSV2000_.h"
#include "DF_Wrapper.d/GRSV.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"

// to change directories
#include <unistd.h>

#include <cmath>
#include <iostream>
#include <string>

using std::cerr;
using std::endl;

namespace TMDGen {
   namespace DF {

      const double g1_GRSV2000_t::x_min = 1e-4;
      const double g1_GRSV2000_t::x_max = 1;
      const double g1_GRSV2000_t::Q2_min = 0.8;
      const double g1_GRSV2000_t::Q2_max = 1e6;

      g1_GRSV2000_t::g1_GRSV2000_t( std::string pT_instructions ) :
         FlavArrayFunc_t( std::string("pT ") + pT_instructions ), GRSV(DF_Wrapper::GRSV_t::Instance())
      {
         // make msg
         constr_msg = "f1: using GRSV98";

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV2000_t::Gluon );

      };

      g1_GRSV2000_t::~g1_GRSV2000_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void g1_GRSV2000_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         pT_func->GetRange( var_min.pT, var_max.pT );
      };

      // flags with 1 each relevant variable used
      void g1_GRSV2000_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double g1_GRSV2000_t::Up( const Var_t& var ) const{
         return GRSV.Eval( u_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Down( const Var_t& var ) const{
         return GRSV.Eval( d_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Strange( const Var_t& var ) const{
         return GRSV.Eval( s_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Charm( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV2000_t::Bottom( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV2000_t::Anti_Up( const Var_t& var ) const{
         return GRSV.Eval( ubar_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Anti_Down( const Var_t& var ) const{
         return GRSV.Eval( dbar_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Anti_Strange( const Var_t& var ) const{
         return GRSV.Eval( s_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Anti_Charm( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV2000_t::Anti_Bottom( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV2000_t::Gluon( const Var_t& var ) const{
         return GRSV.Eval( glue_idx, var.x, var.Q2 ) / var.x;
      };

      double g1_GRSV2000_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double g1_GRSV2000_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };

   };
};
