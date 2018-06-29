/*
  Wrapper for CTEQ 6 times a p_T distribution
*/

#include "DF.d/CTEQ_DF.h"
#include "DF_Wrapper.d/CTEQ.h"
#include "pT_distr.d/pT_distr.h"
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

      const double CTEQ_DF_t::x_min = 1e-5;
      const double CTEQ_DF_t::x_max = 1;
      const double CTEQ_DF_t::Q2_min = 1;
      const double CTEQ_DF_t::Q2_max = 1e4;

      CTEQ_DF_t::CTEQ_DF_t( std::string DF_name, std::string option, std::string pT_instructions ) :
         FlavArrayFunc_t( std::string("pT ") + pT_instructions ), cteq(DF_Wrapper::CTEQ_t::Instance())
      {
         // make msg
         constr_msg = DF_name + ": using CTEQ opt. " + option; 

         double w_init = 0;
         if( DF_name == "f1" )
            w_init = 1;

         for( double *p = &weight[0]; p != &weight[11]; ++p )
            (*p) = w_init;

         //cerr << DF_name << ", w init = " << w_init << endl;

         if( DF_name.at(0) == 'g' || DF_name.at(0) == 'h' ){
            if ( option == "B" ){
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = -0.3;
            } else if ( option == "C" ){
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = -0.3;
               weight[ STR_FLAV ] = -0.3;
               weight[ ANTI_STR_FLAV ] = -0.3;
            } else if ( option == "D" ){
               weight[ UP_FLAV ] = -0.5;
               weight[ DOWN_FLAV ] = 0.5;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = 0.5;
               weight[ STR_FLAV ] = 0.5;
               weight[ ANTI_STR_FLAV ] = 0.5;
            } else if ( option == "E" ){
               weight[ UP_FLAV ] = -0.3;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = 0.5;
               weight[ STR_FLAV ] = 0.5;
               weight[ ANTI_STR_FLAV ] = 0.5;
            } else if ( option == "1" ){
               weight[ UP_FLAV ] = 1.0;
               weight[ DOWN_FLAV ] = 1.0;
               weight[ ANTI_UP_FLAV ] = 1.0;
               weight[ ANTI_DOWN_FLAV ] = 1.0;
               weight[ STR_FLAV ] = 1.0;
               weight[ ANTI_STR_FLAV ] = 1.0;
            } else {
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
            };
         };

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &CTEQ_DF_t::Gluon );

      };

      CTEQ_DF_t::~CTEQ_DF_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void CTEQ_DF_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         pT_func->GetRange( var_min.pT, var_max.pT );
      };

      // flags with 1 each relevant variable used
      void CTEQ_DF_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double CTEQ_DF_t::Up( const Var_t& var ) const{
         return weight[ UP_FLAV ] ? weight[ UP_FLAV ]*cteq.Eval( 1, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Down( const Var_t& var ) const{
         return weight[ DOWN_FLAV ] ? weight[ DOWN_FLAV ]*cteq.Eval( 2, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Strange( const Var_t& var ) const{
         return weight[ STR_FLAV ] ? weight[ STR_FLAV ]*cteq.Eval( 3, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Charm( const Var_t& var ) const{
         return weight[ CHM_FLAV ] ? weight[ CHM_FLAV ]*cteq.Eval( 4, var.x, var.Q2 ) : 0;
      };

      double CTEQ_DF_t::Bottom( const Var_t& var ) const{
         return weight[ BOT_FLAV ] ? weight[ BOT_FLAV ]*cteq.Eval( 5, var.x, var.Q2 ) : 0;
      };

      double CTEQ_DF_t::Anti_Up( const Var_t& var ) const{
         return weight[ ANTI_UP_FLAV ] ? weight[ ANTI_UP_FLAV ]*cteq.Eval( -1, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Anti_Down( const Var_t& var ) const{
         return weight[ ANTI_DOWN_FLAV ] ? weight[ ANTI_DOWN_FLAV ]*cteq.Eval( -2, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Anti_Strange( const Var_t& var ) const{
         return weight[ ANTI_STR_FLAV ] ? weight[ ANTI_STR_FLAV ]*cteq.Eval( -3, var.x, var.Q2 ) : 0; 
      };

      double CTEQ_DF_t::Anti_Charm( const Var_t& var ) const{
         return weight[ ANTI_CHM_FLAV ] ? weight[ ANTI_CHM_FLAV ]*cteq.Eval( -4, var.x, var.Q2 ) : 0;
      };

      double CTEQ_DF_t::Anti_Bottom( const Var_t& var ) const{
         return weight[ ANTI_BOT_FLAV ] ? weight[ ANTI_BOT_FLAV ]*cteq.Eval( -5, var.x, var.Q2 ) : 0;
      };

      double CTEQ_DF_t::Gluon( const Var_t& var ) const{
         return weight[ GLUON_FLAV ] ? weight[ GLUON_FLAV ]*cteq.Eval( 0, var.x, var.Q2 ) : 0;
      };

      double CTEQ_DF_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double CTEQ_DF_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };

   };
};
