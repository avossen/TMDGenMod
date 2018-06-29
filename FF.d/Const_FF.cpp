/*
    Returns a constant times the kT distribution
*/
 

#include "FF.d/Const_FF.h"

#include "Common.d/Exceptions.h"
#include "pT_distr.d/pT_distr.h"

#include <iostream>
#include <string>
#include <sstream>

using std::cerr;
using std::endl;

#include <cmath>

namespace TMDGen {
   namespace FF {

      // constructor
      Const_FF_t::Const_FF_t( std::string name, double val_in, const char* kT_instruction ) : FlavArrayFunc_t( std::string("kT ") + kT_instruction ), val(val_in) {
         // set function pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_FF_t::Gluon );

         // make msg
         {
            std::stringstream ss;
            ss << name << ": const = " << val;
            constr_msg = ss.str();
         }
      };

      Const_FF_t::~Const_FF_t(){
         /* nothing to do */
      };


      void Const_FF_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         pT_func->GetRange( var_min.kT, var_max.kT );
      };

      // flags with 1 each relevant variable used
      void Const_FF_t::GetRelevantVar( Var_t& var ) const{
         var.kT = 1;
      };
      
      double Const_FF_t::Up( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Anti_Up( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Down( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Anti_Down( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Strange( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Anti_Strange( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Charm( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Anti_Charm( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Bottom( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Anti_Bottom( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Quark( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Gluon( const Var_t& var ) const {
         return val;
      };

      double Const_FF_t::Sea( const Var_t& var ) const {
         return val;
      };

   };
};
