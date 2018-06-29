/*
    Returns a constant times the kT distribution
*/
 

#include "DF.d/Const_DF.h"

#include "Common.d/Exceptions.h"
#include "pT_distr.d/pT_distr.h"

#include <iostream>
#include <string>
#include <sstream>
using std::cerr;
using std::endl;

#include <cmath>

namespace TMDGen {
   namespace DF {

      // constructor
      Const_DF_t::Const_DF_t( double val_in, const char* pT_instruction ) : FlavArrayFunc_t( std::string("pT ") + pT_instruction ), val(val_in) {
         // set function pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Const_DF_t::Gluon );

         // make msg
         {
            std::stringstream ss;
            ss << "f1: set to " << val << " x y^2";
            constr_msg = ss.str();
         }
      };

      Const_DF_t::~Const_DF_t(){
         /* nothing to do */
      };


      void Const_DF_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.y = 1e-5;
         var_max.y = 1.0;
         var_min.x = 1e-5;
         var_max.x = 1.0;
         pT_func->GetRange( var_min.pT, var_max.pT );
      };

      // flags with 1 each relevant variable used
      void Const_DF_t::GetRelevantVar( Var_t& var ) const{
         var.pT = 1;
      };
      
      double Const_DF_t::Up( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Anti_Up( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Down( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Anti_Down( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Strange( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Anti_Strange( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Charm( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Anti_Charm( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Bottom( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Anti_Bottom( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Quark( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Gluon( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

      double Const_DF_t::Sea( const Var_t& var ) const {
         return val; //*var.x*var.y*var.y;
      };

   };
};
