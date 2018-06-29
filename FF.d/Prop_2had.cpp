/*
  Kinematic factor times another fragmentation function
*/

#include "Prop_2had.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"


#include <string>
#include <cmath>
#include <sstream>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;


namespace TMDGen {
   namespace FF {

      Prop_2had_t::Prop_2had_t( std::string this_key, int l, int m, std::string prop_to_key_in,
                                const Full_FF_Set_t *FF_set_in, int npar, const double *params_in ) :
         FlavArrayFunc_t("kT NULL"), prop_to_key(prop_to_key_in), FF_set(FF_set_in), is_init(0), FF_val_ptr(0) {

         double params[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
         if( npar < 0 || npar > 9 )
            throw Error::Constructing("FF::Prop_2had_t", "Invalid number of parameters: must be in 0 to 9" );

         if( !FF_set )
            throw Error::Constructing("FF::Prop_2had_t", "Given null pointer for FF_Set." );

         // copy over any specified values, and leave the rest zero
         for( int i = 0; i<npar; ++i )
            params[i] = params_in[i];

         a    = params[0];
         b_z  = params[1];
         c_z  = params[2];
         b_kT = params[3];
         b_Mh = params[4];
         alpha_Mh = params[5];
         alpha_kT = params[6];
         beta_kT  = params[7];
         gamma_kT = params[8];

         if( !a )
            throw Error::Constructing("FF::Prop_2had_t", "Must provide non-zero 'a' value" );

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &Prop_2had_t::Up );

         // make msg
         {
            std::stringstream ss;
            ss << this_key << " | " << l << ", " << m << " > prop to " << prop_to_key;
            constr_msg = ss.str();
         };
      };

      Prop_2had_t::~Prop_2had_t() { /* */ };

      void Prop_2had_t::Prepare_For_Precompute(){
         if( !FF_set->Includes( prop_to_key ) )
            throw Error::Base_t("Runtime Error in", "Prop_2had_t",
                                std::string("supposed to be proportional to '") + prop_to_key + "', but key does not exist in FF set" ); 

         // note: this will only work if the other one happens to be evaluated first
         // may need to fix this later
         FF_val_ptr = FF_set->Get_Val_Ptr( prop_to_key );

         is_init = 1;
      };

      double Prop_2had_t::func( const Var_t& var ) const {
         double val = a;

         if( b_z )
            val *= pow( var.z, b_z );
         if( c_z )
            val *= pow( 1.-var.z, c_z );

         if( b_kT )
            val *= pow( var.kT, b_kT );

         if( b_Mh )
            val *= pow( var.had_0.M, b_Mh );

         if( alpha_Mh )
            val *= exp( alpha_Mh * var.had_0.M );

         if( alpha_kT ){
            double lambda = alpha_kT;
            if( beta_kT )
               lambda *= pow( var.z, beta_kT );
            if( gamma_kT )
               lambda *= pow( 1.-var.z, gamma_kT );

            double arg = var.kT/lambda;
            arg *= arg;
            if( alpha_kT < 0 )
               arg = -arg;

            val *= exp( arg );
         };

         return val;
      };

      void Prop_2had_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const {
         // nothing to do
      };

      // flags with 1 each variable with domain cuts
      void Prop_2had_t::GetRelevantVar( Var_t& var ) const {
         // nothing to do
      };

      double Prop_2had_t::Up( const Var_t& var ) const {
//          if( !var.integrating )
//             cout << "**\t" << func(var) << ' ' << *FF_val_ptr << ' ' << var.kT << ' ' << var.z << endl;

         if( !FF_val_ptr )
            throw Error::SanityCheckFailure("Prop_2had_t", "Null pointer to FF value" );

         return ((*FF_val_ptr) == FF_NOT_SET_CODE) ? FF_NOT_SET_CODE : func( var ) * (*FF_val_ptr);
      };

      double Prop_2had_t::Down( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Anti_Up( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Anti_Down( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Strange( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Charm( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Bottom( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Anti_Strange( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Anti_Charm( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Anti_Bottom( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Quark( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Gluon( const Var_t& var ) const {
         return Up( var );
      };

      double Prop_2had_t::Sea( const Var_t& var ) const {
         return Up( var );
      };



   };
};
