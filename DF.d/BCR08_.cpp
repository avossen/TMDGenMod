/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/


#include "DF.d/BCR08_.h"

#include "Common.d/Consts.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"


#include <string>
#include <iostream>


namespace TMDGen {
   namespace DF {

      const double BCR08_t::x_min = 1e-5;
      const double BCR08_t::x_max = 1;
      const double BCR08_t::Q2_min = 0.1;
      const double BCR08_t::Q2_max = 1e4;

      BCR08_t::BCR08_t() : FlavArrayFunc_t( "pT NULL" ),
                           m( 0.3 ),
                           C_F_alpha_S( 4./3.*ALPHA_S ),
                           M_s     ( 0.822 ),
                           M_a     ( 1.492),
                           M_aprime( 0.890 ),
                           Lambda_s     ( 0.609 ),
                           Lambda_a     ( 0.716 ),
                           Lambda_aprime( 0.376 ),
                           c_sq_s     ( 0.847*0.847 ),
                           c_sq_a     ( 1.061*1.061 ),
                           c_sq_aprime( 0.880*0.880 )
      {

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &BCR08_t::Gluon );

      };

      BCR08_t::~BCR08_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void BCR08_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
      };

      // flags with 1 each relevant variable with domain cuts
      void BCR08_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
      };

      double BCR08_t::Up( const Var_t& var ) const{
         double temp = c_sq_s * func_s( var ) + c_sq_a * func_a( var );
//          if(!var.integrating && !temp){
//             std::cout << "BCR08 Up: " << c_sq_s << ' ' << func_s( var ) << ' ' << c_sq_a << ' ' << func_a( var ) << std::endl;
//          };

         return c_sq_s * func_s( var ) + c_sq_a * func_a( var );
      };

      double BCR08_t::Down( const Var_t& var ) const{
         return c_sq_aprime * func_aprime( var );
      };

      double BCR08_t::Strange( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Charm( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Bottom( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Anti_Up( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Anti_Down( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Anti_Strange( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Anti_Charm( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Anti_Bottom( const Var_t& var ) const{
         return 0;
      };

      double BCR08_t::Gluon( const Var_t& var ) const{
         return Undefined( var );
      };

      double BCR08_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double BCR08_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };

      double BCR08_t::L_sq_s     ( double x ) const {
         return x*M_s*M_s + (1.-x)*Lambda_s*Lambda_s - x*(1.-x)*PROTON_MASS_SQUARED;
      };

      double BCR08_t::L_sq_a     ( double x ) const {
         return x*M_a*M_a + (1.-x)*Lambda_a*Lambda_a - x*(1.-x)*PROTON_MASS_SQUARED;
      };

      double BCR08_t::L_sq_aprime( double x ) const {
         return x*M_aprime*M_aprime + (1.-x)*Lambda_aprime*Lambda_aprime - x*(1.-x)*PROTON_MASS_SQUARED;
      };



   };
};
