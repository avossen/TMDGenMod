 

#include "FF.d/D1_Spec_Ia.h"
#include "FF.d/H1perp_Spec_Ia.h"

#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"
#include "Common.d/LundPID.h"

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

      const double H1perp_Spec_Ia_t::Q2_min = 1.0;
      const double H1perp_Spec_Ia_t::Q2_max = 1e6;
      const double H1perp_Spec_Ia_t::x_min = 0.023;
      const double H1perp_Spec_Ia_t::x_max = 0.4;
      const double H1perp_Spec_Ia_t::y_min = 0.1;
      const double H1perp_Spec_Ia_t::y_max = 0.95;
      const double H1perp_Spec_Ia_t::z_min = 0.01;
      const double H1perp_Spec_Ia_t::z_max = 0.99;
      const double H1perp_Spec_Ia_t::W2_min = 0.4;
      const double H1perp_Spec_Ia_t::W2_max = 1e6;

      // constructor
      H1perp_Spec_Ia_t::H1perp_Spec_Ia_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m ) :
         FlavArrayFunc_t( "kT NULL" ),
         params(0)
      {
         if( hadron_1_PID != PI_PLUS && hadron_1_PID != PI_MINUS && hadron_1_PID != PI_ZERO ){
            throw Error::Constructing( "H1perp_Spec_Ia_t", "Only pions are currently available" );
         };

         if( hadron_2_PID != PI_PLUS && hadron_2_PID != PI_MINUS && hadron_2_PID != PI_ZERO ){
            throw Error::Constructing( "H1perp_Spec_Ia_t", "Only pions are currently available" );
         };

         // not yet programmed changing parameter sets

         params = new params_t( DiHad_Spec_Ia::G10 );

         // check l,m state
         if( fabs(m) > l )
            throw Error::Constructing( "H1perp_Spec_Ia_t", "|m| > l" );

         switch ( l ){
         case 0:
            H1perp_lm_ptr = &H1perp_Spec_Ia_t::ALSO_ZERO;
            break;
         case 1:
            if( m == 0 ){
               H1perp_lm_ptr = &H1perp_Spec_Ia_t::H1perp_10;
            } else if ( m== 1) {
               H1perp_lm_ptr = &H1perp_Spec_Ia_t::H1perp_11;
            } else {
               H1perp_lm_ptr = &H1perp_Spec_Ia_t::H1perp_1m1;
            };
            break;
         case 2:
            H1perp_lm_ptr = &H1perp_Spec_Ia_t::ALSO_ZERO;
            break;
         default:
            throw Error::Constructing( "H1perp_Spec_Ia_t", "Invalid `l' state" );
         };

         // set function pointers

         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::Up );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::Up );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &H1perp_Spec_Ia_t::ZERO );

         // make msg
         {
            std::stringstream ss;
            ss << "H_1^perp | " << l << ", " << m << " > using Spec_Ia";
            constr_msg = ss.str();
         }
      };

      H1perp_Spec_Ia_t::~H1perp_Spec_Ia_t(){
         delete params;
      };


      void H1perp_Spec_Ia_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.x = x_min;
         var_max.x = x_max;
         var_min.y = y_min;
         var_max.y = y_max;
         var_min.z = z_min;
         var_max.z = z_max;
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.W2 = W2_min;
         var_max.W2 = W2_max;
      };

      // flags with 1 each variable with domain cuts
      void H1perp_Spec_Ia_t::GetRelevantVar( Var_t& var ) const{
         var.x = 1;
         var.y = 1;
         var.z = 1;
         var.Q2 = 1;
         var.W2 = 1;
      };
      
      double H1perp_Spec_Ia_t::Up( const Var_t& var ) const {
         return params->alpha_H1perp * prefactor(*params, var) * (this->*H1perp_lm_ptr)( *params, var );
      };

      double H1perp_Spec_Ia_t::Down( const Var_t& var ) const {
         return Up( var );
      };

      double H1perp_Spec_Ia_t::Anti_Up( const Var_t& var ) const {
         return Up( var );
      };

      double H1perp_Spec_Ia_t::Anti_Down( const Var_t& var ) const {
         return Up( var );
      };

      double H1perp_Spec_Ia_t::Strange( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Charm( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Bottom( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Anti_Strange( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Anti_Charm( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Anti_Bottom( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Quark( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Gluon( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::Sea( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::ZERO( const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::ALSO_ZERO( const params_t& par, const Var_t& var ) const {
         return 0;
      };

      double H1perp_Spec_Ia_t::H1perp_11( const params_t& par, const Var_t& var ) const {
         double Im_Fs_star_Fp = par.f_s * D1_Spec_Ia_t::Im_F_p( par, var );

         double k2_ = D1_Spec_Ia_t::k2( par, var );

         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double temp1 = -var.R_3mag / var.kT;
         double temp2 = var.kT;
         temp2 *= temp2;
         temp2 += k2_;

         double temp3 = (1.-var.z)*k2_;
         temp3 -= var.z*var.z*var.kT*var.kT;

         return temp1 * temp2 * temp3 * Im_Fs_star_Fp
            * exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) );
      };

      double H1perp_Spec_Ia_t::H1perp_10( const params_t& par, const Var_t& var ) const{
         double Im_Fs_star_Fp = par.f_s * D1_Spec_Ia_t::Im_F_p( par, var );

         double k2_ = D1_Spec_Ia_t::k2( par, var );

         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double temp1 = var.kT;
         temp1 *= temp1;
         temp1 += k2_;
         temp1 *= var.z;
         temp1 += var.had_0.M*var.had_0.M;
         temp1 *= -2;
         temp1 += var.z*k2_;
         temp1 *= var.had_0.M*var.R_3mag / var.z;
         temp1 *= Im_Fs_star_Fp;
         temp1 *= exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) );

         return temp1;
      };

      double H1perp_Spec_Ia_t::H1perp_1m1( const params_t& par, const Var_t& var ) const{
         double Im_Fs_star_Fp = par.f_s * D1_Spec_Ia_t::Im_F_p( par, var );

         double k2_ = D1_Spec_Ia_t::k2( par, var );

         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double temp1 = var.had_0.M;
         temp1 *= -temp1;
         temp1 *= var.R_3mag;
         temp1 *= var.kT;
         temp1 *= Im_Fs_star_Fp;
         temp1 *= exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) );

         return temp1;
      };

      double H1perp_Spec_Ia_t::prefactor( const params_t& par, const Var_t& var ) const{
         double k2_ = D1_Spec_Ia_t::k2( par, var );
         double denom = 8. * PI_SQ * k2_ * k2_;

         return var.R_3mag / denom;
      };

      // Copy over struct
      void H1perp_Spec_Ia_t::CopyParams( const params_t& params_in ){
         (*params) = params_in;
      };


   };
};
