 

#include "FF.d/D1_Spec_Ia.h"

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

      const double D1_Spec_Ia_t::Q2_min = 1.0;
      const double D1_Spec_Ia_t::Q2_max = 1e6;
      const double D1_Spec_Ia_t::x_min = 0.023;
      const double D1_Spec_Ia_t::x_max = 0.4;
      const double D1_Spec_Ia_t::y_min = 0.1;
      const double D1_Spec_Ia_t::y_max = 0.95;
      const double D1_Spec_Ia_t::z_min = 0.01;
      const double D1_Spec_Ia_t::z_max = 0.99;
      const double D1_Spec_Ia_t::W2_min = 0.4;
      const double D1_Spec_Ia_t::W2_max = 1e6;

      // constructor
      D1_Spec_Ia_t::D1_Spec_Ia_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m ) :
         FlavArrayFunc_t( "kT NULL" ),
         params(0)
      {
         if( hadron_1_PID != PI_PLUS && hadron_1_PID != PI_MINUS && hadron_1_PID != PI_ZERO ){
            throw Error::Constructing( "D1_Spec_Ia_t", "Only pions are currently available" );
         };

         if( hadron_2_PID != PI_PLUS && hadron_2_PID != PI_MINUS && hadron_2_PID != PI_ZERO ){
            throw Error::Constructing( "D1_Spec_Ia_t", "Only pions are currently available" );
         };

         // not yet programmed changing parameter sets on construction
         // can copy in other sets later though
         params = new params_t( DiHad_Spec_Ia::BR06 );

         // check l,m state
         if( fabs(m) > l )
            throw Error::Constructing( "D1_Spec_Ia_t", "|m| > l" );

         // extra sign to take care of d & ubar having opposite sign as u and dbar
         extra_sign = 1.0;
         if( l == 1 )
            extra_sign = -1.0;

         switch ( l ){
         case 0:
            D1_lm_ptr = &D1_Spec_Ia_t::D1_00;
            break;
         case 1:
            if( m == 0 ){
               D1_lm_ptr = &D1_Spec_Ia_t::D1_10;
            } else {
               D1_lm_ptr = &D1_Spec_Ia_t::D1_11;
            };
            break;
         case 2:
            if( m == 0 ){
               D1_lm_ptr = &D1_Spec_Ia_t::D1_20;
            } else if( m == 1 || m == -1 ) {
               D1_lm_ptr = &D1_Spec_Ia_t::D1_21;
            } else {
               D1_lm_ptr = &D1_Spec_Ia_t::D1_22;
            };
            break;
         default:
            cerr << "\tl = " << l << endl;
            throw Error::Constructing( "D1_Spec_Ia_t", "Invalid `l' state" );
         };


         // set function pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &D1_Spec_Ia_t::Gluon );

         // make msg
         {
            std::stringstream ss;
            ss << "D_1 | " << l << ", " << m << " > using Spec_Ia";
            constr_msg = ss.str();
         }
      };

      D1_Spec_Ia_t::~D1_Spec_Ia_t(){
         delete params;
      };


      void D1_Spec_Ia_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
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
      void D1_Spec_Ia_t::GetRelevantVar( Var_t& var ) const{
         var.x = 1;
         var.y = 1;
         var.z = 1;
         var.Q2 = 1;
         var.W2 = 1;
      };
      
      double D1_Spec_Ia_t::Up( const Var_t& var ) const {
         //std::cout << constr_msg << "\t" << prefactor(*params, var) << ' ' << (this->*D1_lm_ptr)( *params, var ) << std::endl;
         return prefactor(*params, var) * (this->*D1_lm_ptr)( *params, var );
      };

      double D1_Spec_Ia_t::Down( const Var_t& var ) const {
         return extra_sign*Up( var );
      };

      double D1_Spec_Ia_t::Anti_Up( const Var_t& var ) const {
         return extra_sign*Up( var );
      };

      double D1_Spec_Ia_t::Anti_Down( const Var_t& var ) const {
         return Up( var );
      };

      double D1_Spec_Ia_t::Strange( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Charm( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Bottom( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Anti_Strange( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Anti_Charm( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Anti_Bottom( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Quark( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Gluon( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::Sea( const Var_t& var ) const {
         return 0;
      };

      double D1_Spec_Ia_t::ZERO( const Var_t& var ) const {
         return 0;
      };

      // inner functions
      double D1_Spec_Ia_t::tau_s( const params_t& par, const Var_t& var ) const{
         double temp1 = var.z*var.kT;
         temp1 *= temp1;

         double temp2 = par.M_s_factor*var.had_0.M;
         temp2 *= temp2;

         temp1 += temp2;
         temp1 /= (1.-var.z);

         return temp1;
      };

      double D1_Spec_Ia_t::tau_sp_0( const Var_t& var ) const {
         double temp1 = var.z*var.kT;
         temp1 *= temp1;

         double temp2 = var.had_0.M;
         temp2 *= temp2;

         temp1 += temp2;
         temp1 *= 2.0 * var.R_3mag / var.z / var.had_0.M;

         return temp1;
      };

      double D1_Spec_Ia_t::tau_sp_1( const Var_t& var ) const{
         return -4.*var.kT*var.R_3mag;
      };

      double D1_Spec_Ia_t::tau_pp_00( const params_t& par, const Var_t& var ) const {
         double tau_sp_0_ = tau_sp_0( var );

         double tau_sp_1_sq_ = tau_sp_1( var );
         tau_sp_1_sq_ *= tau_sp_1_sq_;

         double eta_ = eta( par, var );

         double val( tau_sp_0_ + eta_ );
         val *= tau_sp_0_;
         val += tau_sp_1_sq_;
         val /= 12.;

         val -= var.R_4mag_sq * tau_s( par, var );

         return val;
      };

      double D1_Spec_Ia_t::tau_pp_20( const params_t& par, const Var_t& var ) const {
         double tau_sp_0_ = tau_sp_0( var );
         double tau_sp_1_sq_ = tau_sp_1( var );
         tau_sp_1_sq_ *= tau_sp_1_sq_;
         double eta_ = eta( par, var );

         double val( tau_sp_0_ + eta_ );
         val *= tau_sp_0_;
         val -= 0.5*tau_sp_1_sq_;
         val /= 6.;

         return val;
      };

      double D1_Spec_Ia_t::tau_pp_21( const params_t& par, const Var_t& var ) const {
         double tau_sp_0_ = tau_sp_0( var );
         double tau_sp_1_ = tau_sp_1( var );
         double eta_ = eta( par, var );

         double val( tau_sp_0_ + eta_ );
         val *= tau_sp_1_;
         val += tau_sp_0_*tau_sp_1_;
         val *= 2./24.;

         return val;
      };

      double D1_Spec_Ia_t::tau_pp_22( const Var_t& var ) const {
         double temp = tau_sp_1( var );
         temp *= temp;
         temp /= 24.;
         return temp;
      };

      double D1_Spec_Ia_t::eta( const params_t& par, const Var_t& var ) const {
         return 2. * var.R_3mag * var.z / var.had_0.M * k2( par, var ); 
      };

      double D1_Spec_Ia_t::D1_00( const params_t& par, const Var_t& var ) const{

         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         double temp = tau_s( par, var ) * F_s_sq( par, var ) * exp( -2.*k2_/lambda_s_sq );
         temp += tau_pp_00( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq );

         temp *= kT_slope( par, var );
         return temp;
      };


      double D1_Spec_Ia_t::D1_11( const params_t& par, const Var_t& var ) const {
         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );
         double M_s = par.M_s_factor*var.had_0.M;

         return tau_sp_1( var ) * M_s * Re_F_s_star_F_p( par, var ) * exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) ) * kT_slope( par, var );
      };

      double D1_Spec_Ia_t::D1_10( const params_t& par, const Var_t& var ) const{
         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );
         double M_s = par.M_s_factor*var.had_0.M;

         return tau_sp_0( var ) * M_s * Re_F_s_star_F_p( par, var ) * exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) ) * kT_slope( par, var );
      };

      double D1_Spec_Ia_t::D1_22( const params_t& par, const Var_t& var ) const{
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_22( var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq ) * kT_slope( par, var );
      };

      double D1_Spec_Ia_t::D1_21( const params_t& par, const Var_t& var ) const{
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_21( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq ) * kT_slope( par, var );
      };

      double D1_Spec_Ia_t::D1_20( const params_t& par, const Var_t& var ) const{
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_20( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq ) * kT_slope( par, var );
      };

      double D1_Spec_Ia_t::F_s_sq( const params_t& par, const Var_t& var ) const {
         return par.f_s*par.f_s;
      };

      double D1_Spec_Ia_t::Re_F_p( const params_t& par, const Var_t& var ) const {
         double rho_temp1 = var.had_0.M*var.had_0.M - RHO_MASS_SQ;
         double omega_temp1 = var.had_0.M*var.had_0.M - OMEGA_MASS_SQ;

         double rho_temp2 = rho_temp1*rho_temp1 + RHO_MASS_WIDTH_SQ;
         double omega_temp2 = omega_temp1*omega_temp1 + OMEGA_MASS_WIDTH_SQ;

         //std::cout << endl << rho_temp2 << ' ' << omega_temp2 << ' ' << var.had_0.M << "\t" << par.f_rho << ' ' << par.f_omega << std::endl;

         return (
                 par.f_rho   * rho_temp1   / rho_temp2 + 
                 par.f_omega * omega_temp1 / omega_temp2  );
      };

      double D1_Spec_Ia_t::Im_F_p( const params_t& par, const Var_t& var ) {
         double rho_temp1 = var.had_0.M*var.had_0.M - RHO_MASS_SQ;
         double omega_temp1 = var.had_0.M*var.had_0.M - OMEGA_MASS_SQ;

         double rho_temp2 = rho_temp1*rho_temp1 + RHO_MASS_WIDTH_SQ;
         double omega_temp2 = omega_temp1*omega_temp1 + OMEGA_MASS_WIDTH_SQ;

         double output = (
                          par.f_rho   * RHO_MASS_WIDTH   / rho_temp2 +
                          par.f_omega * OMEGA_MASS_WIDTH / omega_temp2    );

         if( var.had_0.M < OMEGA_MASS_MINUS_PION_MASS ){
            double m1 = var.had_0.M + PION_MASS;
            double m2 = var.had_0.M - PION_MASS;

            double lambda = OMEGA_MASS_SQ - m1*m1;
            lambda *= OMEGA_MASS_SQ - m2*m2;

            double temp = par.f_omega_prime * sqrt(lambda);
            temp /= FOUR_PI * OMEGA_WIDTH * OMEGA_MASS_SQ * pow( 4 * OMEGA_MASS_SQ * PION_MASS_SQUARED + lambda, 0.25 );

            output += temp;
         };

         return -output;
      };

      double D1_Spec_Ia_t::F_p_sq( const params_t& par, const Var_t& var ) const {
         double r = Re_F_p( par, var );
         double i = Im_F_p( par, var );

         //std::cout << "r = " << r << ", i = " << i << endl;

         return r*r + i*i;
      };

      double D1_Spec_Ia_t::Re_F_s_star_F_p( const params_t& par, const Var_t& var ) const {
         return par.f_s*Re_F_p( par, var );
      };

      double D1_Spec_Ia_t::Im_F_s_star_F_p( const params_t& par, const Var_t& var ) const {
         return par.f_s*Im_F_p( par, var );
      };

      double D1_Spec_Ia_t::k2( const params_t& par, const Var_t& var ) {
         //double M_s = par.M_s_factor*pow(var.had_0.M /PROTON_MASS, par.M_s_const);
         double M_s = par.M_s_factor*var.had_0.M;
         return (var.z*var.kT*var.kT + M_s*M_s)/(1-var.z) + var.had_0.M*var.had_0.M/var.z;
      };

      double D1_Spec_Ia_t::prefactor( const params_t& par, const Var_t& var ) const{
         double k2_ = k2( par, var );
         double denom = 16. * PI_SQ * var.had_0.M * k2_ * k2_ * var.z * var.z;

         return var.R_3mag / denom;
      };

      // Copy over struct
      void D1_Spec_Ia_t::CopyParams( const params_t& params_in ){
         (*params) = params_in;
      };

   };
};
