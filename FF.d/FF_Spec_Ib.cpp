
//#define DO_OUT_3

#include "FF.d/FF_Spec_Ib.h"

#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"
#include "Common.d/LundPID.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using std::cout;
using std::cerr;
using std::endl;

#include <cmath>
#include <string>
#include <cstdlib>


namespace TMDGen {
   namespace FF {

      const double FF_Spec_Ib_t::Q2_min = 1.0;
      const double FF_Spec_Ib_t::Q2_max = 1e6;
      const double FF_Spec_Ib_t::x_min = 0.023;
      const double FF_Spec_Ib_t::x_max = 0.4;
      const double FF_Spec_Ib_t::y_min = 0.1;
      const double FF_Spec_Ib_t::y_max = 0.95;
      const double FF_Spec_Ib_t::z_min = 0.01;
      const double FF_Spec_Ib_t::z_max = 0.99;
      const double FF_Spec_Ib_t::W2_min = 0.4;
      const double FF_Spec_Ib_t::W2_max = 1e6;

      // constructor
      FF_Spec_Ib_t::FF_Spec_Ib_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m, bool is_pol ) :
         FlavArrayFunc_t( "kT NULL" ),
         par_0(0), par_1(0), par_2(0)
      {
         for( int i=0; i<GMC_TRANS_N_FLAVORS; ++i )
            par_array[i] = 0;

         if( hadron_1_PID == PI_PLUS && hadron_2_PID == PI_MINUS ){
            par_0 = new params_t( "rho0" );

            par_array[ UP_FLAV ] = par_0;
            par_array[ DOWN_FLAV ] = par_0;
            par_array[ ANTI_UP_FLAV ] = par_0;
            par_array[ ANTI_DOWN_FLAV ] = par_0;

            Re_F_p_ptr = &FF_Spec_Ib_t::Re_F_p_rho;
            Im_F_p_ptr = &FF_Spec_Ib_t::Im_F_p_rho;

         } else if ( hadron_1_PID == PI_PLUS && hadron_2_PID == PI_ZERO ){
            par_0 = new params_t( "rho+ u" );
            par_1 = new params_t( "rho+ d" );

            par_array[ UP_FLAV ] = par_0;
            par_array[ DOWN_FLAV ] = par_1;
            par_array[ ANTI_UP_FLAV ] = par_1;
            par_array[ ANTI_DOWN_FLAV ] = par_0;

            Re_F_p_ptr = &FF_Spec_Ib_t::Re_F_p_rho;
            Im_F_p_ptr = &FF_Spec_Ib_t::Im_F_p_rho;

         } else if ( hadron_1_PID == PI_MINUS && hadron_2_PID == PI_ZERO ){
            par_0 = new params_t( "rho- ubar" );
            par_1 = new params_t( "rho- dbar" );

            par_array[ UP_FLAV ] = par_1;
            par_array[ DOWN_FLAV ] = par_0;
            par_array[ ANTI_UP_FLAV ] = par_0;
            par_array[ ANTI_DOWN_FLAV ] = par_1;

            Re_F_p_ptr = &FF_Spec_Ib_t::Re_F_p_rho;
            Im_F_p_ptr = &FF_Spec_Ib_t::Im_F_p_rho;

         } else if ( hadron_1_PID == K_PLUS && hadron_2_PID == K_MINUS ){
            par_0 = new params_t( "phi u" );
            par_1 = new params_t( "phi d" );
            par_2 = new params_t( "phi s" );

            par_array[ UP_FLAV ]   = par_0;
            par_array[ DOWN_FLAV ] = par_1;
            par_array[ STR_FLAV ]  = par_2;
            par_array[ ANTI_UP_FLAV ]   = par_0;
            par_array[ ANTI_DOWN_FLAV ] = par_1;
            par_array[ ANTI_STR_FLAV ]  = par_2;

            Re_F_p_ptr = &FF_Spec_Ib_t::Re_F_p_phi;
            Im_F_p_ptr = &FF_Spec_Ib_t::Im_F_p_phi;
         } else {
            cerr << "FF_Spec_Ib::FF_Spec_Ib(...) Unsuported type: " << hadron_1_PID << " " << hadron_2_PID << endl;
            cerr << "\tAvailable combinations are:" << endl;
            cerr << "\t" << PI_PLUS  << ' ' << PI_ZERO  << " (rho^+)" << endl;
            cerr << "\t" << PI_PLUS  << ' ' << PI_MINUS << " (rho^0)" << endl;
            cerr << "\t" << PI_MINUS << ' ' << PI_ZERO  << " (rho^-)" << endl;
            cerr << "\t" << K_PLUS   << ' ' << K_MINUS  << " (phi)" << endl;
            throw Error::Constructing( "FF_Spec_Ib_t", "Invalid PID" );
         };

         // check l,m state
         if( fabs(m) > l )
            throw Error::Constructing( "FF_Spec_Ib_t", "|m| > l" );

         switch ( l ){
         case 0:
            if( !is_pol )
               FF_lm_ptr = &FF_Spec_Ib_t::D1_00;

            break;
         case 1:
            if( m == 0 ){
               if( is_pol )
                  FF_lm_ptr = &FF_Spec_Ib_t::H1perp_10;
               else
                  FF_lm_ptr = &FF_Spec_Ib_t::D1_10;
            } else {
               if( is_pol ){
                  if ( m==1 )
                     FF_lm_ptr = &FF_Spec_Ib_t::H1perp_11;
                  else
                     FF_lm_ptr = &FF_Spec_Ib_t::H1perp_1m1;
               } else {
                  FF_lm_ptr = &FF_Spec_Ib_t::D1_11;
               };
            };
            break;
         case 2:
            if( !is_pol ){
               if( m == 0 ){
                  FF_lm_ptr = &FF_Spec_Ib_t::D1_20;
               } else if( m == 1 || m == -1 ) {
                  FF_lm_ptr = &FF_Spec_Ib_t::D1_21;
               } else {
                  FF_lm_ptr = &FF_Spec_Ib_t::D1_22;
               };
            };
            break;
         default:
            cerr << "\tl = " << l << endl;
            throw Error::Constructing( "FF_Spec_Ib_t", "Invalid `l' state" );
         };

         if( is_pol )
            prefactor_ptr = &FF_Spec_Ib_t::H1perp_prefactor;
         else
            prefactor_ptr = &FF_Spec_Ib_t::D1_prefactor;


         // set function pointers
         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &FF_Spec_Ib_t::Gluon );

         // make msg
         {
            std::stringstream ss;
            if( is_pol )
               ss << "H_1^perp";
            else
               ss << "D_1";

            ss << " | " << l << ", " << m << " > using Spec_Ib";
            constr_msg = ss.str();
         }
      };

      FF_Spec_Ib_t::~FF_Spec_Ib_t(){
         if( par_0 )
            delete par_0;
         if( par_1 )
            delete par_1;
         if( par_2 )
            delete par_2;
      };


      void FF_Spec_Ib_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
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
      void FF_Spec_Ib_t::GetRelevantVar( Var_t& var ) const{
         var.x = 1;
         var.y = 1;
         var.z = 1;
         var.Q2 = 1;
         var.W2 = 1;
      };
      
      double FF_Spec_Ib_t::Up( const Var_t& var ) const {
         params_t *par = par_array[ UP_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0;
      };

      double FF_Spec_Ib_t::Down( const Var_t& var ) const {
         params_t *par = par_array[ DOWN_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Anti_Up( const Var_t& var ) const {
         params_t *par = par_array[ ANTI_UP_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Anti_Down( const Var_t& var ) const {
         params_t *par = par_array[ ANTI_DOWN_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Strange( const Var_t& var ) const {
         params_t *par = par_array[ STR_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Charm( const Var_t& var ) const {
         params_t *par = par_array[ CHM_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Bottom( const Var_t& var ) const {
         params_t *par = par_array[ BOT_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Anti_Strange( const Var_t& var ) const {
         params_t *par = par_array[ ANTI_STR_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Anti_Charm( const Var_t& var ) const {
         params_t *par = par_array[ ANTI_CHM_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Anti_Bottom( const Var_t& var ) const {
         params_t *par = par_array[ ANTI_BOT_FLAV ];

#ifdef DO_OUT_3
         if( !var.integrating )
            cout << "++\t" << flavor_string[ var.flavor ] << ' ' << (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0) << ' ' << constr_msg << endl;
#endif

         return (par ? (*prefactor_ptr)( *par, var) * (this->*FF_lm_ptr)( *par, var ) : 0);
      };

      double FF_Spec_Ib_t::Quark( const Var_t& var ) const {
         return 0;
      };

      double FF_Spec_Ib_t::Gluon( const Var_t& var ) const {
         return 0;
      };

      double FF_Spec_Ib_t::Sea( const Var_t& var ) const {
         return 0;
      };

      double FF_Spec_Ib_t::ZERO( const Var_t& var ) const {
         return 0;
      };

      // inner functions
      double FF_Spec_Ib_t::tau_s( const params_t& par, const Var_t& var ) {
         double temp1 = var.z*var.kT;
         temp1 *= temp1;

         double temp2 = par.M_s_factor*var.had_0.M;
         temp2 *= temp2;

         temp1 += temp2;
         temp1 /= (1.-var.z);

         return temp1;
      };

      double FF_Spec_Ib_t::tau_sp_0( const Var_t& var ) {
         double temp1 = var.z*var.kT;
         temp1 *= temp1;

         double temp2 = var.had_0.M;
         temp2 *= temp2;

         temp1 += temp2;
         temp1 *= 2.0 * var.R_3mag / var.z / var.had_0.M;

         return temp1;
      };

      double FF_Spec_Ib_t::tau_sp_1( const Var_t& var ) {
         return -4.*var.kT*var.R_3mag;
      };

      double FF_Spec_Ib_t::tau_pp_00( const params_t& par, const Var_t& var ) {
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

      double FF_Spec_Ib_t::tau_pp_20( const params_t& par, const Var_t& var ) {
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

      double FF_Spec_Ib_t::tau_pp_21( const params_t& par, const Var_t& var ) {
         double tau_sp_0_ = tau_sp_0( var );
         double tau_sp_1_ = tau_sp_1( var );
         double eta_ = eta( par, var );

         double val( tau_sp_0_ + eta_ );
         val *= tau_sp_1_;
         val += tau_sp_0_*tau_sp_1_;
         val *= 2./24.;

         return val;
      };

      double FF_Spec_Ib_t::tau_pp_22( const Var_t& var ) {
         double temp = tau_sp_1( var );
         temp *= temp;
         temp /= 24.;
         return temp;
      };

      double FF_Spec_Ib_t::eta( const params_t& par, const Var_t& var ) {
         return 2. * var.R_3mag * var.z / var.had_0.M * k2( par, var ); 
      };

      double FF_Spec_Ib_t::D1_00( const params_t& par, const Var_t& var ) const {

         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         double temp = tau_s( par, var ) * F_s_sq( par, var ) * exp( -2.*k2_/lambda_s_sq );
         temp += tau_pp_00( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq );

         return temp;
      };


      double FF_Spec_Ib_t::D1_11( const params_t& par, const Var_t& var ) const {
         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );
         double M_s = par.M_s_factor*var.had_0.M;

         return tau_sp_1( var ) * M_s * Re_F_s_star_F_p( par, var ) * exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) );
      };

      double FF_Spec_Ib_t::D1_10( const params_t& par, const Var_t& var ) const {
         double lambda_s_sq = par.alpha_s*pow( var.z, par.beta_s)*pow( 1-var.z, par.gamma_s );
         lambda_s_sq *= lambda_s_sq;

         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );
         double M_s = par.M_s_factor*var.had_0.M;

         return tau_sp_0( var ) * M_s * Re_F_s_star_F_p( par, var ) * exp( -k2_*(1./lambda_s_sq + 1./lambda_p_sq) );
      };

      double FF_Spec_Ib_t::D1_22( const params_t& par, const Var_t& var ) const {
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_22( var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq );
      };

      double FF_Spec_Ib_t::D1_21( const params_t& par, const Var_t& var ) const {
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_21( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq );
      };

      double FF_Spec_Ib_t::D1_20( const params_t& par, const Var_t& var ) const {
         double lambda_p_sq = par.alpha_p*pow( var.z, par.beta_p)*pow( 1-var.z, par.gamma_p );
         lambda_p_sq *= lambda_p_sq;

         double k2_ = k2( par, var );

         return tau_pp_20( par, var ) * F_p_sq( par, var ) * exp( -2.*k2_/lambda_p_sq );
      };


      double FF_Spec_Ib_t::H1perp_11( const params_t& par, const Var_t& var ) const {
         double Im_Fs_star_Fp = par.f_s * (*Im_F_p_ptr)( par, var );

         double k2_ = FF_Spec_Ib_t::k2( par, var );

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

      double FF_Spec_Ib_t::H1perp_10( const params_t& par, const Var_t& var ) const {
         double Im_Fs_star_Fp = par.f_s * (*Im_F_p_ptr)( par, var );

         double k2_ = FF_Spec_Ib_t::k2( par, var );

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

      double FF_Spec_Ib_t::H1perp_1m1( const params_t& par, const Var_t& var ) const {
         double Im_Fs_star_Fp = par.f_s * (*Im_F_p_ptr)( par, var );

         double k2_ = FF_Spec_Ib_t::k2( par, var );

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

      double FF_Spec_Ib_t::F_s_sq( const params_t& par, const Var_t& var ) {
         return par.f_s*par.f_s;
      };

      double FF_Spec_Ib_t::Re_F_p_rho( const params_t& par, const Var_t& var ) {
         double rho_temp1 = var.had_0.M*var.had_0.M - RHO_MASS_SQ;
         double omega_temp1 = var.had_0.M*var.had_0.M - OMEGA_MASS_SQ;

         double rho_temp2 = rho_temp1*rho_temp1 + RHO_MASS_WIDTH_SQ;
         double omega_temp2 = omega_temp1*omega_temp1 + OMEGA_MASS_WIDTH_SQ;

         //std::cout << endl << rho_temp2 << ' ' << omega_temp2 << ' ' << var.had_0.M << "\t" << par.f_rho << ' ' << par.f_omega << std::endl;

         return (
                 par.f_rho   * rho_temp1   / rho_temp2 + 
                 par.f_omega * omega_temp1 / omega_temp2  );
      };

      double FF_Spec_Ib_t::Im_F_p_rho( const params_t& par, const Var_t& var ) {
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


      double FF_Spec_Ib_t::Re_F_p_phi( const params_t& par, const Var_t& var ) {
         double temp1 = var.had_0.M*var.had_0.M - PHI_MASS_SQ;
         double temp2 = temp1*temp1 + PHI_MASS_WIDTH_SQ;
         //std::cout << endl << rho_temp2 << ' ' << omega_temp2 << ' ' << var.had_0.M << "\t" << par.f_rho << ' ' << par.f_omega << std::endl;

         // for the case of phi, f_rho is used for the f_phi coupling
         // as never need both
         return ( par.f_rho  * temp1 / temp2 );
      };

      double FF_Spec_Ib_t::Im_F_p_phi( const params_t& par, const Var_t& var ) {
         double temp1 = var.had_0.M*var.had_0.M - PHI_MASS_SQ;
         double temp2 = temp1*temp1 + PHI_MASS_WIDTH_SQ;

         // for the case of phi, f_rho is used for the f_phi coupling
         // as never need both
         return -par.f_rho*PHI_MASS_WIDTH/temp2;
      };

      double FF_Spec_Ib_t::F_p_sq( const params_t& par, const Var_t& var ) const {
         double r = (*Re_F_p_ptr)( par, var );
         double i = (*Im_F_p_ptr)( par, var );

         return r*r + i*i;
      };

      double FF_Spec_Ib_t::Re_F_s_star_F_p( const params_t& par, const Var_t& var ) const {
         return par.f_s*(*Re_F_p_ptr)( par, var );
      };

      double FF_Spec_Ib_t::Im_F_s_star_F_p( const params_t& par, const Var_t& var ) const {
         return par.f_s*(*Im_F_p_ptr)( par, var );
      };

      double FF_Spec_Ib_t::k2( const params_t& par, const Var_t& var ) {
         double M_s = par.M_s_factor*var.had_0.M;
         return (var.z*var.kT*var.kT + M_s*M_s)/(1-var.z) + var.had_0.M*var.had_0.M/var.z;
      };

      double FF_Spec_Ib_t::D1_prefactor( const params_t& par, const Var_t& var ) {
         double k2_ = k2( par, var );
         double denom = 16. * PI_SQ * var.had_0.M * k2_ * k2_ * var.z * var.z;

         return var.R_3mag / denom * kT_slope( par, var );
      };

      double FF_Spec_Ib_t::H1perp_prefactor( const params_t& par, const Var_t& var ) {
         double k2_ = FF_Spec_Ib_t::k2( par, var );
         double denom = 8. * PI_SQ * k2_ * k2_;

         return par.alpha_H1perp * var.R_3mag / denom * kT_slope( par, var );
      };

      // Copy over struct
      void FF_Spec_Ib_t::CopyParams( int idx, const params_t& params_in ){
         if( idx == 0 ){
            (*par_0) = params_in;
         } else if ( idx == 1 ){
            (*par_1) = params_in;
         } else if ( idx == 2 ){
            (*par_2) = params_in;
         };
      };


   };
};

