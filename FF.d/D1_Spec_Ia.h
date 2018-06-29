/*

--> Gliske's 2010 reworking of Bacchetta, Radici 2006 Spec. Model

    So far only programmed for pions, so need to add
    - terms dependent on the mass difference
    - strange quark flavors

*/

#ifndef _D1_SPEC_IA_H_
#define _D1_SPEC_IA_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#include "FF.d/Spec_Ia_Params.h"

#include <cmath>

namespace TMDGen {
   namespace FF {

      class D1_Spec_Ia_t : public FlavArrayFunc_t {
         friend class H1perp_Spec_Ia_t;

         NO_EQ_OP( D1_Spec_Ia_t );
         NO_COPY_CONSTR_W_CONST( D1_Spec_Ia_t ) : FlavArrayFunc_t("") { /* */ };

      protected:
         typedef DiHad_Spec_Ia::Spec_Ia_Params_t params_t;
         static const double Q2_min, Q2_max, W2_min, W2_max,   // W2, Q2 in units of GeV^2/c^4
            x_min, x_max, y_min, y_max, z_min, z_max;

         double extra_sign;

         params_t* params;

         // inner functions
         double tau_s( const params_t& par, const Var_t& ) const; 
         double tau_sp_0( const Var_t& ) const; 
         double tau_sp_1( const Var_t& ) const; 
         double tau_pp_00( const params_t& par, const Var_t& ) const; 
         double tau_pp_20( const params_t& par, const Var_t& ) const; 
         double tau_pp_21( const params_t& par, const Var_t& ) const; 
         double tau_pp_22( const Var_t& ) const; 
         double eta( const params_t& par, const Var_t& ) const;

         // put here so inlined
         static double kT_slope( const params_t& par, const Var_t& var ){
            double val = par.kT_alpha*var.kT*pow( var.z, par.kT_beta )*pow( 1.-var.z, par.kT_gamma );
            return exp( -val*val );
         };

         // par struct allows different settings for different flavors
         double D1_00( const params_t& par, const Var_t& var ) const;
         double D1_11( const params_t& par, const Var_t& var ) const;
         double D1_10( const params_t& par, const Var_t& var ) const;
         double D1_22( const params_t& par, const Var_t& var ) const;
         double D1_21( const params_t& par, const Var_t& var ) const;
         double D1_20( const params_t& par, const Var_t& var ) const;

         // vertex functions
         double F_s_sq( const params_t& par, const Var_t& var ) const;
         double F_p_sq( const params_t& par, const Var_t& var ) const;
         double Re_F_s_star_F_p( const params_t& par, const Var_t& var ) const;
         double Im_F_s_star_F_p( const params_t& par, const Var_t& var ) const;
         double Re_F_p( const params_t& par, const Var_t& var ) const;
         static double Im_F_p( const params_t& par, const Var_t& var );
         static double k2( const params_t& par, const Var_t& var );
         double prefactor( const params_t& par, const Var_t& var ) const;

         // pointer to choose correct partial wave
         double (D1_Spec_Ia_t::*D1_lm_ptr)( const params_t& par, const Var_t& var ) const;

         virtual double Up( const Var_t& var ) const;
         virtual double Down( const Var_t& var ) const;
         virtual double Strange( const Var_t& var ) const;
         virtual double Charm( const Var_t& var ) const;
         virtual double Bottom( const Var_t& var ) const;

         virtual double Quark( const Var_t& var ) const;
         virtual double Gluon( const Var_t& var ) const;
         virtual double Sea( const Var_t& var ) const;

         virtual double Anti_Up( const Var_t& var ) const;
         virtual double Anti_Down( const Var_t& var ) const;
         virtual double Anti_Strange( const Var_t& var ) const;
         virtual double Anti_Charm( const Var_t& var ) const;
         virtual double Anti_Bottom( const Var_t& var ) const;

         virtual double ZERO( const Var_t& var ) const;

      public:
         D1_Spec_Ia_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m );
         virtual ~D1_Spec_Ia_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;


         // Copy struct into this class
         void CopyParams( const params_t& params_in );

      };
   };
};


#endif
