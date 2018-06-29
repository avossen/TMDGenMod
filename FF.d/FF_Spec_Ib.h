/*

--> Gliske's 2010 reworking of Bacchetta, Radici 2006 Spec. Model

    Allows three parameter sets, for various flavor dependence
    Only is for dihadrons decaying into particles of identical mass. (i.e. no K*s)

*/

#ifndef _FF_SPEC_IB_H_
#define _FF_SPEC_IB_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#include "FF.d/Spec_Ia_Params.h"
#include <cmath>

namespace TMDGen {
   namespace FF {

      class FF_Spec_Ib_t : public FlavArrayFunc_t {

      protected:
         typedef DiHad_Spec_Ia::Spec_Ia_Params_t params_t;

         static const double Q2_min, Q2_max, W2_min, W2_max,   // W2, Q2 in units of GeV^2/c^4
            x_min, x_max, y_min, y_max, z_min, z_max;

         params_t *par_0, *par_1, *par_2, *par_array[ GMC_TRANS_N_FLAVORS ];

         // inner functions
         static double tau_s( const params_t& par, const Var_t& ); 
         static double tau_sp_0( const Var_t& ); 
         static double tau_sp_1( const Var_t& ); 
         static double tau_pp_00( const params_t& par, const Var_t& ); 
         static double tau_pp_20( const params_t& par, const Var_t& ); 
         static double tau_pp_21( const params_t& par, const Var_t& ); 
         static double tau_pp_22( const Var_t& ); 
         static double eta( const params_t& par, const Var_t& );

         // put here so inlined
         static double kT_slope( const params_t& par, const Var_t& var ){
            double val = 1;
            if( par.kT_alpha ){
               double temp = par.kT_alpha*pow( var.z, par.kT_beta )*pow( 1.-var.z, par.kT_gamma );

               if( temp ){
                  val = var.kT;
                  val *= var.z/(1.-var.z);
                  val /= temp;
                  val *= val;
                  val = exp( -2.*val );
               };
            };
            return val;
         };

         // vertex functions
         static double F_s_sq( const params_t& par, const Var_t& var );
         double F_p_sq( const params_t& par, const Var_t& var ) const;
         double Re_F_s_star_F_p( const params_t& par, const Var_t& var ) const;
         double Im_F_s_star_F_p( const params_t& par, const Var_t& var ) const;
         static double Re_F_p_rho( const params_t& par, const Var_t& var );
         static double Im_F_p_rho( const params_t& par, const Var_t& var );
         static double Re_F_p_phi( const params_t& par, const Var_t& var );
         static double Im_F_p_phi( const params_t& par, const Var_t& var );
         static double k2( const params_t& par, const Var_t& var );

         static double D1_prefactor( const params_t& par, const Var_t& var );
         static double H1perp_prefactor( const params_t& par, const Var_t& var );

         // par struct allows different settings for different flavors
         double D1_00( const params_t& par, const Var_t& var ) const;
         double D1_11( const params_t& par, const Var_t& var ) const;
         double D1_10( const params_t& par, const Var_t& var ) const;
         double D1_22( const params_t& par, const Var_t& var ) const;
         double D1_21( const params_t& par, const Var_t& var ) const;
         double D1_20( const params_t& par, const Var_t& var ) const;

         double H1perp_11( const params_t& par, const Var_t& var ) const;
         double H1perp_10( const params_t& par, const Var_t& var ) const;
         double H1perp_1m1( const params_t& par, const Var_t& var ) const;

         // pick correct vertex function
         double (*Re_F_p_ptr)( const params_t& par, const Var_t& var );
         double (*Im_F_p_ptr)( const params_t& par, const Var_t& var );

         // pointer to choose correct partial wave
         double (FF_Spec_Ib_t::*FF_lm_ptr)( const params_t& par, const Var_t& var ) const;

         // pick correct prefactor
         double (*prefactor_ptr)( const params_t& par, const Var_t& var );


         // flavor functions
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
         // is_pol == 0 for D1
         // is_pol == 1 for H1perp
         FF_Spec_Ib_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m, bool is_pol ); 
         virtual ~FF_Spec_Ib_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

         // Copy struct into this class
         void CopyParams( int idx, const params_t& params_in );

      };
   };
};


#endif
