/*

--> Gliske's 2010 reworking of Bacchetta, Radici 2006 Spec. Model

    So far only programmed for pions, so need to add
    - terms dependent on the mass difference
    - strange quark flavors

*/

#ifndef _H1PERP_SPEC_IA_H_
#define _H1PERP_SPEC_IA_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#include "FF.d/Spec_Ia_Params.h"

namespace TMDGen {
   namespace FF {

      class H1perp_Spec_Ia_t : public FlavArrayFunc_t {
         NO_EQ_OP( H1perp_Spec_Ia_t );
         NO_COPY_CONSTR_W_CONST( H1perp_Spec_Ia_t ) : FlavArrayFunc_t("") { /* */ };

      protected:
         typedef DiHad_Spec_Ia::Spec_Ia_Params_t params_t;
         static const double Q2_min, Q2_max, W2_min, W2_max,   // W2, Q2 in units of GeV^2/c^4
            x_min, x_max, y_min, y_max, z_min, z_max;

         params_t* params;


         // par struct allows different settings for different flavors
         double H1perp_11( const params_t& par, const Var_t& var ) const;
         double H1perp_10( const params_t& par, const Var_t& var ) const;
         double H1perp_1m1( const params_t& par, const Var_t& var ) const;

         double prefactor( const params_t& par, const Var_t& var ) const;

         // pointer to choose correct partial wave
         double (H1perp_Spec_Ia_t::*H1perp_lm_ptr)( const params_t& par, const Var_t& var ) const;

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
         virtual double ALSO_ZERO( const params_t& par, const Var_t& var ) const;

      public:
         H1perp_Spec_Ia_t( LundPID_t hadron_1_PID, LundPID_t hadron_2_PID, int l, int m );
         virtual ~H1perp_Spec_Ia_t();

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
