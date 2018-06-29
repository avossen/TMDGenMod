/*
  Kinematic factor times another fragmentation function
*/

#ifndef _PROP_2HAD_H_
#define _PROP_2HAD_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/Var.h"
#include "FF.d/Full_FF_Set.h"

#include <string>

namespace TMDGen {
   namespace FF {

      class Prop_2had_t : public FlavArrayFunc_t {

      protected:
         std::string prop_to_key;

         const Full_FF_Set_t *FF_set;
         int is_init;
         const double *FF_val_ptr;

         double a, b_z, c_z, b_kT, b_Mh, alpha_Mh, alpha_kT, beta_kT, gamma_kT;

         double func( const Var_t& var ) const;

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

      public:
         // is_pol == 0 for D1
         // is_pol == 1 for H1perp
         Prop_2had_t( std::string this_key, int l, int m, std::string prop_to_key, const Full_FF_Set_t *FF_set_in, int npar, const double *params_in );
         virtual ~Prop_2had_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

         virtual void Prepare_For_Precompute();
      };
   };
};



#endif
