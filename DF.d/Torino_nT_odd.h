/*
   Torino Group's Sivers and Boer Mulers functions

   http://arxiv.org/abs/0807.0166
   http://arxiv.org/abs/0805.2677
   http://arxiv.org/abs/0912.5194

*/

#ifndef _TORINO_NT_ODD_H_
#define _TORINO_NT_ODD_H_

#include "DF.d/Torino_nT_odd_Params.h"

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Flavor.h"
#include "DF_Wrapper.d/GRV.h"
#include "DF_Wrapper.d/GRSV.h"

namespace TMDGen {
   namespace DF {

      class Torino_nT_odd_t : public FlavArrayFunc_t {
         NO_EQ_OP( Torino_nT_odd_t );
         NO_COPY_CONSTR_W_CONST( Torino_nT_odd_t ) : FlavArrayFunc_t(""), f1_GRV(DF_Wrapper::GRV_t::Instance()) { /* */ };

      protected:

         Torino_nT_odd_Params_t params;
         DF_Wrapper::GRV_t &f1_GRV;
         double slope, prefactor;

         double inner_function( flavor_t flav1, GRV_flavor_t flav2, const Var_t& var ) const;

         static const double Q2_min, Q2_max, x_min, x_max, pT_min, pT_max;  // Q2 in units of GeV^2/c^4

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

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
         Torino_nT_odd_t( std::string DF_name, std::string param_code );
         virtual ~Torino_nT_odd_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;
      };

   };
};


#endif
