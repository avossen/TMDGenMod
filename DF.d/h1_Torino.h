/*
   Torino Group's transversity function
   arXiv:hep-ph/0701006v3
   arXiv:0812.4366v1
*/

#ifndef _h1_Torino_H_
#define _h1_Torino_H_

#include "DF.d/h1_Torino_Params.h"

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Flavor.h"
#include "DF_Wrapper.d/GRV.h"
#include "DF_Wrapper.d/GRSV.h"

namespace TMDGen {
   namespace DF {

      class h1_Torino_t : public FlavArrayFunc_t {
         NO_EQ_OP( h1_Torino_t );
         NO_COPY_CONSTR_W_CONST( h1_Torino_t ) :
            FlavArrayFunc_t(""),
            f1_GRV(DF_Wrapper::GRV_t::Instance()),
            g1_GRSV(DF_Wrapper::GRSV_t::Instance())
            { /* */ };

      protected:

         h1_Torino_Params_t params;
         DF_Wrapper::GRV_t &f1_GRV;
         DF_Wrapper::GRSV_t &g1_GRSV;

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
         h1_Torino_t( std::string param_code );
         virtual ~h1_Torino_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;
      };

   };
};


#endif
