/*
  DF from GRV98 wrapper times a p_T distribution
*/

#ifndef _f1_GRV98_H_
#define _f1_GRV98_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "DF_Wrapper.d/GRV.h"

namespace TMDGen {
   namespace DF {

      class f1_GRV98_t : public FlavArrayFunc_t {
      private:
         NO_EQ_OP( f1_GRV98_t );
         NO_COPY_CONSTR_W_CONST( f1_GRV98_t ) : FlavArrayFunc_t(""), GRV(DF_Wrapper::GRV_t::Instance()) { /* */ };

      protected:
         DF_Wrapper::GRV_t& GRV;

         static const double Q2_min, Q2_max, x_min, x_max;  // Q2 in units of GeV^2/c^4

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
         f1_GRV98_t( std::string pT_instructions );
         virtual ~f1_GRV98_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

      };
   };
};


#endif
