/*
    Returns a constant times the kT distribution
*/

#ifndef _CONST_FF_H_
#define _CONST_FF_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include <string>

namespace TMDGen {
   namespace FF {

      class Const_FF_t : public FlavArrayFunc_t {
         NO_EQ_OP( Const_FF_t );
         NO_COPY_CONSTR_W_CONST( Const_FF_t ) : FlavArrayFunc_t("") { /* */ };

      protected:
         double val;

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
         Const_FF_t( std::string name, double val, const char* kT_instruction );
         virtual ~Const_FF_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

      };
   };
};


#endif
