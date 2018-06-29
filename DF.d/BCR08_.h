/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

#ifndef _BCR08_H_
#define _BCR08_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

namespace TMDGen {
   namespace DF {

      class BCR08_t : public FlavArrayFunc_t {
         NO_EQ_OP( BCR08_t );
         NO_COPY_CONSTR_W_CONST( BCR08_t ) : FlavArrayFunc_t("") { /* */ };

      protected:

         static const double Q2_min, Q2_max, x_min, x_max;  // Q2 in units of GeV^2/c^4

         double m, C_F_alpha_S;
         double M_s, M_a, M_aprime;
         double Lambda_s, Lambda_a, Lambda_aprime;
         double c_sq_s, c_sq_a, c_sq_aprime;


         virtual double func_s( const Var_t& var ) const = 0;
         virtual double func_a( const Var_t& var ) const = 0;
         virtual double func_aprime( const Var_t& var ) const = 0;

         double L_sq_s     ( double x ) const;
         double L_sq_a     ( double x ) const;
         double L_sq_aprime( double x ) const;

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
         BCR08_t();
         virtual ~BCR08_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;
      };

   };
};


#endif
