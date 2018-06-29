/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

#ifndef _G1L_BCR08_H_
#define _G1L_BCR08_H_

#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#include "DF.d/BCR08_.h"

namespace TMDGen {
   namespace DF {

      class g1L_BCR08_t : public BCR08_t {
         NO_EQ_OP( g1L_BCR08_t );
         NO_COPY_CONSTR( g1L_BCR08_t );

      protected:
         virtual double func_s( const Var_t& var ) const;
         virtual double func_a( const Var_t& var ) const;
         virtual double func_aprime( const Var_t& var ) const;

      public:
         g1L_BCR08_t();
         virtual ~g1L_BCR08_t(){ /* */ };
      };

   };
};


#endif
