/*
   Implements BCR06 Spec. Model from
   http://arXiv.org/abs/0807.0323v2
*/

#ifndef _H1LPERP_BCR08_H_
#define _H1LPERP_BCR08_H_

#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#include "DF.d/BCR08_.h"

namespace TMDGen {
   namespace DF {

      class h1Lperp_BCR08_t : public BCR08_t {
         NO_EQ_OP( h1Lperp_BCR08_t );
         NO_COPY_CONSTR( h1Lperp_BCR08_t );

      protected:
         virtual double func_s( const Var_t& var ) const;
         virtual double func_a( const Var_t& var ) const;
         virtual double func_aprime( const Var_t& var ) const;

      public:
         h1Lperp_BCR08_t();
         virtual ~h1Lperp_BCR08_t(){ /* */ };
      };

   };
};


#endif
