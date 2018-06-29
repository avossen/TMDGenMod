/*
  Simple exponential function, with constant mean and sigma
*/


#ifndef _NULL_pT_H_
#define _NULL_pT_H_

#include "pT_distr.d/pT_distr.h"
#include "Common.d/Exceptions.h"

namespace TMDGen {
   namespace pT_distr {
      class Null_pT_t : public pT_distr_t {
      public:
         Null_pT_t() {
            min = 0;
            max = 1e300;
         };

         virtual double Eval(const TMDGen::Var_t&) const {
            throw Error::SanityCheckFailure("pT_distr::Null_pT", "Eval function called" );
            return 0;
         };

         virtual bool NonNull() const{
            return 0;
         };
      };
      
   };
};


#endif
