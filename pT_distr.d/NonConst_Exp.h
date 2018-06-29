/*
  Simple exponential function, with constant mean = x^alpha (1-x)^beta
*/


#ifndef _NONCONST_EXP_H_
#define _NONCONST_EXP_H_

#include "pT_distr.d/pT_distr.h"

namespace TMDGen {
   namespace pT_distr {
      class NonConst_Exp_t : public pT_distr_t {
         bool is_kT;
         double C, alpha, beta;

      public:
         NonConst_Exp_t( int model, double C, bool is_kT_ );
         virtual double Eval(const TMDGen::Var_t&) const;
      };
   };
};


#endif
