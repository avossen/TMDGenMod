/*
  Simple exponential function, with constant mean and sigma
*/


#ifndef _CONST_EXP_H_
#define _CONST_EXP_H_

#include "pT_distr.d/pT_distr.h"

namespace TMDGen {
   namespace pT_distr {
      class Const_Exp_t : public pT_distr_t {
         bool is_kT;
         double minus_R_sq, normfactor;

      public:
         Const_Exp_t( double alpha, bool is_kT );
         virtual double Eval(const TMDGen::Var_t&) const;

      };
      
   };
};


#endif
