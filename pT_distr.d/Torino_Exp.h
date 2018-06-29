/*
  Simple exponential function, with constant mean and sigma
*/


#ifndef _TORINO_EXP_H_
#define _TORINO_EXP_H_

#include "pT_distr.d/pT_distr.h"

namespace TMDGen {
   namespace pT_distr {
      class Torino_Exp_t : public pT_distr_t {
         bool is_kT;
         double ave_pT_sq, normfactor;

      public:
         Torino_Exp_t( double alpha, bool is_kT );
         virtual double Eval(const TMDGen::Var_t&) const;

      };
      
   };
};


#endif
