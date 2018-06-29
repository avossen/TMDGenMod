// Template calls to LHAPDF and allow weighting
// and one step memory of called values

#ifndef NO_LHAPDF

#ifndef _LHAPDF_fx_HPP_
#define _LHAPDF_fx_HPP_

#include "LHAPDF/LHAPDF.h"
#include "Common.d/Var.h"
#include <cmath>

namespace TMDGen {
   namespace DF {

      template < LHAPDF::Flavour lha_flav >
      double LHAPDF_fx( const double& weight, const Var_t& var ){
         static double Q2 = var.Q2;
         static double Q = sqrt(Q2);
         static double x = var.x;
         static double val = (weight && x) ? (weight*LHAPDF::xfx(x, Q, lha_flav)/x) : 0;

         if( var.x != x || var.Q2 != Q2 ){
            Q2 = var.Q2;
            Q = sqrt(Q);
            x = var.x;
            val = (weight && x) ? (weight*LHAPDF::xfx(x, Q, lha_flav)/x) : 0;
         };

         return val;
      };

   };
};

#endif
#endif
