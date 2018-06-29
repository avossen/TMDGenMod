
// the A,B,C and D functions of y

#include "Common.d/Var.h"

#include <cmath>

namespace TMDGen {
   namespace yFunctions {

      // define the function pointer
      typedef double (*yFunction_ptr)( const Var_t& );

      // the functions
      double A( const Var_t& );
      double B( const Var_t& );
      double C( const Var_t& );
      double V( const Var_t& );
      double W( const Var_t& );

      // order of "y functions", i.e. A(x,y) B(x,y), etc.
      enum y_func_order_t { A_idx, B_idx, C_idx, V_idx, W_idx };

   };
};
