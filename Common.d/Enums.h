
#ifndef _ENUMS_H_
#define _ENUMS_H_

namespace TMDGen {
   namespace Enum {

      // beam and target polarization
      enum pol_state_t { LONG, TRANS, UNPOL };

      // enum final_state_t { SINGLE_HADRON, HADRON_PAIR, VECTOR_MESON };

      // order of cos functions
      enum cos_order_t { P_00, P_10, P_11, P_20, P_21, P_22 };

   };
};

#endif
