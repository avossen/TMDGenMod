// Abstract base class for QED Radiative Corrections

#ifndef RADCOR_H_
#define RADCOR_H_

#include "Common.d/ParseInput.h"

namespace TMDGen {
   namespace RadCor {

      class RadCor_t {

      public:
         RadCor_t();
         virtual int Initialize(  const sgUtil::ParseInputReturn_t &parsed_input ) = 0;
         virtual int Cor_InitialState( int& did_radiate, const FourMom_t& lep_0, FourMom_t& lep_1 ) = 0;
         virtual int Cor_FinalState( int& did_radiate, const FourMom_t& lep_2, FourMom_t& lep_3 ) = 0;
         virtual int FinalCorrection( Var_t& var ) = 0;
      };

   };
};

#endif
