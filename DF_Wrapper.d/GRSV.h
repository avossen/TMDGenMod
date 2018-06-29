/*
  Wrapper for grsv2000.F
  Has optimization if call GRSV multiple times at same point
*/

#ifndef _GRSV_WRAPPER_H_
#define _GRSV_WRAPPER_H_

#include "Wrapper.h"
#include <string>

namespace TMDGen {
   namespace DF_Wrapper {

      class GRSV_t : public Wrapper_t {
         GRSV_t();
         virtual double Child_Eval( int iparton, double x, double Q2 );

      public:
         static GRSV_t& Instance();

         static int iset;
         static std::string path;
      };

   };
};

#endif
