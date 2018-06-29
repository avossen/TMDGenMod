/*
  Wrapper for grv98.F
  Has optimization if call GRV multiple times at same point
*/

#ifndef _GRV_WRAPPER_H_
#define _GRV_WRAPPER_H_

#include "Wrapper.h"
#include <string>

namespace TMDGen {
   namespace DF_Wrapper {

      class GRV_t : public Wrapper_t {
         GRV_t();
         virtual double Child_Eval( int iparton, double x, double Q2 );

      public:
         static GRV_t& Instance();

         static int iset;
         static std::string path;
      };

   };
};

#endif
