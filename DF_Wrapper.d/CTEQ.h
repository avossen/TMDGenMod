/*
  Includes ctq6pdf.F interface, so can call PartonX6 directly
  Also initializes CTEQ
  Has optimization if call CTEQ multiple times at same point
*/

#ifndef _CTEQ_WRAPPER_H_
#define _CTEQ_WRAPPER_H_

#include "Wrapper.h"
#include <string>

namespace TMDGen {
   namespace DF_Wrapper {

      class CTEQ_t : public Wrapper_t {
         CTEQ_t();
         virtual double Child_Eval( int iparton, double x, double Q2 );

      public:
         static CTEQ_t& Instance();

         static int iset;
         static std::string path;
      };

   };
};

#endif
