/*
  Includes ctq6pdf.F interface, so can call PartonX6 directly
  Also initializes CTEQ
  Has optimization if call CTEQ multiple times at same point
*/

#ifndef _CTEQ_WRAPPER_H_
#define _CTEQ_WRAPPER_H_

#include <string>

namespace TMDGen {
   namespace CTEQ {


      class Wrapper_t {

         double last_iparton;
         double last_x;
         double last_Q2;
         double last_val;

         Wrapper_t();

      public:
         static Wrapper_t& Instance();
         double Eval( int iparton, double x, double Q2 );

         static int iset;
         static std::string path;
      };

      extern int& iset;
      extern std::string& path; 

   };
};

#endif
