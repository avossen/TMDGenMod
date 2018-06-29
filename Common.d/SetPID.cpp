/*
  Interpret string and return the Lund PID code
*/

#include "Common.d/LundPID.h"
//#include "Common.d/Exceptions.h"

#include <iostream>

namespace TMDGen {

   int SetPID( std::string pid_s, LundPID_t& pid ){
      int ierr = 0;

      if( pid_s == "pi+" ){
         pid = PI_PLUS;
      } else if ( pid_s == "pi-" ){
         pid = PI_MINUS;
      } else if ( pid_s == "pi0" ){
         pid = PI_ZERO;
      } else if( pid_s == "K+" ){
         pid = K_PLUS;
      } else if ( pid_s == "K-" ){
         pid = K_MINUS;
      } else {
         std::cerr << "Unknown pid designation '" << pid_s << "'" << std::endl;
         ierr = 1;
      };

      return ierr;
   };

};
