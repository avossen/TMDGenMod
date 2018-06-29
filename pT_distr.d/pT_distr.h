

// abstract base class for pT^2 distributions for angular integrated cross sections

#ifndef _PT_DISTR_H_
#define _PT_DISTR_H_

#include "Common.d/Var.h"
#include <string>

namespace TMDGen {
   namespace pT_distr {

      class pT_distr_t {
      protected:
         double min, max;

      public:
         virtual ~pT_distr_t() { /* */ };

         virtual double Eval( const Var_t& var ) const = 0;

         void GetRange( double& min_pT, double& max_pT ) const {
            min_pT = min;
            max_pT = max;
         };

         virtual bool NonNull() const{
            return 1;
         };
      };

      pT_distr_t* Alloc( const std::string& );
   };
};

#endif
