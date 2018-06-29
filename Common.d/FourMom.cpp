
#include "Common.d/Exceptions.h"
#include "Common.d/FourMom.h"
#include <cmath>
#include <iostream>
using std::endl;
using std::cout;

namespace TMDGen {

   FourMom_t::FourMom_t() : M(0), E(0), P(0), theta(0), phi(0){ /* */ };

   bool FourMom_t::InBoxAcc(){
      int is_in = (theta > min_theta) && (theta < max_theta);

      if( !theta )
         throw Error::SanityCheckFailure( "FourMom_t::InBoxAcc()", "theta == 0" );

      //      std::cout << theta << ' ' << min_theta << ' ' << max_theta << " | " << is_in << endl;

      if( is_in ){
         double tan_theta = tan( theta );
         double tan_theta_X = fabs(tan( cos(phi)*tan_theta ));
         double tan_theta_Y = fabs(tan( sin(phi)*tan_theta ));

         is_in =
            (tan_theta_X > min_tan_theta_X) && (tan_theta_X < max_tan_theta_X ) &&
            (tan_theta_Y > min_tan_theta_Y) && (tan_theta_Y < max_tan_theta_Y );

//           std::cout << "\t";
//           std::cout << tan_theta_X << ' ' << min_tan_theta_X << ' ' << max_tan_theta_X << " | ";
//           std::cout << tan_theta_Y << ' ' << min_tan_theta_Y << ' ' << max_tan_theta_Y << endl;
      };

      return is_in;
   };

   double FourMom_t::min_tan_theta_X = 0.0;         // = tan( 0.   );
   double FourMom_t::max_tan_theta_X = 0.18196953;  // = tan( 0.18 );
   double FourMom_t::min_tan_theta_Y = 0.030009003; // = tan( 0.03 );
   double FourMom_t::max_tan_theta_Y = 0.15113522;  // = tan( 0.15 );
   double FourMom_t::min_theta = 0.03;
   double FourMom_t::max_theta = 0.30;

};
