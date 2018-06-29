
#ifndef _FOURMOM_H_
#define _FOURMOM_H_

namespace TMDGen {

   struct FourMom_t {
      double M,   // mass squared
         E,       // energy
         P,       // 3-momentum magnitude
         theta,   // theta
         phi;     // phi

      static double min_tan_theta_X;
      static double max_tan_theta_X;
      static double min_tan_theta_Y;
      static double max_tan_theta_Y;
      static double min_theta;
      static double max_theta;

      FourMom_t();

      bool InBoxAcc();
   };

};

#endif
