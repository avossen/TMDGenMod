
#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include "DF.d/f1_GRV_LO.h"
#include "Common.d/Var.h"
#include "Common.d/Consts.h"

using namespace TMDGen;

extern "C" {
   void pystpr_();
};

extern struct pystpr_block_t {
   double x, Q2, val;
} pystpr_block_;


int main( int argc, char* argv[] ){

   if( argc != 2 ){
      cerr << "Usage: " << argv[0] << " <max_Q2>" << endl;
      return 127;
   };


   DF::f1_GRV_LO_t f1;

   Var_t min, max, var;

   f1.GetVarRange( min, max );

   max.Q2 = atof( argv[1] );

   var.Q2 = 5.01;
   double E_beam_lab = 27.57;

   double xstep = (max.x - min.x)/100;
   double Q2step = (max.Q2 - min.Q2)/100;

   for( var.x = min.x; var.x < max.x; var.x += xstep ){
      double aveval = 0;
      double aveval2 = 0;

      for( var.Q2 = min.Q2; var.Q2 < max.Q2; var.Q2 += Q2step ){
         double val = f1( DOWN_FLAV, var );

         double val2 = 0;
         pystpr_block_.x = var.x;
         pystpr_block_.Q2 = var.Q2;

         pystpr_();

         val2 = pystpr_block_.val;

         // compute other needed quantities
         var.nu         = var.Q2/var.x/TWICE_PROTON_MASS;
         var.y          = var.nu/E_beam_lab;

         // energy of scattered beam
         double e2_E = E_beam_lab - var.nu;

         double min_physical_Q2 = ELECTRON_MASS_SQUARED * var.y * var.y / (1. - var.y );

         double temp_theta_arg = (var.Q2 - min_physical_Q2) / 4. / E_beam_lab / e2_E;


         if( var.y > 1.0 || e2_E*e2_E < ELECTRON_MASS_SQUARED && (temp_theta_arg < 0 || temp_theta_arg > 1) ){
            val = 0;
            val2 = 0;
         };

         double factor = ( 1 - var.y + var.y*var.y/2 ) / var.x*var.y*var.y;
         val *= factor;
         val2 *= factor;

         aveval += val;
         aveval2 += val2;
      };

      aveval *= Q2step;
      aveval2 *= Q2step;

      std::cout << var.x << ' ' << aveval << ' ' << aveval2 << endl;
   };

   cerr << "Q2 step = " << Q2step << endl;

};
