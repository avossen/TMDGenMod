/*
  Constants needed for various things
*/

#ifndef _CONSTS_H_
#define _CONSTS_H_

#ifndef PI
#define PI 3.14159265358979323846
#endif

#ifndef TWO_PI
#define TWO_PI 6.283185307179586232
#endif

#ifndef FOUR_PI
#define FOUR_PI 12.5663706143592
#endif

#ifndef PI_SQ
#define PI_SQ 9.86960440108936
#endif

// in GeV^2
#define PROTON_MASS           0.938272029
#define TWICE_PROTON_MASS     1.876544058
#define PROTON_MASS_SQUARED   0.880354400
#define PION_MASS             0.13957018
#define PION_MASS_SQUARED     0.0194798351
#define KAON_MASS_SQUARED     0.243717
#define ELECTRON_MASS         0.000510998903
#define ELECTRON_MASS_SQUARED 0.000000261119894197170688

#define RHO_MASS      0.776
#define RHO_MASS_SQ   0.602176
#define RHO_WIDTH     0.150
#define OMEGA_MASS    0.783
#define OMEGA_MASS_SQ 0.613089
#define OMEGA_WIDTH   0.008
#define PHI_MASS      1.019455
#define PHI_MASS_SQ   1.0392885
#define PHI_WIDTH     0.00426


#define OMEGA_MASS_MINUS_PION_MASS 0.643
#define RHO_MASS_WIDTH             0.1164
#define OMEGA_MASS_WIDTH           0.006264
#define RHO_MASS_WIDTH_SQ          0.01354896
#define OMEGA_MASS_WIDTH_SQ        3.923769e-5

#define PHI_MASS_WIDTH             0.004343
#define PHI_MASS_WIDTH_SQ          1.886e-05

// EM coupling strength
#define ALPHA_EM              0.00729735257

// Strong coupling constant at 2 GeV
// from http://www-theory.lbl.gov/~ianh/alpha/alpha.html
#define ALPHA_S 0.29942


// classical electron radius, in microbarns
#define R0 2.818e19

// in microbarns GeV^2
#define HBARC2 389.379

#endif
