/*

--> Ported from the file pkhff.F, which begins with description

***********************************************************************
*                                                                     *
*  LO and NLO FRAGMENTATION FUNCTIONS                                 *
*  for charged pions, kaons and the inclusive sum of charged hadrons  *
*                                                                     *
*  as in S. Kretzer (kretzer@pa.msu.edu):                             *
*                                                                     * 
*  `Fragmentaion Functions from Flavour-inclusive and Flavour-tagged  *
*  e^+ e^- Annihilations';                                            *
*  Phys. Rev. D 62, 054001 (2000)                                     *
*                                                                     *
*  See above reference for details!                                   *  
*								      *	
...
*                                                                     *
*  This interpolation routine returs FFs in the range:                *
*             Z                    (between  0.01   and  1.0)         *
*             Q2 = scale in GeV^2  (between  0.8    and  1.D6)        *
*             Attention: Z \lesssim 0.05   is (very!) delicate        *
*                        ( as discussed in above reference )          *
*                                                                     *
*  The quality of the interpolation (as compared to the exact         * 
*  evolution) is about:                                               *
*  | z < 0.75  | z < 0.9 | z < 1.0 |                                  *                       
*   -------------------------------                                   *
*  |  < 3%     |  < 10%  |  > 10%  |                                  *
*                                                                     *
...
*                                                                     *
***********************************************************************

--> Ported to C++ by S. Gliske, Jan, 2010.

    Actually outputs D1( z, Q2 ) times a kT distribution

*/

#ifndef _D1_KRETZER_H_
#define _D1_KRETZER_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#define D1_KRETZER_NQ 27
#define D1_KRETZER_NZ 49
#define D1_KRETZER_NPART 9

#ifndef _FORTRAN_FINT_
#define _FORTRAN_FINT_
extern "C" {
   extern double fint_( const int*, double (*)[2], const int (*)[2], const double (*)[ D1_KRETZER_NZ + D1_KRETZER_NQ ], const double (*)[D1_KRETZER_NZ*D1_KRETZER_NQ] ); 
};
#endif

namespace TMDGen {
   namespace FF {

      class D1_Kretzer_t : public FlavArrayFunc_t {
         NO_EQ_OP( D1_Kretzer_t );
         NO_COPY_CONSTR_W_CONST( D1_Kretzer_t ) : FlavArrayFunc_t("") { /* */ };

      protected:

         static const double Q2_min, Q2_max, z_min, z_max;  // Q2 in units of GeV^2/c^4

         int shift, do_ave;

         enum flavor_identifier_t { u_idx, u_bar_idx, d_idx, d_bar_idx, s_idx, s_bar_idx, c_idx, b_idx, glue_idx };

         // grid for interpolation
         static const double zq_grid[ D1_KRETZER_NZ + D1_KRETZER_NQ ];
         static const int N_var;
         static const int N_per[2];

         double grid[ D1_KRETZER_NPART ][ D1_KRETZER_NZ * D1_KRETZER_NQ ];

         // inner fuctions for evalutating the grid
         inline double inner_function_1( const int& iflav, const Var_t& var ) const;
         inline double inner_function_2( const flavor_identifier_t& flav, const Var_t& var, int a ) const;

         virtual double Up( const Var_t& var ) const;
         virtual double Down( const Var_t& var ) const;
         virtual double Strange( const Var_t& var ) const;
         virtual double Charm( const Var_t& var ) const;
         virtual double Bottom( const Var_t& var ) const;

         virtual double Quark( const Var_t& var ) const;
         virtual double Gluon( const Var_t& var ) const;
         virtual double Sea( const Var_t& var ) const;

         virtual double Anti_Up( const Var_t& var ) const;
         virtual double Anti_Down( const Var_t& var ) const;
         virtual double Anti_Strange( const Var_t& var ) const;
         virtual double Anti_Charm( const Var_t& var ) const;
         virtual double Anti_Bottom( const Var_t& var ) const;

      public:
         D1_Kretzer_t( LundPID_t hadron_PID, int order, const char *path, const char* kT_instruction );
         virtual ~D1_Kretzer_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

      };
   };
};


#endif
