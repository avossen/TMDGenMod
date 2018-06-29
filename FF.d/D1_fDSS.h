/*

--> From the file fDSS.F, which begins with the description

********************************************************************
*                                                                  *
*        fDSS  UNPOLARIZED FRAGMENTATION FUNCTIONS                 *
*  D.de Florian, R.Sassot, M.Stratmann   Phys.Rev.D75 114010 2007  *
*                                 *and*  Phys.Rev.D76 074033 2007  *
*                                                                  *
...
*                                                                  *
*            X                    (between  0.05   and  1.0)       *
*            Q2 = scale in GeV**2 (between  1.0    and  1.D5)      *
*                                                                  *
...
********************************************************************

--> Ported to C++ by S. Gliske, Jan, 2010.

    Actually is product of fDSS with a kT distribution

*/

#ifndef _D1_DFSS_H_
#define _D1_DFSS_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/LundPID.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"

#define D1_DFSS_NQ 24
#define D1_DFSS_NX 35
#define D1_DFSS_NPART 9


#ifndef _FORTRAN_FINT_
#define _FORTRAN_FINT_
extern "C" {
   extern double fint_( const int*, double (*)[2], const int (*)[2], const double (*)[ D1_DFSS_NX + D1_DFSS_NQ ], const double (*)[D1_DFSS_NX*D1_DFSS_NQ] ); 
};
#endif

namespace TMDGen {
   namespace FF {

      class D1_fDSS_t : public FlavArrayFunc_t {
         NO_EQ_OP( D1_fDSS_t );
         NO_COPY_CONSTR_W_CONST( D1_fDSS_t ) : FlavArrayFunc_t("") { /* */ };

      protected:

         static const double Q2_min, Q2_max, z_min, z_max;  // Q2 in units of GeV^2/c^4

         double quark_factor, antiquark_factor;

         enum flavor_identifier_t { u_idx, d_idx, s_idx, c_idx, b_idx, glue_idx, u_val_idx, d_val_idx, s_val_idx };

         // grid for interpolation
         static const double xq_grid[ D1_DFSS_NX + D1_DFSS_NQ ];
         static const double log_xq_grid[ D1_DFSS_NX + D1_DFSS_NQ ];
         static const int N_var;
         static const int N_per[2];

         double grid[ D1_DFSS_NPART ][ D1_DFSS_NX*D1_DFSS_NQ ];

         // inner fuctions for evalutating the grid
         inline double inner_function_1( const Var_t& var, double factor, flavor_identifier_t flav1, flavor_identifier_t flav2 ) const;
         inline double inner_function_2( const Var_t& var, flavor_identifier_t flav1 ) const;

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
         D1_fDSS_t( LundPID_t hadron_PID, int order, const char *path, const char* kT_instruction );
         virtual ~D1_fDSS_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

      };
   };
};


#endif
