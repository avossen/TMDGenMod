/*
  GRSV98 for Delta q, Delta g, standard scenario

  Journal reference: Phys.Rev. D63 (2001) 094005
  DOI: 10.1103/PhysRevD.63.094005
  Report number: DO-TH 2000/14
  Cite as: arXiv:hep-ph/0011215v1
*/

#ifndef _g1_GRSV98_H_
#define _g1_GRSV98_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "pT_distr.d/pT_distr.h"


namespace TMDGen {
   namespace DF {

      class g1_GRSV98_t : public FlavArrayFunc_t {
         friend class h1_Torino_t;

         NO_EQ_OP( g1_GRSV98_t );
         NO_COPY_CONSTR_W_CONST( g1_GRSV98_t ) : FlavArrayFunc_t("") { /* */ };

      protected:

         int nset, member;

         enum param_idx_t { UP_VAL = 0, DOWN_VAL = 1, SEA = 2, GLUON = 3 };

         static const double params_LO[4][3];
         static const double params_NLO[4][3];
         const double (*params_ptr)[3];

         static double prop_function( double x, const double ptr[3] );

         static const double Q2_min, Q2_max, x_min, x_max;  // Q2 in units of GeV^2/c^4

         // flags with 1 each relevant variable used
         virtual void GetRelevantVar( Var_t& var ) const;

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
         g1_GRSV98_t( int nset_in, int member_in, int nlo, std::string pT_instructions );
         virtual ~g1_GRSV98_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;
      };

   };
};


#endif
