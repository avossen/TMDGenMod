/*
  Wrapper for LHAPDF 6 times a p_T distribution
*/

#ifndef NO_LHAPDF

#ifndef _LHAPDF_DF_H_
#define _LHAPDF_DF_H_

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "pT_distr.d/pT_distr.h"


namespace TMDGen {
   namespace DF {

      class LHAPDF_DF_t : public FlavArrayFunc_t {
         NO_EQ_OP( LHAPDF_DF_t );
         NO_COPY_CONSTR_W_CONST( LHAPDF_DF_t ) : FlavArrayFunc_t("") { /* */ };

      protected:

         int nset, member;
         double weight[11];

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
         LHAPDF_DF_t( std::string DF_name, std::string weight_option, int nset_in, int member_in,
                      std::string pT_instructions );
         virtual ~LHAPDF_DF_t();

         // To query allowed values for variables
         virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const;
      };

   };
};


#endif

#endif
