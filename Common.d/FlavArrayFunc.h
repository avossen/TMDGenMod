/*

Generic parent class for fragmentation functions & distribution function that are
actually arrays of functions, one per quark flavor with possibly a gluon type

*/

#ifndef _FLAVARRAYFUNC_H_
#define _FLAVARRAYFUNC_H_

#include "Common.d/Var.h"
#include "Common.d/Flavor.h"
#include "Common.d/NoCopy.h"
#include "pT_distr.d/pT_distr.h"

#include <string>

namespace TMDGen {

   class FlavArrayFunc_t;

   extern const double FF_NOT_SET_CODE;

   typedef double (FlavArrayFunc_t::*FlavFunc_ptr_t)( const Var_t& var ) const;

   // alias this base type as the base type for
   // distribution and fragmentation functions
   typedef FlavArrayFunc_t DF_t;
   typedef FlavArrayFunc_t FF_t;

   class FlavArrayFunc_t {
      NO_COPY_CONSTR_W_CONST( FlavArrayFunc_t ) : pT_func(0) { /* */ };
      NO_EQ_OP( FlavArrayFunc_t );

   protected:
      // pT (or kT for FFs) distribution
      const pT_distr::pT_distr_t* pT_func;
      bool pT_func_is_nonnull;

      // array per flavor
      FlavFunc_ptr_t func_ptr[ GMC_TRANS_N_FLAVORS ];

      // message about the model constructed
      std::string constr_msg;

      // undefined function--throws an exception
      double Undefined( const Var_t& var ) const;

      // possible flavor functions in the array
      virtual double Up( const Var_t& var ) const = 0;
      virtual double Down( const Var_t& var ) const = 0;
      virtual double Strange( const Var_t& var ) const = 0;
      virtual double Charm( const Var_t& var ) const = 0;
      virtual double Bottom( const Var_t& var ) const = 0;

      virtual double Quark( const Var_t& var ) const = 0;
      virtual double Gluon( const Var_t& var ) const = 0;
      virtual double Sea( const Var_t& var ) const = 0;

      virtual double Anti_Up( const Var_t& var ) const = 0;
      virtual double Anti_Down( const Var_t& var ) const = 0;
      virtual double Anti_Strange( const Var_t& var ) const = 0;
      virtual double Anti_Charm( const Var_t& var ) const = 0;
      virtual double Anti_Bottom( const Var_t& var ) const = 0;

   public:
      FlavArrayFunc_t( std::string pT_instruction );
      virtual ~FlavArrayFunc_t();

      // evaluate correction function for given flavor
      double operator() ( const Var_t& var ) const;

      // To query allowed values for variables
      virtual void GetVarRange( Var_t& var_min, Var_t& var_max ) const = 0;

      // flags with 1 each relevant variable used
      // does not flag with 0 the variables not used
      // so that the Var_t can be extended without changing
      // all the children of this class
      // i.e. WARNING: makes sure var is initialized to zero
      // before calling this function
      virtual void GetRelevantVar( Var_t& var ) const = 0;

      // in case children need access to FlavArrayFuncSet's 
      // array of precomputed values
      virtual void Prepare_For_Precompute();

      const std::string& ConstrMessage() const;
   };
};

#endif
