/*
  base class to hold a set (actually an associative array) of FlavArrayFunc_t*
*/


#ifndef _FLAVARRAYFUNCSET_H_
#define _FLAVARRAYFUNCSET_H_

#include <map>
#include <string>
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/Var.h"
#include "Common.d/NoCopy.h"

namespace TMDGen {

   class FlavArrayFuncSet_t {
      NO_COPY_CONSTR( FlavArrayFuncSet_t );
      NO_EQ_OP( FlavArrayFuncSet_t );

   protected:
      std::map< std::string, FlavArrayFunc_t* > func_map;

      // free memory
      virtual void Free();

      // inner function for getting intersection of valid domains over whole set
      static void GetVarRange_inner( double& min, double& max, double model_min, double model_max );

      void DisplayConstrMsg() const;

      size_t precomputed_value_size;
      double* precomputed_value;
      void Prepare_for_Precompute();  // to be called from children class constructors

   public:
      FlavArrayFuncSet_t();
      virtual ~FlavArrayFuncSet_t();

      virtual bool Includes( std::string key ) const;

      virtual FlavArrayFunc_t* GetPtr( std::string key ) const;

      virtual int GetVarRange( Var_t& min, Var_t& max ) const;

      void Precompute( const Var_t& var);

      const double* Get_Val_Ptr( std::string key ) const;

   };

};


#endif
