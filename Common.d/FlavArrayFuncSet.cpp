/*
  base class to hold a set (actually an associative array) of FlavArrayFunc_t*
*/

#include "Common.d/Exceptions.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/FlavArrayFuncSet.h"
#include "Common.d/Var.h"

#include <map>
#include <string>
#include <iostream>

#define NTRIALS 3

namespace TMDGen {

   FlavArrayFuncSet_t::FlavArrayFuncSet_t() : precomputed_value_size(0), precomputed_value(0) {
      // not needed, but just to be overly careful
      func_map.clear();
   };

   FlavArrayFuncSet_t::~FlavArrayFuncSet_t(){
      Free();
      if( precomputed_value )
         delete[] precomputed_value;
   };

   bool FlavArrayFuncSet_t::Includes( std::string key ) const {
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.find( key );
      return it != func_map.end();
   };

   FlavArrayFunc_t* FlavArrayFuncSet_t::GetPtr( std::string key ) const {
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.find( key );
      return it == func_map.end() ? 0 : it->second;
   };

   void FlavArrayFuncSet_t::Free(){
      std::map< std::string, FlavArrayFunc_t* >::iterator it;

      for( it = func_map.begin(); it != func_map.end(); ++it ){
         if( it->second ){
            delete it->second;
            it->second = 0;
         };
      };

      func_map.clear();
   };

   // only implemented for x, y, Q2, z, P_hperp, M_hh
   int FlavArrayFuncSet_t::GetVarRange( Var_t& min, Var_t& max ) const{
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it;

      Var_t model_min("FlavArrayFuncSet_t::model_min"),
         model_max("FlavArrayFuncSet_t::model_max"),
         model_used("FlavArrayFuncSet_t::model_used");

      // zero the same set as will be checked
      model_used.Q2 = 0;
      model_min.Q2 = 0;
      model_max.Q2 = 0;
      model_used.x = 0;
      model_min.x = 0;
      model_max.x = 0;
      model_used.y = 0;
      model_min.y = 0;
      model_max.y = 0;
      model_used.z = 0;
      model_min.z = 0;
      model_max.z = 0;
      model_used.pT = 0;
      model_min.pT = 0;
      model_max.pT = 0;
      model_used.kT = 0;
      model_min.kT = 0;
      model_max.kT = 0;
      model_used.had_0.M = 0;
      model_min.had_0.M = 0;
      model_max.had_0.M = 0;

      for( it = func_map.begin(); it != func_map.end(); ++it ){
         it->second->GetVarRange( model_min, model_max );
         it->second->GetRelevantVar( model_used );

         if( model_used.Q2 )
            GetVarRange_inner( min.Q2, max.Q2, model_min.Q2, model_max.Q2 );
         if( model_used.x )
            GetVarRange_inner( min.x, max.x, model_min.x, model_max.x );
         if( model_used.y )
            GetVarRange_inner( min.y, max.y, model_min.y, model_max.y );
         if( model_used.z )
            GetVarRange_inner( min.z, max.z, model_min.z, model_max.z );
         if( model_used.pT )
            GetVarRange_inner( min.pT, max.pT, model_min.pT, model_max.pT );
         if( model_used.kT )
            GetVarRange_inner( min.kT, max.kT, model_min.kT, model_max.kT );
         if( model_used.had_0.M )
            GetVarRange_inner( min.had_0.M, max.had_0.M, model_min.had_0.M, model_max.had_0.M );
      };

      return 0;
   };


   // inner function, to be called for each variable checked
   void FlavArrayFuncSet_t::GetVarRange_inner( double& min, double& max, double model_min, double model_max ){
      if( min == max ){
         // not yet set
         min = model_min;
         max = model_max;
      } else {
         if( min < model_min )
            min = model_min;
         if( max > model_max )
            max = model_max;
      };
   };

   // display constructor messages
   void FlavArrayFuncSet_t::DisplayConstrMsg() const {
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.begin();
      for( ; it != func_map.end(); ++it )
         std::cerr << "\t\t" << it->second->ConstrMessage() << std::endl;
      //std::cerr << "\t\t" << it->first << " : " << it->second->ConstrMessage() << std::endl;
   };

   // to be called from children class constructors
   void FlavArrayFuncSet_t::Prepare_for_Precompute(){
      if( precomputed_value )
         delete precomputed_value;

      precomputed_value_size = func_map.size();
      precomputed_value = new double [precomputed_value_size];
      //std::cout << "precomputed_value at " << precomputed_value << std::endl;

      // prepare the children
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.begin();
      for( ; it != func_map.end(); ++it )
         it->second->Prepare_For_Precompute();

   };

   void FlavArrayFuncSet_t::Precompute( const Var_t& var ){
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.begin();

      if( precomputed_value_size != func_map.size() )
         throw Error::SanityCheckFailure("FlavArrayFuncSet_t::Precompute(...)", "Function map changed size" );

      for( double *p = precomputed_value; p != &precomputed_value[precomputed_value_size]; ++p )
         (*p) = FF_NOT_SET_CODE;

      int all_set = 0;

      for( int itrial = 0; itrial < NTRIALS && !all_set; ++itrial ){
         all_set = 1;

         it = func_map.begin();
         for( int i=0; it != func_map.end(); ++it, ++i ){
            if( precomputed_value[i] == FF_NOT_SET_CODE ){
               precomputed_value[i] = (*it->second)( var );

               if( precomputed_value[i] != precomputed_value[i] )
                  throw Error::SanityCheckFailure( "FlavArrayFuncSet_t::Precompute(...)",
                                                   std::string( "function '") + it->first + "' returned NaN." );

               all_set &= ( precomputed_value[i] != FF_NOT_SET_CODE );
            };
            //            cout << "trial " << itrial << ", i " << i << ", val " << precomputed_value[i] << ", all set = " << all_set << endl;
         };
      };

      if( !all_set )
         throw Error::SanityCheckFailure("FlavArrayFuncSet_t::Precompute(...)", "After max iterations, still not all FFs set." );

   };

   const double* FlavArrayFuncSet_t::Get_Val_Ptr( std::string key ) const{
      std::map< std::string, FlavArrayFunc_t* >::const_iterator it = func_map.find( key );
      if( it == func_map.end() )
         throw Error::SanityCheckFailure( "FlavArrayFuncSet_t::Get_Val_Ptr( std::string key )",
                                          std::string("\n\t\t\tAttempting to access non-mapped key: '") +
                                          key + "'.\n\t\t\tShould have checked Includes(...) first" );

      size_t i = std::distance( func_map.begin(), it);
      if( precomputed_value_size <= i )
         throw Error::SanityCheckFailure( "FlavArrayFuncSet_t::Get_Val_Ptr( std::string key)",
                                          "Map changed since last prepared to precompute" );

      return &precomputed_value[i];
   };


};
