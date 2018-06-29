/*

Generic parent class for fragmentation functions & distribution function that are
actually arrays of functions, one per quark flavor with possibly a gluon type

TODO inline the parent class

*/

#include "Common.d/FlavArrayFunc.h"
#include "Common.d/Var.h"
#include "Common.d/Flavor.h"
#include "Common.d/Exceptions.h"
#include "pT_distr.d/pT_distr.h"

#include <iostream>
#include <string>
using std::cerr;
using std::cout;
using std::endl;


namespace TMDGen {

   const double FF_NOT_SET_CODE = -9999;

   FlavArrayFunc_t::FlavArrayFunc_t( std::string pT_instruction ) : pT_func(0) {

      // first set up pT distribution
      pT_func = pT_distr::Alloc( pT_instruction );
      if( !pT_func )
         throw Error::Constructing( "FlavArrayFunc_t", std::string("error allocating pT (kT) distribution for directive '") + pT_instruction + "'" );

      pT_func_is_nonnull = pT_func->NonNull();

      for( int i=0; i<GMC_TRANS_N_FLAVORS; ++i )
         func_ptr[i] = &FlavArrayFunc_t::Undefined;
   };

   FlavArrayFunc_t::~FlavArrayFunc_t(){
      // always check is non-zero, just to be safe
      if( pT_func )
         delete pT_func;
   };

   double FlavArrayFunc_t::Undefined( const Var_t& ) const{
      throw Error::UndefinedFunction( "FlavArrayFunc_t" );
      return 0;
   };

   double FlavArrayFunc_t::operator() ( const Var_t& var ) const{
      // don't need to check if pT_func pointer is non zero, checked earlier
      // std::cerr << "func_ptr " << this << ' ' << func_ptr << ' ' << var.flavor << std::endl;

     if( var.flavor < UP_FLAV || var.flavor > GLUON_FLAV ){
         std::cerr << var.flavor << std::endl;
         throw Error::SanityCheckFailure("FlavArrayFunc_t::operator(...)", "Invalid flavor" );
     };

//      if( !var.integrating )
//         cerr << "Going to evaluate '" << constr_msg << "' for flavor " << flavor_string[var.flavor] << "\t" << func_ptr[ var.flavor] << endl;

      if( !func_ptr[ var.flavor ] )
         throw Error::SanityCheckFailure("FlavArrayFunc_t::operator(...)",
                                  std::string("Undefined function for flavor '") + flavor_string[ var.flavor ] + "'" );
      //if( !pT_func )
      //   throw Error::SanityCheckFailure("FlavArrayFunc_t::operator(...)", "Undefined pT/kT function." );

//       if( !var.integrating )
//          std::cout << "Eval " << constr_msg << ' ' << (this->*func_ptr[ var.flavor ])( var ) << ' ' << flavor_string[ var.flavor ] << std::endl;

      return pT_func_is_nonnull ? (this->*func_ptr[ var.flavor ])( var ) * pT_func->Eval( var ) : (this->*func_ptr[ var.flavor ])( var );
   };

   const std::string& FlavArrayFunc_t::ConstrMessage() const {
      return constr_msg;
   };

   void FlavArrayFunc_t::Prepare_For_Precompute() {
      // nothing to do
      return;
   };


};
