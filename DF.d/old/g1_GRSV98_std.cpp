/*
  GRSV98 for Delta q, Delta g, standard scenario

  Journal reference: Phys.Rev. D63 (2001) 094005
  DOI: 10.1103/PhysRevD.63.094005
  Report number: DO-TH 2000/14
  Cite as: arXiv:hep-ph/0011215v1

  DEPRICATED IN FAVOR OF WRAPPING THE FORTRAN CLASSES

*/

#include "DF.d/LHAPDF_fx.hpp"
#include "DF.d/g1_GRSV98_std.h"
#include "pT_distr.d/pT_distr.h"
#include "Common.d/FlavArrayFunc.h"
#include "Common.d/NoCopy.h"
#include "Common.d/Var.h"
#include "Common.d/Exceptions.h"

#include <cmath>
#include <iostream>
#include <string>
#include <sstream>

using std::cerr;
using std::endl;

namespace TMDGen {
   namespace DF {

      const double g1_GRSV98_t::params_LO[4][3] = {
         {  0.851, 0.45, 0.00 },
         { -0.734, 0.49, 0.03 },
         { -0.587, 0.68, 0.00 },
         {  1.669, 1.79, 0.15 }
      };

      const double g1_GRSV98_t::params_NLO[4][3] = {
         {  1.019, 0.52, 0.12 },
         { -0.669, 0.43, 0.00 },
         { -0.272, 0.38, 0.00 },
         {  1.419, 1.43, 0.15 }
      };

      double g1_GRSV98_t::prop_function( double x, const double ptr[3] ){
         return ptr[0]*pow( x, ptr[1] )*pow( 1-x, ptr[2] );
      };

      // will set in constructor
      const double g1_GRSV98_t::x_min = 0;
      const double g1_GRSV98_t::x_max = 0;
      const double g1_GRSV98_t::Q2_min = 0;
      const double g1_GRSV98_t::Q2_max = 0;

      g1_GRSV98_t::g1_GRSV98_t( int nset_in, int member_in, int nlo, std::string pT_instructions ) :
         FlavArrayFunc_t( std::string("pT ") + pT_instructions ), nset(nset_in), member(member_in)
      {
         // make msg
         {
            std::stringstream ss;
            ss << "g_1: using GRSV98, prop. to LHAPDF " << nset << ", " << member << ", nlo = " << nlo << "." << endl;
            constr_msg = ss.str();
         };

         LHAPDF::usePDFMember (nset, member);
         if( member > LHAPDF::numberPDF() || member < 0 )
            throw Error::Constructing( "g1_GRSV98_t", "Requesting invalid member of LHAPDF set");

         *const_cast< double* >( &x_min ) = LHAPDF::getXmin(nset, member);
         *const_cast< double* >( &x_max ) = LHAPDF::getXmax(nset, member);
         *const_cast< double* >( &Q2_min ) = LHAPDF::getQ2min(nset, member);
         *const_cast< double* >( &Q2_max ) = LHAPDF::getQ2max(nset, member);

         if( x_min >= x_max )
            throw Error::Constructing( "g1_GRSV98_t", "x_min >= x_max.  Did you forget to initialize LHA_PDF?");

         if( Q2_min >= Q2_max )
            throw Error::Constructing( "g1_GRSV98_t", "Q2_min >= Q2_max.  Did you forget to initialize LHA_PDF?");

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &g1_GRSV98_t::Gluon );

         params_ptr = ( nlo ? params_NLO : params_LO );
      };

      g1_GRSV98_t::~g1_GRSV98_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void g1_GRSV98_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         pT_func->GetRange( var_min.pT, var_max.pT );
      };

      // flags with 1 each relevant variable used
      void g1_GRSV98_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double g1_GRSV98_t::Up( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return prop_function( var.x, params_ptr[ UP_VAL ] ) * LHAPDF_fx< LHAPDF::UP >( 1, var ) 
            + prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::UBAR >( 0.5, var ) + LHAPDF_fx< LHAPDF::DBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Down( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return prop_function( var.x, params_ptr[ DOWN_VAL ] ) * LHAPDF_fx< LHAPDF::DOWN >( 1, var ) 
            + prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::UBAR >( 0.5, var ) + LHAPDF_fx< LHAPDF::DBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Strange( const Var_t& var ) const{
         return prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::STRANGE >( 0.5, var ) + LHAPDF_fx< LHAPDF::SBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Charm( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV98_t::Bottom( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV98_t::Anti_Up( const Var_t& var ) const{
         return prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::UBAR >( 0.5, var ) + LHAPDF_fx< LHAPDF::DBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Anti_Down( const Var_t& var ) const{
         return prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::UBAR >( 0.5, var ) + LHAPDF_fx< LHAPDF::DBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Anti_Strange( const Var_t& var ) const{
         return prop_function( var.x, params_ptr[ SEA ] ) * ( LHAPDF_fx< LHAPDF::STRANGE >( 0.5, var ) + LHAPDF_fx< LHAPDF::SBAR >( 0.5, var ) );
      };

      double g1_GRSV98_t::Anti_Charm( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV98_t::Anti_Bottom( const Var_t& var ) const{
         return 0;
      };

      double g1_GRSV98_t::Gluon( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return prop_function( var.x, params_ptr[ GLUON ] ) * LHAPDF_fx< LHAPDF::GLUON >( 1, var );
      };

      double g1_GRSV98_t::Quark( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

      double g1_GRSV98_t::Sea( const Var_t& var ) const {
         // this should never be called
         return Undefined( var );
      };

   };
};
