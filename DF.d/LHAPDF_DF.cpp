/*
  Wrapper for LHAPDF times a p_T distribution
*/

#ifndef NO_LHAPDF

#include "DF.d/LHAPDF_fx.hpp"
#include "DF.d/LHAPDF_DF.h"
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

      // will set in constructor
      const double LHAPDF_DF_t::x_min = 0;
      const double LHAPDF_DF_t::x_max = 0;
      const double LHAPDF_DF_t::Q2_min = 0;
      const double LHAPDF_DF_t::Q2_max = 0;

      LHAPDF_DF_t::LHAPDF_DF_t( std::string DF_name, std::string weight_option, int nset_in, int member_in,
                                std::string pT_instructions ) :
         FlavArrayFunc_t( std::string("pT ") + pT_instructions ), nset(nset_in), member(member_in)
      {
         // make msg
         {
            std::stringstream ss;
            ss << DF_name << ": using LHAPDF " << nset << ", " << member << " w/ weight opt. '" << weight_option << "'"; 
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
            throw Error::Constructing( "LHAPDF_DF_t", "x_min >= x_max.  Did you forget to initialize LHA_PDF?");

         if( Q2_min >= Q2_max )
            throw Error::Constructing( "LHAPDF_DF_t", "Q2_min >= Q2_max.  Did you forget to initialize LHA_PDF?");

         double w_init = 0;
         if( DF_name == "f1" || weight_option == "1" )
            w_init = 1;

         for( double *p = &weight[0]; p != &weight[11]; ++p )
            (*p) = w_init;

         //cerr << DF_name << ", w init = " << w_init << endl;

         if( DF_name.at(0) == 'g' || DF_name.at(0) == 'h' ){
            if ( weight_option == "B" ){
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = -0.3;
            } else if ( weight_option == "C" ){
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = -0.3;
               weight[ STR_FLAV ] = -0.3;
               weight[ ANTI_STR_FLAV ] = -0.3;
            } else if ( weight_option == "D" ){
               weight[ UP_FLAV ] = -0.5;
               weight[ DOWN_FLAV ] = 0.5;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = 0.5;
               weight[ STR_FLAV ] = 0.5;
               weight[ ANTI_STR_FLAV ] = 0.5;
            } else if ( weight_option == "E" ){
               weight[ UP_FLAV ] = -0.3;
               weight[ DOWN_FLAV ] = -0.3;
               weight[ ANTI_UP_FLAV ] = 0.5;
               weight[ ANTI_DOWN_FLAV ] = 0.5;
               weight[ STR_FLAV ] = 0.5;
               weight[ ANTI_STR_FLAV ] = 0.5;
            } else {
               weight[ UP_FLAV ] = 0.5;
               weight[ DOWN_FLAV ] = -0.3;
            };
         };

         // set pointers
         func_ptr[ UP_FLAV ]        = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Up );
         func_ptr[ DOWN_FLAV ]      = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Down );
         func_ptr[ STR_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Strange );
         func_ptr[ CHM_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Charm );
         func_ptr[ BOT_FLAV ]       = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Bottom );
         func_ptr[ ANTI_UP_FLAV ]   = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Anti_Up );
         func_ptr[ ANTI_DOWN_FLAV ] = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Anti_Down );
         func_ptr[ ANTI_STR_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Anti_Strange );
         func_ptr[ ANTI_CHM_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Anti_Charm );
         func_ptr[ ANTI_BOT_FLAV ]  = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Anti_Bottom );
         func_ptr[ GLUON_FLAV ]     = static_cast< double (FlavArrayFunc_t::*)( const Var_t& var ) const >( &LHAPDF_DF_t::Gluon );

      };

      LHAPDF_DF_t::~LHAPDF_DF_t(){
         // nothing to do
      };

      // To query allowed values for variables
      void LHAPDF_DF_t::GetVarRange( Var_t& var_min, Var_t& var_max ) const{
         var_min.Q2 = Q2_min;
         var_max.Q2 = Q2_max;
         var_min.x = x_min;
         var_max.x = x_max;
         pT_func->GetRange( var_min.pT, var_max.pT );
      };

      // flags with 1 each relevant variable used
      void LHAPDF_DF_t::GetRelevantVar( Var_t& var ) const {
         var.x = 1;
         var.Q2 = 1;
         var.pT = 1;
      };

      double LHAPDF_DF_t::Up( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::UP >( weight[UP_FLAV], var );
      };

      double LHAPDF_DF_t::Down( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::DOWN >( weight[DOWN_FLAV], var );
      };

      double LHAPDF_DF_t::Strange( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::STRANGE >( weight[STR_FLAV], var );
      };

      double LHAPDF_DF_t::Charm( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::CHARM >( weight[CHM_FLAV], var );
      };

      double LHAPDF_DF_t::Bottom( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::BOTTOM >( weight[BOT_FLAV], var );
      };

      double LHAPDF_DF_t::Anti_Up( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::UBAR >( weight[ANTI_UP_FLAV], var );
      };

      double LHAPDF_DF_t::Anti_Down( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::DBAR >( weight[ANTI_DOWN_FLAV], var );
      };

      double LHAPDF_DF_t::Anti_Strange( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::SBAR >( weight[ANTI_STR_FLAV], var );
      };

      double LHAPDF_DF_t::Anti_Charm( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::CBAR >( weight[ANTI_CHM_FLAV], var );
      };

      double LHAPDF_DF_t::Anti_Bottom( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::BBAR >( weight[ANTI_BOT_FLAV], var );
      };

      double LHAPDF_DF_t::Gluon( const Var_t& var ) const{
         LHAPDF::usePDFMember (nset, member);
         return LHAPDF_fx< LHAPDF::GLUON >( weight[GLUON_FLAV], var );
      };

      double LHAPDF_DF_t::Quark( const Var_t& var ) const {
         return Undefined( var );
      };

      double LHAPDF_DF_t::Sea( const Var_t& var ) const {
         return Undefined( var );
      };

   };
};

#endif
