#include "XSec.d/SIDIS_1had.h"

#include "Common.d/Consts.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/GetMassSquared.h"
#include "Common.d/SetPID.h"
#include "Common.d/Yfunctions.h"

#include "XSec_Term.d/SIDIS_1had_f1D1_.h"
#include "XSec_Term.d/SIDIS_1had_cos2phi.h"
#include "XSec_Term.d/SIDIS_1had_cosphi.h"
//#include "XSec_Term.d/Hermes_Auu_1had.h"

// Independent variables are: x, y, z, P_hperp, phi_h, pT, phi_pT, psi in var

namespace TMDGen {
   namespace XSec {

      // Note: does not handle the unusual case where a specific given DF or FF
      // restricts the range, yet is not actually used in the cross section
      void SIDIS_1had_t::DetermineDomain( const sgUtil::ParseInputReturn_t& parsed_input ){
         sgUtil::ParseInputReturn_const_iterator_t it;

         // Variable Ranges

         if( DF_Set->GetVarRange( min, max ) )
            throw Error::Constructing( "SIDIS_1had", "error determining variable range for DF set" );


         if( FF_Set->GetVarRange( min, max ) ){
            throw Error::Constructing( "SIDIS_1had", "error determining variable range for FF set" );
         };

         GetRange( parsed_input, "Q2", min.Q2, max.Q2 );
         GetRange( parsed_input, "W2", min.W2, max.W2 );
         GetRange( parsed_input, "x", min.x, max.x );
         GetRange( parsed_input, "y", min.y, max.y );
         GetRange( parsed_input, "z", min.z, max.z );
         GetRange( parsed_input, "pT", min.pT, max.pT );
         GetRange( parsed_input, "kT", min.kT, max.kT );
         GetRange( parsed_input, "P_hperp", min.P_hperp, max.P_hperp );
      };

      // compute other quantities in var
      int SIDIS_1had_t::ComputeOtherVar( Var_t& var ){
         var.had.M = had_M;

         return var.ComputeOtherVar_SIDIS_1had( min, max );
      };

      double SIDIS_1had_t::Overall_Kinematic_Factor( const Var_t& var ) const {
         double denom =  var.x * var.y * var.Q2;
         double paren = 0.5 * var.gamma_sq / var.x + 1;

         return var.P_hperp * var.pT * paren / denom;
      };

      void SIDIS_1had_t::MakeTerms( const sgUtil::ParseInputReturn_t& parsed_input ){
         sgUtil::ParseInputReturn_const_iterator_t it;

         // sanity check
         // clear all terms
         while( !xsec_term.empty() ){
            delete xsec_term.back();
            xsec_term.pop_back();
         };

         // angular integrated
         angular_integrated_xsec_term = new XSec_Term::SIDIS_1had_f1D1_t( pol_UU, &yFunc_array[yFunctions::A_idx], DF_Set, FF_Set ); 
         yFunc_used[yFunctions::A_idx] = 1;

         if( angular_integrated_xsec_term->NonZero() ){
            xsec_term.push_back( angular_integrated_xsec_term );
         } else {
            delete angular_integrated_xsec_term;
            throw Error::SanityCheckFailure("SIDIS_1had_t::MakeTerms(...)",
                                            "Angular integrated term in cross section is zero");
         };

         XSec_Term::XSec_Term_t* temp_term;

         it = parsed_input.find("No_Angular_Dependence");
         if( it == parsed_input.end() || it->second == "false" || it->second == "False" || it->second == "FALSE" ){
            // check first for overall model
            it = parsed_input.find("Auu");
            if( it != parsed_input.end() ){
               if( it->second == "HERMES" ){
                  // try to put DF and FF terms together to make unpolarized moments
                  throw Error::Constructing( "SIDIS_1had_t", "Hermes Auu moments not yet fully implemented." );

//                   temp_term = new XSec_Term::Hermes_Auu_1had_t( had_PID );
//                   if( temp_term->NonZero() )
//                      xsec_term.push_back( temp_term );
               } else {
                  throw Error::Constructing( "SIDIS_1had_t", "Invalid 'Auu' directive" );
               };
            } else {
               // cos2phi moment (twist 2 contribution)
               temp_term = new XSec_Term::SIDIS_1had_cos2phi_t( pol_UU, &yFunc_array[yFunctions::B_idx], DF_Set, FF_Set ); 
               if( temp_term->NonZero() ){
                  yFunc_used[yFunctions::B_idx] = 1;
                  xsec_term.push_back( temp_term );
               } else {
                  delete temp_term;
               };

               if( twist > 2 ){
                  // twist three cos phi term
                  temp_term = new XSec_Term::SIDIS_1had_cosphi_t( pol_UU, &yFunc_array[yFunctions::C_idx], DF_Set, FF_Set ); 
                  if( temp_term->NonZero() ){
                     yFunc_used[yFunctions::C_idx] = 1;
                     xsec_term.push_back( temp_term );
                  } else {
                     delete temp_term;
                  };
               };
            };

            if( beam_pol != Enum::UNPOL )
               throw Error::Constructing( "SIDIS_1had_t", "Not yet programmed polarized beam cross section for single hadrons" );

            if( target_pol != Enum::UNPOL )
               throw Error::Constructing( "SIDIS_1had_t", "Not yet programmed polarized target cross section for single hadrons" );
         };
      };

      SIDIS_1had_t::SIDIS_1had_t( const sgUtil::ParseInputReturn_t& parsed_input, Var_t& min_in, Var_t& max_in ) :
         XSec_t( parsed_input, min_in, max_in )
      {
         sgUtil::ParseInputReturn_const_iterator_t it;

         // FINAL STATE

         it = parsed_input.find("Process");
         if( it == parsed_input.end() ){
            throw Error::SanityCheckFailure( "XSec::SIDIS_1had", "no 'Process' directive given" );
         } else if ( it->second != "SIDIS" ){
            throw Error::SanityCheckFailure( "XSec::SIDIS_1had", "'Process' directive is not 'SIDIS'" );
         };

         it = parsed_input.find("Final_State");
         if( it == parsed_input.end() )
            throw Error::SanityCheckFailure( "XSec::SIDIS_1had", "no 'Final_State' directive given" );

         if( it->second == "Single Hadron" ){ 
            //final_state = SINGLE_HADRON;

            it = parsed_input.find("Hadron_PID");
            if( it == parsed_input.end() )
               throw Error::Constructing( "XSec::SIDIS_1had", "no 'Hadron_PID' directive given" );

            if( SetPID( it->second, had_PID ) )
               throw Error::Constructing( "XSec::SIDIS_1had", std::string("invalid Hadron_PID: '") + it->second + "'" );

            // set mass squared
            had_M = sqrt( GetMassSquared( had_PID ) );

         } else {
            throw Error::SanityCheckFailure( "XSec::SIDIS_1had", "final state is not 'Single Hadron'" );
         };

         // OVERALL FACTOR

         // overall factor is based on Alessandro, Deihl, et al paper, w/o the A(y) function and times 1/(hbar c^2) for units
         *const_cast< double* >( &overall_const_factor ) = ALPHA_EM * ALPHA_EM * HBARC2 * TWO_PI;

         FinishConstructing( parsed_input );
      };

      SIDIS_1had_t::~SIDIS_1had_t(){
         // nothing to do
      };

      LundPID_t SIDIS_1had_t::Get_Had1_PID() const {
         return had_PID;
      };

   };
};
