#include "XSec.d/SIDIS_2had.h"

#include "Common.d/Consts.h"
#include "Common.d/Exceptions.h"
#include "Common.d/GetMassSquared.h"
#include "Common.d/SetPID.h"
#include "Common.d/Yfunctions.h"

#include "XSec_Term.d/XSec_Term.h"
#include "XSec_Term.d/SIDIS_2had_f1D1_.h"
#include "XSec_Term.d/SIDIS_2had_cosphi.h"
#include "XSec_Term.d/SIDIS_2had_cos2phi.h"
#include "XSec_Term.d/SIDIS_2had_Sivers.h"
#include "XSec_Term.d/SIDIS_2had_Collins.h"
//#include "XSec_Term.d/SIDIS_2had_Pretzel.h"
#include "XSec_Term.d/SIDIS_2had_Naive.h"

#ifndef MAX_L
#define MAX_L 2
#endif

#include <iostream>
#include <sstream>

// Independent variables are: x, y, z, P_hperp, phi_h, pT, phi_pT, psi, M_h, cos_vartheta, phi_R

namespace TMDGen {
   namespace XSec {

      // Note: does not handle the unusual case where a specific given DF or FF
      // restricts the range, yet is not actually used in the cross section
      void SIDIS_2had_t::DetermineDomain( const sgUtil::ParseInputReturn_t& parsed_input ){
         sgUtil::ParseInputReturn_const_iterator_t it;

         // Variable Ranges

         if( DF_Set->GetVarRange( min, max ) )
            throw Error::Constructing( "XSec::SIDIS_2had_t(...)", "error determining variable range for DF set" );


         if( FF_Set->GetVarRange( min, max ) ){
            throw Error::Constructing( "XSec::SIDIS_2had_t(...)", "error determining variable range for FF set" );
         };

         GetRange( parsed_input, "Q2", min.Q2, max.Q2 );
         GetRange( parsed_input, "W2", min.W2, max.W2 );
         GetRange( parsed_input, "x", min.x, max.x );
         GetRange( parsed_input, "y", min.y, max.y );
         GetRange( parsed_input, "z", min.z, max.z );
         GetRange( parsed_input, "pT", min.pT, max.pT );
         GetRange( parsed_input, "kT", min.kT, max.kT );
         GetRange( parsed_input, "P_hperp", min.P_hperp, max.P_hperp );
         GetRange( parsed_input, "M_h", min.had_0.M, max.had_0.M );

         if( min.had_0.M < had_1_M + had_2_M ){
            std::cerr << "\tWARNING: minimum dihadron invarient mass below threshold." << std::endl;
            std::cerr << "\tResetting to threshold: " << had_1_M + had_2_M << "." << std::endl << std::endl;
            min.had_0.M = had_1_M + had_2_M;
         };

      };

      // compute other quantities in var
      int SIDIS_2had_t::ComputeOtherVar( Var_t& var ) {
         var.had_1.M = had_1_M;
         var.had_2.M = had_2_M;


         int ierr = var.ComputeOtherVar_SIDIS_2had( min, max );

         // cos array in order |0,0>, |1,1>, |1,0>, |2,2>, |2,1>, |2,0>

         cos_array[Enum::P_00] = 1;
         cos_array[Enum::P_11] = var.sin_vartheta;   // sin(vartheta)
         cos_array[Enum::P_10] = var.cos_vartheta;
         cos_array[Enum::P_22] = 3*var.sin_vartheta*var.sin_vartheta;
         cos_array[Enum::P_21] = 3*var.cos_vartheta*var.sin_vartheta;  // i.e. 3 cos(theta)sin(theta) = 3/2 sin(2 theta) 
         cos_array[Enum::P_20] = 0.5*(3.*var.cos_vartheta*var.cos_vartheta - 1.);


//          std::cout << (var.vartheta) << ' ';
//          for( int i=0; i<6; ++i )
//             std::cout << cos_array[i] << ' ';
//          std::cout << std::endl;

         return ierr;

      };

      double SIDIS_2had_t::Overall_Kinematic_Factor( const Var_t& var ) const {
         double denom =  var.x * var.y * var.Q2;
         double paren = 0.5 * var.gamma_sq / var.x + 1;

         // Sept 20--tried taking out the var.P_hperp to see the effect -- subtle, but hurts shape of P_hperp distr.
         return var.P_hperp * var.pT * var.had_0.M * paren / denom; // 
      };

      void SIDIS_2had_t::MakeTerms( const sgUtil::ParseInputReturn_t& parsed_input ){
         sgUtil::ParseInputReturn_const_iterator_t it;

         // sanity check
         // clear all terms
         while( !xsec_term.empty() ){
            delete xsec_term.back();
            xsec_term.pop_back();
         };

         XSec_Term::XSec_Term_t* temp_term;

         // angular integrated
         // already checked the sets include "f1" and "D1"
         angular_integrated_xsec_term = new XSec_Term::SIDIS_2had_f1D1_t( 0, 0, pol_UU,
                                                                          &yFunc_array[yFunctions::A_idx],
                                                                          cos_array, DF_Set, FF_Set ); 
         yFunc_used[yFunctions::A_idx] = 1;

         if( angular_integrated_xsec_term->NonZero() ){
            xsec_term.push_back( angular_integrated_xsec_term );
         } else {
            delete angular_integrated_xsec_term;
            throw Error::SanityCheckFailure("SIDIS_2had_t::MakeTerms(...)",
                                            "Angular integrated term in dihadron cross section is zero");
         };

         if( twist > 3 ){
            throw Error::NotYetProgrammed("XSec::SIDIS_2had with twist > 3");
         };

         // todo: check first for overall model
         it = parsed_input.find("No_Angular_Dependence");
         if( it == parsed_input.end() || it->second == "false" || it->second == "False" || it->second == "FALSE" ){

            // make unpolarized angular terms
            for( int l = 0; l <= MAX_L; ++l ){
               for( int m = -l; m <= l; ++m ){

                  // don't duplicate the angular integrated term, nor include negative m states
                  // as states are really real combinations anyhow
                  if( l && m > -1 ){
                     temp_term = new XSec_Term::SIDIS_2had_f1D1_t( l, m, pol_UU,
                                                                   &yFunc_array[yFunctions::A_idx],
                                                                   cos_array, DF_Set, FF_Set );

                     //std::cerr << l << ' ' << m << " : f1 D1 term zero " << !temp_term->NonZero() << ' ' << temp_term << std::endl;
                     //throw Error::SanityCheckFailure("SIDIS_2had_t::MakeTerms(...)", "UGABOOGA" );

                     if( temp_term->NonZero() ){
                        xsec_term.push_back( temp_term );
                        yFunc_used[yFunctions::A_idx] = 1;
                     } else {
                        delete temp_term;
                     };
                  };

                  temp_term = new XSec_Term::SIDIS_2had_cos2phi_t( l, m, pol_UU, &yFunc_array[yFunctions::B_idx], cos_array, DF_Set, FF_Set );

                  //std::cerr << l << ' ' << m << " : cos2phi term zero " << !temp_term->NonZero() << '\t' << temp_term << std::endl;

                  if( temp_term->NonZero() ){
                     xsec_term.push_back( temp_term );
                     yFunc_used[yFunctions::B_idx] = 1;
                  } else {
                     delete temp_term;
                  };

                  if( twist > 2 ){
                     temp_term = new XSec_Term::SIDIS_2had_cosphi_t( l, m, pol_UU,
                                                                     &yFunc_array[yFunctions::V_idx], cos_array, DF_Set, FF_Set );

                     //std::cerr << l << ' ' << m << " : cosphi term zero " << !temp_term->NonZero() << '\t' << temp_term << std::endl;

                     if( temp_term->NonZero() ){
                        xsec_term.push_back( temp_term );
                        yFunc_used[yFunctions::V_idx] = 1;
                     } else {
                        delete temp_term;
                     };
                  };
               };
            };

            if( beam_pol != Enum::UNPOL )
               throw Error::Constructing( "XSec::SIDIS_2had_t(...)", "Not yet programmed polarized beam cross section for dihadrons" );

            if( target_pol == Enum::TRANS ){

               if( twist > 2 )
                  std::cerr << "WARNING: Twist > 2 not yet programmed for UT moments" << std::endl;


               // make trans target angular terms
               for( int l = 0; l <= MAX_L; ++l ){
                  for( int m = -l; m <= l; ++m ){
                     temp_term = new XSec_Term::SIDIS_2had_Sivers_t( l, m, pol_UT,
                                                                     &yFunc_array[yFunctions::A_idx], cos_array, DF_Set, FF_Set );
                     if( temp_term->NonZero() ){
                        xsec_term.push_back( temp_term );
                        yFunc_used[yFunctions::A_idx] = 1;
                     } else {
                        delete temp_term;
                     };

                     temp_term = new XSec_Term::SIDIS_2had_Collins_t( l, m, pol_UT,
                                                                      &yFunc_array[yFunctions::B_idx], cos_array, DF_Set, FF_Set );
                     if( temp_term->NonZero() ){
                        xsec_term.push_back( temp_term );
                        yFunc_used[yFunctions::B_idx] = 1;
                     } else {
                        delete temp_term;
                     };

                     /*
                     temp_term = new XSec_Term::SIDIS_2had_Pretzel( l, m, pol_UT, &yFunc_array[yFunctions::B_idx], cos_array, DF_Set, FF_Set );
                     if( temp_term->NonZero() ){
                        xsec_term.push_back( temp_term );
                        yFunc_used[yFunctions::B_idx] = 1;
                     } else {
                        delete temp_term;
                     };
                     */

                     //if( twist > 2 )
                     //   throw Error::Constructing( "XSec::SIDIS_2had_t(...)",
                     //                              "Not yet programmed trans. target cross section for dihadrons for twist > 2" );

                     //std::cerr << "WARNING: Twist > 2 not yet programmed for UT moments" << std::endl;
                  };
               };

               if( beam_pol != Enum::UNPOL )
                  throw Error::Constructing( "XSec::SIDIS_2had_t(...)",
                                             "Not yet programmed polarized beam, trans. target cross section for dihadrons" );

            } else if (target_pol == Enum::LONG ){
               throw Error::Constructing( "XSec::SIDIS_2had_t(...)",
                                          "Not yet programmed long. polarized target cross section for dihadrons" );

               if( beam_pol != Enum::UNPOL )
                  throw Error::Constructing( "XSec::SIDIS_2had_t(...)",
                                             "Not yet programmed polarized beam, long. target cross section for dihadrons" );
            } else {
               if( beam_pol != Enum::UNPOL )
                  throw Error::Constructing( "XSec::SIDIS_2had_t(...)",
                                             "Not yet programmed polarized beam, unpolarized target cross section for dihadrons" );

            }; // beam pol

            // Search for Naive model type terms

            // unpolarized
            for( int i=1; i<24; ++i ){
               std::stringstream ss;
               ss << "SIDIS_2had_Naive_XSec_Term_UU_" << i;

               it = parsed_input.find(ss.str());
               if( it != parsed_input.end() ){
                  std::string pol = "UU";
                  std::string params = it->second;
                  temp_term = new XSec_Term::SIDIS_2had_Naive_t( pol, i, params, pol_UU, yFunc_array, cos_array );

                  if( temp_term->NonZero() ){
                     xsec_term.push_back( temp_term );
                     yFunc_used[yFunctions::B_idx] = 1;  // don't know which, so just set both on
                     yFunc_used[yFunctions::V_idx] = 1;
                  } else {
                     delete temp_term;
                  };
               };
            };

            // polarized
            for( int i=0; i<27; ++i ){
               std::stringstream ss;
               ss << "SIDIS_2had_Naive_XSec_Term_UT_" << i;

               it = parsed_input.find(ss.str());
               if( it != parsed_input.end() ){
                  std::string pol = "UT";
                  std::string params = it->second;
                  temp_term = new XSec_Term::SIDIS_2had_Naive_t( pol, i, params, pol_UU, yFunc_array, cos_array );

                  if( temp_term->NonZero() ){
                     xsec_term.push_back( temp_term );
                     yFunc_used[yFunctions::A_idx] = 1;  // don't know which, so just set both on
                     yFunc_used[yFunctions::B_idx] = 1;
                  } else {
                     delete temp_term;
                  };
               };
            };

         };  // whether to include angular dependence
      };




      SIDIS_2had_t::SIDIS_2had_t( const sgUtil::ParseInputReturn_t& parsed_input, Var_t& min_in, Var_t& max_in ) :
         XSec_t( parsed_input, min_in, max_in )
      {
         sgUtil::ParseInputReturn_const_iterator_t it;

         // FINAL STATE

         it = parsed_input.find("Process");
         if( it == parsed_input.end() ){
            throw Error::SanityCheckFailure( "XSec::SIDIS_2had", "no 'Process' directive given" );
         } else if ( it->second != "SIDIS" ){
            throw Error::SanityCheckFailure( "XSec::SIDIS_2had", "'Process' directive is not 'SIDIS'" );
         };

         it = parsed_input.find("Final_State");
         if( it == parsed_input.end() )
            throw Error::SanityCheckFailure( "XSec::SIDIS_2had", "no 'Final_State' directive given" );

         if( it->second == "Dihadron" ){ 
            //final_state = Enum::DIHADRON;

            it = parsed_input.find("Hadron_1_PID");
            if( it == parsed_input.end() )
               throw Error::Constructing( "XSec::SIDIS_2had", "no 'Hadron_1_PID' directive given" );

            if( SetPID( it->second, had_1_PID ) )
               throw Error::Constructing( "XSec::SIDIS_2had", std::string("invalid Hadron_1_PID: '") + it->second + "'" );

            it = parsed_input.find("Hadron_2_PID");
            if( it == parsed_input.end() )
               throw Error::Constructing( "XSec::SIDIS_2had", "no 'Hadron_2_PID' directive given" );

            if( SetPID( it->second, had_2_PID ) )
               throw Error::Constructing( "XSec::SIDIS_2had", std::string("invalid Hadron_2_PID: '") + it->second + "'" );

            // set mass squared
            had_1_M = sqrt( GetMassSquared( had_1_PID ) );
            had_2_M = sqrt( GetMassSquared( had_2_PID ) );

         } else {
            throw Error::SanityCheckFailure( "XSec::SIDIS_2had", "final state is not 'Dihadron'" );
         };

         // OVERALL FACTOR
         // This is the constant part.
         // The kinematic part is in Overall_Kinematic_Factor

         // one two pi cancels 2pi in hbar, one 2pi is in cross section, and the other comes with the dphi_pT
         *const_cast< double* >( &overall_const_factor ) = ALPHA_EM * ALPHA_EM * HBARC2 / TWO_PI / TWO_PI / TWO_PI;


         FinishConstructing( parsed_input );
      };

      SIDIS_2had_t::~SIDIS_2had_t(){
         // nothing to do
      };

      LundPID_t SIDIS_2had_t::Get_Had1_PID() const {
         return had_1_PID;
      };

      LundPID_t SIDIS_2had_t::Get_Had2_PID() const {
         return had_2_PID;
      };

   };
};
