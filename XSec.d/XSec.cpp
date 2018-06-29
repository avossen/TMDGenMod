/*
  Main class for constructing Cross Sections
*/

#define DO_OUT   0
#define DO_OUT_2 0

#include "XSec.d/XSec.h"

#include "Common.d/Consts.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/LundPID.h"
#include "Common.d/ParseInput.h"
#include "Common.d/QuarkChargeSquared.h"
#include "Common.d/SetPID.h"
#include "Common.d/Yfunctions.h"

#include "XSec_Term.d/XSec_Term.h"

#include "DF.d/Full_DF_Set.h"
#include "FF.d/Full_FF_Set.h"

#include <vector>
#include <string>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <cstdlib>

namespace TMDGen {
   namespace XSec {

      // constructs all DF and FF models and calls function to determine range and make terms
      XSec_t::XSec_t( const sgUtil::ParseInputReturn_t& parsed_input, Var_t& min_in, Var_t& max_in ) :
         overall_const_factor(0), min(min_in), max(max_in), skip_overall_factors(0)
      {
         sgUtil::ParseInputReturn_const_iterator_t it;

         yFunc_used[0] = 0;
         yFunc_used[1] = 0;
         yFunc_used[2] = 0;
         yFunc_used[3] = 0;
         yFunc_used[4] = 0;

         // TARGET STATE

         it = parsed_input.find("Target_Pol");
         if( it == parsed_input.end() ){
            throw Error::Constructing( "XSec_t", "no 'Target_Pol' directive given" );
         } else if( it->second == "L" ){
            target_pol = Enum::LONG;
         } else if( it->second == "T" ){
            target_pol = Enum::TRANS;
         } else if( it->second == "U" ){
            target_pol = Enum::UNPOL;
         } else {
            throw Error::Constructing( "XSec_t", "invalid 'Target_Pol' directive: must be 'U', 'L', or 'T'" );
         };

         // BEAM STATE

         it = parsed_input.find("Beam_Pol");
         if( it == parsed_input.end() ){
            throw Error::Constructing( "XSec_t", "no 'Beam_Pol' directive given" );
         } else if( it->second == "L" ){
            beam_pol = Enum::LONG;
         } else if( it->second == "T" ){
            beam_pol = Enum::TRANS;
         } else if( it->second == "U" ){
            beam_pol = Enum::UNPOL;
         } else {
            throw Error::Constructing( "XSec_t", "invalid 'Beam_Pol' directive: must be 'U', 'L', or 'T'" );
         };


         // TWIST

         it = parsed_input.find("Twist");
         if( it == parsed_input.end() ){
            throw Error::Constructing( "XSec_t", "no 'Twist' directive given" );
         };

         twist = atoi( it->second.data() );
         if( twist < 2 )
            throw Error::Constructing( "XSec_t", "Given 'Twist' directive < 2, yet leading order is twist 2" );
         if( twist > 3 )
            throw Error::Constructing( "XSec_t", "Have only programmed twist 2 and 3" );


         // NOW START ALLOCATING MODELS

         // try allocating all distribution functions
         try{
            DF_Set = new DF::Full_DF_Set_t( parsed_input );
         }
         catch( std::exception& e ){
            throw Error::Constructing( "XSec_t", std::string("error in required class DF::Full_DF_Set_t:\n\t\t") + e.what() );
         };
      
         // make sure "f1" exists
         if( !DF_Set->Includes("f1") ){
            delete DF_Set;
            throw Error::Constructing( "XSec_t", "No valid 'f1' directive given" );
         };

         // try allocating all fragmentation functions
         try{
            FF_Set = new FF::Full_FF_Set_t( parsed_input );
         }
         catch( std::exception& e ){
            delete DF_Set;
            throw Error::Constructing( "XSec_t", std::string("error in required class FF::Full_FF_Set_t:\n\t\t") + e.what() );
         };

         // make sure "D1" exists
         if( !FF_Set->Includes("D1") && !FF_Set->Includes("Re_D1_00") ){
            delete DF_Set;
            delete FF_Set;
            throw Error::Constructing( "XSec_t", "No valid 'D1' directive given" );
         };

         it = parsed_input.find("Skip_Overall_Factors");
         if( it != parsed_input.end() )
            skip_overall_factors = 1;

      };

      // Would like to have this in the constructor, but
      // it needs to call pure virtual function `DetermineDomain'
      void XSec_t::FinishConstructing( const sgUtil::ParseInputReturn_t& parsed_input ){

         sgUtil::ParseInputReturn_const_iterator_t it;

         // set min and max values, based on models and input
         try{ 
            DetermineDomain( parsed_input );
         }
         catch( std::exception& e ){
            throw Error::Constructing( "XSec_t", std::string("while determing domain, caught error:\n\t\t") + e.what() );
         };

         // Make all terms
         try{ 
            MakeTerms( parsed_input );
         }
         catch( std::exception& e ){
            while( !xsec_term.empty() ){
               delete xsec_term.back();
               xsec_term.pop_back();
            };
            throw Error::Constructing( "XSec_t", std::string("while making terms, caught error:\n\t\t") + e.what() );
         };

         it = parsed_input.find("Remove_All_Moments_Except");
         if( it != parsed_input.end() ){
            cerr << "\t\tRemoving All Moments Except ang. int. and '" << it->second << "'" << endl;

            std::string term_to_keep_name = it->second;
            XSec_Term::XSec_Term_t* term_to_keep = 0;

            for( size_t i = 0; i < xsec_term.size(); ++i )
               if( xsec_term[i]->Message() == term_to_keep_name )
                  term_to_keep = xsec_term[i];
            
            if( !term_to_keep )
               throw Error::Constructing( "XSec_t", std::string("Told to only keep term '") + term_to_keep_name + "', but this term not in cross section" );

            xsec_term.resize(2);
            xsec_term[0] = angular_integrated_xsec_term;
            xsec_term[1] = term_to_keep;
         };

         cerr << "\tTERMS IN CROSS SECTION" << endl;
         for( size_t i = 0; i < xsec_term.size(); ++i )
            cerr << '\t' << i << ' ' << xsec_term[i]->Message() << endl;


         // OTHER PARAMETERS


      };

      // deconstructor
      XSec_t::~XSec_t(){
         delete DF_Set;
         delete FF_Set;

         // clear all terms
         while( !xsec_term.empty() ){
            delete xsec_term.back();
            xsec_term.pop_back();
         };
      };

      // evaluate 
      double XSec_t::Eval_Inner( Var_t& var ) {
         double val = 0;

//          if( !var.integrating )
//                cout << "A" << endl;

         // another routine should have already set correct var.lep_E_in

         // compute dependent variables
         if( !ComputeOtherVar( var ) ){

//             if( !var.integrating )
//                cout << "B " << val << endl;

            if( yFunc_used[0] )
               yFunc_array[0] = yFunctions::A( var );
            if( yFunc_used[1] )
               yFunc_array[1] = yFunctions::B( var );
            if( yFunc_used[2] )
               yFunc_array[2] = yFunctions::C( var );
            if( yFunc_used[3] )
               yFunc_array[3] = yFunctions::V( var );
            if( yFunc_used[4] )
               yFunc_array[4] = yFunctions::W( var );

//             if( !var.integrating )
//                cerr << "**\tCalling Precompute" << endl;

            FF_Set->Precompute( var );
            DF_Set->Precompute( var );

            pol_UU = QuarkChargeSquared[ var.flavor ]; 
            pol_UT = pol_UU * var.S_T;                      // todo: the rest

//             if( !var.integrating )
//                cout << "pol_UT = " << pol_UU << " * " << var.S_T << endl;

            ang_int_xsec = angular_integrated_xsec_term->Eval( var );
            val = ang_int_xsec;

//             if( !var.integrating )
//                cout << "C " << val << ' ' << angular_integrated_xsec_term->Eval( var ) << endl;

            if( val != val ){
               cerr << "Q2 = " << var.Q2 << endl;
               throw Error::SanityCheckFailure( "XSec_t::Eval_inner(...)", "angular integrated is NaN" );
            };


            if( val < 0 ){
               cerr << "WARNING: Angular integrated cross section is negative: " << ang_int_xsec << endl;
               //cerr << var.x << ' ' << var.y << ' ' << var.z << ' ' << var.P_hperp << endl;
               //cerr << var.had_0.M << ' ' << var.pT << ' ' << var.kT << endl;
               //throw Error::SanityCheckFailure("XSec_t::operator()", "angular integrated cross section is negative");

               return 0;
            };

            //std::cout << val << ' ' << angular_integrated_xsec_term << endl;

            if( DO_OUT && !var.integrating )
               std::cout << var.flavor << ' ' << val;

            if( DO_OUT_2 && !var.integrating )
               std::cout << angular_integrated_xsec_term->Message() << ' ' << val << std::endl;

            // iterate over other terms
            std::vector< XSec_Term::XSec_Term_t* >::const_iterator it = xsec_term.begin();
            for( ; it != xsec_term.end(); ++it ){

               // don't do angular integrated twice!
               if( *it != angular_integrated_xsec_term ){
                  double last_val = val;

                  //double temp = (*it)->Factored() ? ang_int_xsec * (*it)->Eval( var ) : (*it)->Eval( var );
                  //val += temp;
                  //std::cout << "\t___ " << temp << " ___" << endl;

                  double temp = (*it)->Eval( var );

                  if( (*it)->Factored() )
                     temp *= ang_int_xsec;

                  val += temp;

                  if( DO_OUT_2 && !var.integrating )
                     std::cout << '\t' << (*it)->Message() << ' ' << val-last_val << std::endl;

                  if( DO_OUT && !var.integrating ){
                     std::cout << ' ' << val-last_val;

                     //if( !((it-xsec_term.begin()) % 6 ) )
                     //   std::cout << std::endl;
                  };
               };
            };

            if( DO_OUT && !var.integrating )
               std::cout << " | " << val;
         };

//          if( !var.integrating )
//             cout << "D " << val << endl;

         if( val != val ){
            cerr << "Q2 = " << var.Q2 << endl;
            throw Error::SanityCheckFailure( "XSec_t::Eval_inner(...)", "returning NaN" );
         };

         return val;
      };

      // evaluate 
      double XSec_t::operator() ( Var_t& var ) {
         Eval( var );

         return var.XS;
      };

      // evaluate 
      void XSec_t::Eval( Var_t& var ) {
         double val     = Eval_Inner( var );

         var.XS = 0;
         var.XS_int = 0;

         if( DO_OUT && !var.integrating && val )
            std::cout << " kinfactor = " << Overall_Kinematic_Factor( var ) << std::endl;

         if( skip_overall_factors )
            throw 0;

//          if( !var.integrating )
//             cout << "y " << val << ' ' << skip_overall_factors << endl;

         if( val && !skip_overall_factors ){
            double prefactor = overall_const_factor;
            prefactor *= Overall_Kinematic_Factor( var );

//             if( !var.integrating )
//                cout << "x " << val << ' ' << var.XS_int << ' ' << prefactor << endl;

            var.XS = val;
            var.XS *= prefactor;

            var.XS_int = ang_int_xsec;
            var.XS_int *= prefactor;
         };
      };

      // evaluate, but divide by angular integrated cross section
      double XSec_t::Eval_only_Angular( Var_t& var, double& ang_int_xsec_out ){
         double val = Eval_Inner( var );
         if( !ang_int_xsec )
            val = 0;

         //std::cout << val;

         if( val )
            val /= ang_int_xsec;

         ang_int_xsec_out = ang_int_xsec;

         //std::cout << " / " << ang_int_xsec << " = " << val << endl;

         return val;
      };


      void XSec_t::GetRange( const sgUtil::ParseInputReturn_t& parsed_input, std::string var_name, double& var_min, double& var_max ){
         double given_min = -999;
         double given_max = -999;

         sgUtil::ParseInputReturn_const_iterator_t it;

         it = parsed_input.find( var_name + "_min");
         if( it != parsed_input.end() ){
            given_min = atof( it->second.data() );
         };

         it = parsed_input.find( var_name + "_max");
         if( it != parsed_input.end() ){
            given_max = atof( it->second.data() );
         };

         //cerr << "\tChecking range for " << var_name << ' ' << var_min << ' ' << var_max << ' ' << given_min << ' ' << given_max << endl;
         if( var_min == var_max ){
            // i.e. has not yet been set
            var_min = given_min;
            var_max = given_max;

            if( var_min == -999 || var_max == -999 )
               throw Error::Constructing( "XSec_t", std::string("neither directives nor models determine allowable ") + var_name + " range" );

            if( var_min >= var_max )
               throw Error::Constructing( "XSec_t", std::string("invalid given range for variable ") + var_name );

         } else {
            if( given_min == -999 )
               given_min = var_min;

            if( given_max == -999 )
               given_max = var_max;

            if( var_min > given_min  ){
               // i.e. function range min is greater than given min
               cerr << "\tWARNING: given min " << var_name << " (" << given_min << ") less than allowed range of models (" << var_min << ")" << endl;
               cerr << "\t\tmin " << var_name << " <- " << var_min << endl;
            } else {
               var_min = given_min;
            };

            if( var_max < given_max  ){
               // i.e. function range min is greater than given max
               cerr << "\tWARNING: given max " << var_name << " (" << given_max << ") less than allowed range of models (" << var_max << ")" << endl;
               cerr << "\t\tmax " << var_name << " <- " << var_max << endl;
            } else {
               var_max = given_max;
            };

            if( var_min >= var_max )
               throw Error::Constructing( "XSec_t", std::string("unable to find valid range for variable ") + var_name );
         };
      };
 
      Enum::pol_state_t XSec_t::GetTargetPolState() const{
         return target_pol;
      };

      Enum::pol_state_t XSec_t::GetBeamPolState() const{
         return beam_pol;
      };

//       final_state_t XSec_t::GetFinalStateType() const{
//          return final_state;
//       };

      LundPID_t XSec_t::Get_Had1_PID() const {
         return INVALID_LUND_PID;
      };

      LundPID_t XSec_t::Get_Had2_PID() const {
         return INVALID_LUND_PID;
      };

      /*
      void XSec_t::DetermineDomain( const sgUtil::ParseInputReturn_t& parsed_input ){
         throw Error::SanityCheckFailure( "XSec_t", "``DetermineDomain'' should have been overwritten by child class" );
      };
      */

      void XSec_t::MakeTerms( const sgUtil::ParseInputReturn_t& parsed_input ){
         throw Error::SanityCheckFailure( "XSec_t", "``MakeTerms'' should have been overwritten by child class" );
      };

   };
};
