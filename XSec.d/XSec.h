/*
  Main class for constructing SIDIS Cross Sections
*/

#ifndef _XSEC_H_
#define _XSEC_H_

#include "Common.d/Enums.h"
#include "Common.d/ParseInput.h"
#include "Common.d/LundPID.h"

#include "XSec_Term.d/XSec_Term.h"

#include "DF.d/Full_DF_Set.h"
#include "FF.d/Full_FF_Set.h"

#include <vector>
#include <cmath>

// predeclare namespace for friend classes
namespace Opt_D1_Spec_Ia {
   class Opt_Thread_t;
};

namespace TMDGen {
   namespace XSec {

      class XSec_t {
         friend class Opt_D1_Spec_Ia::Opt_Thread_t;

         // no eq operator
         NO_EQ_OP( XSec_t );

      protected:

         XSec_t( Var_t& min_in, Var_t& max_in ) : overall_const_factor(0), min(min_in), max(max_in)
            { /* */ };

         // final_state_t final_state;
         Enum::pol_state_t target_pol, beam_pol;

         // Twist and PIDs
         int twist;

         // Terms in the cross section
         std::vector< XSec_Term::XSec_Term_t* > xsec_term;
         XSec_Term::XSec_Term_t* angular_integrated_xsec_term;

         // to hold models
         DF::Full_DF_Set_t *DF_Set;
         FF::Full_FF_Set_t *FF_Set;

         // precomputed y-functions
         bool yFunc_used[5];          // flag if need to compute particular y-function
         double yFunc_array[5];       // in order A, B, C, V, W

         // save angular integrated value
         double ang_int_xsec;

         // precomputed polarization factors (combination of quark charge, beam pol. and target pol.)
         double pol_UU, pol_UT;       // todo: add the rest

         // Determine valid domain
         virtual void DetermineDomain( const sgUtil::ParseInputReturn_t& parsed_input ) = 0;
         void GetRange( const sgUtil::ParseInputReturn_t& parsed_input, std::string var_name, double& var_min, double& var_max );

         // Computing dependent variables
         virtual int ComputeOtherVar( Var_t& ) = 0;

         // overall factors
         const double overall_const_factor;

         // references to variable ranges
         Var_t &min, &max;

         // for overall kinematic and constant factors for the cross section
         virtual double Overall_Kinematic_Factor( const Var_t& var ) const = 0;

         // make needed terms
         virtual void MakeTerms( const sgUtil::ParseInputReturn_t& parsed_input );

         XSec_t( const sgUtil::ParseInputReturn_t& input, Var_t& min_in, Var_t& max_in );
         void FinishConstructing( const sgUtil::ParseInputReturn_t& parsed_input );

         double Eval_Inner( Var_t& var );

         int skip_overall_factors;

      public:
         virtual ~XSec_t();

         double operator() ( Var_t& var );
         double Eval_only_Angular( Var_t& var, double& ang_int_xsec_out  );
         void Eval( Var_t& var );

         Enum::pol_state_t GetTargetPolState() const;
         Enum::pol_state_t GetBeamPolState() const;
         // final_state_t GetFinalStateType() const;

         // copies domain from this class into the passed structs
         // virtual void CopyDomain( Var_t& min_, Var_t& max_ ) const;

         virtual LundPID_t Get_Had1_PID() const;
         virtual LundPID_t Get_Had2_PID() const;
      };
   };
};


#endif
