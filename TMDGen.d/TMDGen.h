/*
   The main public interface for TMDGen, including the tasks of
   - parsing input file
   - instatiating XSec & variable throwers
   - generating events
   - weighting pre-generated events
   - etc.

   Is a virtual class, with children for various processes
   - SIDIS 1had
   - SIDIS dihadron
   - etc.

*/

#ifndef _GMC_TRANS_H_
#define _GMC_TRANS_H_

#include "Common.d/ParseInput.h"
#include "Common.d/Flavor.h"
#include "Common.d/Var.h"

#include "Decay.d/pi0_.h"

#include "RNG.d/GSL_RNG.h"

#include "Thrower.d/VarThrower.h"
#include "Thrower.d/PolThrower.h"

#include "XSec.d/XSec.h"

#ifdef USE_ROOT
#include <TFile.h>
#include <TTree.h>
#endif

#include <string>

// predeclare namespace for friend classes
namespace Opt_D1_Spec_Ia {
   class Opt_Thread_t;
};

namespace TMDGen {

   class TMDGen_t {
      friend class Opt_D1_Spec_Ia::Opt_Thread_t;

   protected:
      Var_t min, max, var;
      std::string N_hadrons;
      int run, itrials, check_Box_acceptance;

      // constant kinematic variables
      double max_weight_factor, min_XS, Target_P_T, Target_P_L, Beam_P_L, Beam_P_T;
      const double E_beam_lab;
      const int beam_charge;

      // beam direction and vertex position (i.e. offset)
      const double vertexX, vertexY, vertexZ, beamTheta, beamPhi;

      // for integrating the cross section per flavor
      int N_integration_warmup, N_integration_calls, N_max_integration_calls;
      double integrated_CX[ GMC_TRANS_N_FLAVORS ];
      double integrated_CX_abserr[ GMC_TRANS_N_FLAVORS ];

      // GSL random number generator
      GSL_RNG_t *r;
      Thrower::VarThrower_t *var_thrower;
      Thrower::PolThrower_t *polarization_thrower;

      // Cross section
      XSec::XSec_t *xsec;

      // common constructor
      void Construct( sgUtil::ParseInputReturn_t input );

      // Determine maximum value of the cross section
      int max_weight_iters;
      int DetermineMaxWeight( unsigned int N );

      Decay::pi0 *decay_pi0;

      // for output to ROOT
      int output_to_root, root_init;
#ifdef USE_ROOT
      TFile* rootfile;
      TTree* tree;
#endif
      virtual int Init_Root(const std::string& filename ) = 0;
      virtual int Output_to_Root();

      // Generating with or without weights
      int GenEvent_w_Weight();
      int GenEvent_wo_Weight();
      int (TMDGen_t::*genevent_ptr)();

      // to ensure correct process and final state for given child
      virtual void Construct_Child( const sgUtil::ParseInputReturn_t& input ) = 0;

      // for specific details of reconstructing event
      // which depend on specifics of children classes
      virtual int ReconstructEvent_details() = 0;

      // for radiative effects
      bool radcor_on;


      // only children can call this
      // redundant--as class already has purely virtual functions
      TMDGen_t();

   public:
      virtual ~TMDGen_t();

      // Initialize (includes precomputation of integrals, if needed)
      int Initialize();

      // Generate an event
      int GenerateEvent();

      // determine weight for given event
      void WeightEvent( Var_t& var );

      // fill in all event variables, as the other
      // procedure fills in just minimal set of values
      int ReconstructEvent();

      // write out results to file
      int WriteOutEvent();

      // Get codes of what formats are set to output
      int OutputCodes() const;

      // get pointer to variables
      const Var_t* GetEventVariables();

      // Copy min max from local variables
      virtual void CopyMinMax( Var_t& min, Var_t& max );

   };
};


#endif
