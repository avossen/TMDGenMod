
#include "TMDGen.d/TMDGen.h"

#include "Common.d/ParseInput.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/Consts.h"
#include "Common.d/LundPID.h"

#include "Thrower.d/PolThrower_UU.h"
#include "Thrower.d/PolThrower_UT.h"

#include "RNG.d/GSL_RNG.h"

#include "XSec.d/XSec.h"

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include <cmath>
#include <sstream>
#include <iostream>
using std::cerr;
using std::cout;
using std::endl;

#include <cstdlib>
 
namespace TMDGen {

   // all instantiation done through ``Construct'' member function
   TMDGen_t::TMDGen_t() : min("TMDGen_t::min"), max("TMDGen_t::max"), var("TMDGen_t::var"),
                                min_XS(0), E_beam_lab(27.57), beam_charge(1), vertexX(0), vertexY(0), vertexZ(0), beamTheta(0), beamPhi(0),
                                xsec(0), decay_pi0(0) { /* */ };

   TMDGen_t::~TMDGen_t(){
      delete polarization_thrower;
      delete r;

#ifdef USE_ROOT
      if( rootfile ){
         rootfile->Write();
         rootfile->Close();
         rootfile->Delete(); // deletes everything in the TFile, including the TTree
         delete rootfile;
      };
#endif

   };

   void TMDGen_t::Construct( sgUtil::ParseInputReturn_t parsed_input ) {

      sgUtil::ParseInputReturn_iterator_t it;

      // check of process already done in children classes

      // xsec allocated in child class
      if( !xsec ){
         throw Error::SanityCheckFailure( "TMDGen_t", "child class did not set cross section pointer" );
      };


      // set to default values
      Target_P_T = 1;
      Target_P_L = 1;

      // change, if available
      it = parsed_input.find("Target_P_T");
      if( it != parsed_input.end() )
         Target_P_T = atof( it->second.data() );
      if( fabs(Target_P_T) > 1 ){
         delete xsec;
         throw Error::Constructing("TMDGen_t", "Magnitude of transverse target polarization > 1");
      };

      it = parsed_input.find("Target_P_L");
      if( it != parsed_input.end() )
         Target_P_L = atof( it->second.data() );
      if( fabs(Target_P_L) > 1 ){
         delete xsec;
         throw Error::Constructing("TMDGen_t", "Magnitude of longitudinal target polarization > 1");
      };

      check_Box_acceptance = 0;
      it = parsed_input.find("Apply_Box_Acceptance");
      if( it != parsed_input.end() ){
         if( it->second.empty() || it->second == "1" || it->second[0] == 'T' || it->second[0] == 't' || 
             it->second[0] == 'Y' || it->second[0] == 'y' )
            check_Box_acceptance = 1;
         cerr << "\t\tApply_Box_Acceptance = " << check_Box_acceptance << endl;
      };




      // run number
      it = parsed_input.find("Run");
      if( it != parsed_input.end() )
         run = atoi( it->second.data() );

      // Alloc RNG
      try{
         r = new GSL_RNG_t( gsl_rng_ranlxs2 );
      }
      catch( std::exception& e ){
         delete xsec;
         throw Error::Constructing( "TMDGen_t", std::string("error allocating RNG:\n\t\t") + e.what() );
      };

      switch( xsec->GetBeamPolState() ){
      case Enum::UNPOL:
         Beam_P_T = 0;
         Beam_P_L = 0;

         switch( xsec->GetTargetPolState() ){
         case Enum::UNPOL:
            polarization_thrower = new Thrower::PolThrower_UU_t;
            Target_P_T = 0;
            Target_P_L = 0;
            break;
         case Enum::TRANS:
            polarization_thrower = new Thrower::PolThrower_UT_t( r, Target_P_T );
            Target_P_L = 0;
            break;
         case Enum::LONG:
            Target_P_T = 0;
            delete xsec;
            delete r;
            throw Error::Constructing("TMDGen_t", "Not yet programmed unpolarized beam, longitudinally polarized target");
            break;
         };
         break;
      case Enum::TRANS:
         Beam_P_T = 0;
         delete r;
         delete xsec;
         throw Error::Constructing("TMDGen_t", "Not yet programmed polarized beam");
         break;
      case Enum::LONG:
         Beam_P_T = 0;
         delete r;
         delete xsec;
         throw Error::Constructing("TMDGen_t", "Not yet programmed polarized beam");
         break;
      };

      // beam energy
      it = parsed_input.find("E_beam_lab");
      if( it != parsed_input.end() ){
         double temp_E = atof(it->second.data());
         if( temp_E > 0 && temp_E < 1e6 ){
            *const_cast< double* >( &E_beam_lab ) = temp_E;
         } else {
            cerr << "WARNING in TMDGen_t::TMDGen_t(...): Invalid 'E_beam_lab' directive '" << it->second << "'.  Using default value of " << E_beam_lab << "." << endl;
         };
      };
      var.lep_0.E = E_beam_lab;
      var.lep_0.P = E_beam_lab;
      var.lep_0.M = ELECTRON_MASS;
      var.lep_0.theta = 0;
      var.lep_0.phi = 0;
      var.lep_1 = var.lep_0;   // initialize to some values
      //cerr << "E_beam_lab = " << E_beam_lab << ' ' << &var.lep_0.E << ' ' << &var.lep_1.E << ' ' << &var << endl;

      // beam charge
      it = parsed_input.find("Beam_Charge");
      if( it != parsed_input.end() ){
         int charge = atoi(it->second.data());
         if( charge != 1 && charge != -1 ){
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t", std::string("Invalid beam charge: '") + it->second + "'" );
         };

         *const_cast< int* >( &beam_charge ) = charge;
      };

      // vertex position and beam direction at the vertex
      it = parsed_input.find("Vertex_Offset_X");
      if( it != parsed_input.end() ){
         double vertexX_ = atof(it->second.data());
         if( vertexX_ > -3 && vertexX_ < 3 ){  // to avoid crazy values
            *const_cast< double* >( &vertexX ) = vertexX_;
         } else {
            cout << "WARNING: vertexX position at" << vertexX_ << " is out of range." << endl;
         };
      };

      // vertex position and beam direction at the vertex
      it = parsed_input.find("Vertex_Offset_Y");
      if( it != parsed_input.end() ){
         double vertexY_ = atof(it->second.data());
         if( vertexY_ > -3 && vertexY_ < 3 ){  // to avoid crazy values
            *const_cast< double* >( &vertexY ) = vertexY_;
         } else {
            cout << "WARNING: vertexY position at" << vertexY_ << " is out of range." << endl;
         };
      };

      // vertex position and beam direction at the vertex
      it = parsed_input.find("Vertex_Offset_Z");
      if( it != parsed_input.end() ){
         double vertexZ_ = atof(it->second.data());
         if( vertexZ_ > -300 && vertexZ_ < 300 ){  // to avoid crazy values
            *const_cast< double* >( &vertexZ ) = vertexZ_;
         } else {
            cout << "WARNING: vertexZ position at" << vertexZ_ << " is out of range." << endl;
         };
      };

      // vertex position and beam direction at the vertex
      it = parsed_input.find("Beam_Theta_at_Vertex");
      if( it != parsed_input.end() ){
         double theta = atof(it->second.data());
         if( theta >= 0 && theta < 3.1416 ){  // to avoid crazy values
            *const_cast< double* >( &beamTheta ) = theta;
         } else {
            cout << "WARNING: beam theta at vertex (" << theta << ") is not in range." << endl;
         };
      };

      // vertex position and beam direction at the vertex
      it = parsed_input.find("Beam_Phi_at_Vertex");
      if( it != parsed_input.end() ){
         double phi = atof(it->second.data());
         if( phi >= 0 && phi < 6.29 ){  // to avoid crazy values
            *const_cast< double* >( &beamPhi ) = phi;
         } else {
            cout << "WARNING: beam phi at vertex (" << phi << ") is not in range." << endl;
         };
      };

      // to determine max weight
      max_weight_iters = -999;
      it = parsed_input.find("MaxWeight");
      if( it == parsed_input.end() ){
         max_weight_iters = 100000;
      } else if (it->second.substr(0,3) == "Est"){
         max_weight_iters = 100000;
         if( it->second.size() > 4 )
            max_weight_iters = atoi( it->second.substr(3).data() );

         if( max_weight_iters < 1 ){
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t", 
                                      std::string("Invalid number of attempts for estimating the max. weight: '") +
                                      it->second.substr(3) + "'" );
         };
      } else if (it->second.substr(0,3) == "Set"){
         if( it->second.size() < 5 ){
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t", 
                                      std::string("Invalid given max. weight: '") +
                                      it->second.substr(3) + "'" );
         };
         max_weight_factor = atof( it->second.substr(3).data() );
         if( max_weight_factor <= 0 ){
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t", 
                                      std::string("Invalid given max. weight: '") +
                                      it->second.substr(3) + "'" );
         };
         max_weight_factor = 1./max_weight_factor;

         cerr << "\tSetting maximum weight before factor to " << 1./max_weight_factor << endl;
      } else if (it->second == "Skip"){
	max_weight_iters = 0;
      } else {
         delete r;
         delete xsec;
         throw Error::Constructing("TMDGen_t", std::string("invalid option for 'MaxWeight' directive: '")
                                   + it->second + "'" );
      };

      // output options
      output_to_root = 0;
#ifdef USE_ROOT
      rootfile = 0;
      it = parsed_input.find("Output_To_Root_File");
      if( it != parsed_input.end() && !it->second.empty() ){
         if( Init_Root(it->second) ){  // it->second is filename
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t", "cannot create root file");
         };

         output_to_root = 1;
      };
#endif

      genevent_ptr = &TMDGen_t::GenEvent_wo_Weight;
      it = parsed_input.find("Generate_with_Weights");
      if( it != parsed_input.end() ){
         if( it->second.empty() || it->second[0] == 'y' || it->second[0] == 'Y' ){
            genevent_ptr = &TMDGen_t::GenEvent_w_Weight;
         } else if ( it->second[0] != 'n' && it->second[0] != 'N' ){
            delete r;
            delete xsec;

            throw Error::Constructing("TMDGen_t",
                                      std::string("invalid option for 'Generate_with_Weights' directive: '")
                                      + it->second + "'.  Must be start with 'Y','y','N', or 'n'" );
         };
      };


      // determine if radiative effects are turned on
      // defaults to off
      it = parsed_input.find("RadCor");
      radcor_on = 0;

      if( it != parsed_input.end() ){
         if( it->second == "TRUE" || it->second == "true" || it->second == "True" ){
            radcor_on = 1;
            throw Error::Constructing("TMDGen_t", "not yet programmed Radiative Corrections" );
         } else if ( it->second == "FALSE" || it->second == "false" || it->second == "False" ){
            radcor_on = 0;
         } else {
            delete r;
            delete xsec;
            throw Error::Constructing("TMDGen_t",
                                      "Invalid 'RadCor' directive: must be 'TRUE', 'True', 'true', 'FALSE', 'False', or 'false'" );
         };
      };

      if( radcor_on ){
         // insert code to construct child of RadCor_t class
         // RadCor_t* radcor = new ???

         // insert code to do any initializations neccessary
         // if( radcor->Init( parsed_input ) ){
         // throw error
         //};
      };
   };


   int TMDGen_t::Initialize(){
      // Prepare variable thrower
      int ierr = var_thrower->Initialize( var, *xsec );

      if( !ierr )
         if( max_weight_iters > 0 )          // determine max weight, if needed
            ierr = DetermineMaxWeight( max_weight_iters );

      var.integrating = 0;

      return ierr;
   };


   // generate an event with a weight
   int TMDGen_t::GenEvent_w_Weight(){
      double thrown_pdf_val = 1;
      int ierr = 0;

      var.lep_0.E = E_beam_lab;
      var.lep_0.P = E_beam_lab;
      var.lep_0.M = ELECTRON_MASS;
      var.XS = 0;
      var.weight = 0;

      // Dec. 2013: set the direction based on input parameters
      var.vertex[0] = vertexX;
      var.vertex[1] = vertexY;
      var.vertex[2] = vertexZ;
      var.lep_0.theta = beamTheta;
      var.lep_0.phi = beamPhi;

      if( radcor_on ){
         // call initial state radiation
         // if( radcor->Cor_InitialState( var.ISR, lep_0, lep_1 ) ){
         // throw error
         //};
         // If radcor model doesn't involve separate initial state radiation
         // then it also just copies var.lep_1 = var.lep_0
      } else { 
         // let the kinematics of the lepton in the vertex equal those for the lepton in the beam
         var.ISR = 0;
         var.lep_1 = var.lep_0;
      };

      if( !ierr ){
         int N_max_tries = 2;   // todo: allow this to be changed by instr. file
         int i = 0;

         do { 
            // set polarization to be +/- for beam and/or target 
            polarization_thrower->Throw( var );

            // determine variables
            var_thrower->Throw( var, thrown_pdf_val );

            // evaluate cross section
            xsec->Eval( var );

            ++i;
         } while( fabs(var.XS) <= min_XS && i < N_max_tries );


         if( i == N_max_tries && fabs(var.XS) < min_XS ){
            std::stringstream ss;
            ss << "Absolute value of the cross section less than " << min_XS << " for " << N_max_tries << " attempts.  Flavor is " << flavor_string[var.flavor];
            throw Error::Base_t( "Error in", "TMDGen_t::GenEvent_w_Weight(...)", ss.str() );
         };

         ierr = 1;
         if( fabs(var.XS) > min_XS ){
            var.weight = var.XS / thrown_pdf_val * max_weight_factor;
            //std::cout << "b " << var.XS << ' ' << thrown_pdf_val << ' ' << max_weight_factor << " = " << var.weight <<  std::endl;
            var.thrown_weight = 1;
            itrials = i;
            ierr = 0;
         };
      };

      //cerr << "\tGenerating with weights done..." << endl;

      return ierr;
   };


   // generate an event without a weight
   int TMDGen_t::GenEvent_wo_Weight(){

      int ierr = 0;
      do {
         do { /* */ } while ( GenEvent_w_Weight() );
      } while ( fabs(var.weight) < r->EvalUnif() );

      var.thrown_weight = var.weight;
      var.weight = (var.weight) > 0 ? 1 : -1;
      //std::cout << var.thrown_weight << ' ' << ierr << endl;

      return ierr;
   };

   // Determine maximum value of the cross section
   int TMDGen_t::DetermineMaxWeight( unsigned int N ){

      cerr << "\tDetermining maximum weight" << endl;

      // decide somehow if need to turn off radcor while determining max weight

      // generate few events for the first guess at the max
      //std::cout << "! Determine max weight" << endl;
      max_weight_factor = 1;
      double max_weight = 0;
      unsigned int ierr = 0;
      int n_non_zero = 0;
      for( int i=0; i<100000 && n_non_zero<10; ++i ){
         ierr = GenEvent_w_Weight();

         if( var.weight <= 0 && !ierr )
            ierr = 1;
      
         if( !ierr && check_Box_acceptance ){
            ierr = ReconstructEvent();
            if( ierr )
               var.weight = 0;
         };

         if( var.weight > max_weight )
            max_weight = var.weight;

         if( var.weight > 0 )
            ++n_non_zero;

         //std::cout << "c " << var.weight << ' ' << i << ' ' << n_non_zero << endl;
      };

      ierr = 0;
      if( max_weight <= 0 ){
         cerr << "\tERROR: Max_Weight <= 0" << endl;
         ierr = 1;
      };

      max_weight_factor = ( ierr ? -1 : 1./max_weight );

      if( !ierr ){
         for( unsigned int i=0; i<N; ++i ){
            ierr += GenEvent_w_Weight();

            if( check_Box_acceptance ){
               int ierr2 = ReconstructEvent();
               if( ierr2 )
                  var.weight = 0;
            };
            //std::cout << "c " << var.weight << endl;

            if( var.weight > 1 )
               max_weight_factor /= var.weight;
         };

         if( ierr < N )
            ierr = 0;
      };

      if( !ierr ){
         cerr << "\t\tMax weight before factor is " << 1./max_weight_factor << endl;
      } else {
         cerr << "\t\tError determining max. weight" << endl;
      };

      //std::cout << "! Determine max weight DONE" << endl;

      double min_XS = 0; //1e-18 / max_weight_factor;
      std::cerr << "\t\tMin cross section is " << min_XS << std::endl;

      //cerr << "\t\tDetermine max weight returning ierr = " << ierr << endl;
      return ierr;
   };

   int TMDGen_t::WriteOutEvent(){
      int ierr = 0;
#ifdef USE_ROOT
      if( output_to_root ){
         if( Output_to_Root() ){
            cerr << "\t\tERROR writting out to ROOT" << endl;
            ierr += 1;
         };
      };
#endif
#if !( defined(USE_ROOT) )
      ierr = 1;
#endif

      return ierr;
   };

   int TMDGen_t::GenerateEvent(){
      //std::cout << "\tGenerating next event..." << endl;
      return (this->*genevent_ptr)();
   };

   int TMDGen_t::OutputCodes() const {
      return output_to_root*2;
   };

   int TMDGen_t::ReconstructEvent(){
      int ierr = ReconstructEvent_details();

      if( !ierr ){

         var.FSR = 0;
         if( radcor_on ){
            // call final state radiation
            // if( radcor->Cor_FinalState( var.FSR, lep_2, lep_3 ) ){
            // throw error
            //};
            // If radcor model doesn't involve separate final state radiation
            // then it also just copies var.lep_3 = var.lep_2
         } else {
            // let kinematics of final lepton equal lepton coming out of primary vertex
            var.FSR = 0;
            var.lep_3 = var.lep_2;
         };
      };

      return ierr;
   };

   int TMDGen_t::Output_to_Root(){
#ifndef USE_ROOT
      cerr << "\t\tERROR: Root not available.  Was not compiled with '-DUSE_ROOT' option" << endl;
      return 127;
#else
      int ierr = 1;
      if( rootfile ){
         // convert enum to int
         int iflavor = var.flavor;
         tree->SetBranchAddress("quark_flavor", &iflavor );

         // fill tree
         tree->Fill();
         ierr = 0;
      };
      return ierr;
#endif
   };

   // Copy min max from local variables
   void TMDGen_t::CopyMinMax( Var_t& min_out, Var_t& max_out ){
      min_out.SmallCopy( min );
      max_out.SmallCopy( max );
   };


   // determine weight for given event
   // assume all variables are accurate in the event
   void TMDGen_t::WeightEvent( Var_t& var_in ){

      double thrown_pdf_val = 1.;

      var_thrower->ThrowFlavor( var_in, thrown_pdf_val );
      var_thrower->Throw_pT( var_in, thrown_pdf_val );

      // evaluate cross section, sets XS
      xsec->Eval( var_in );

      //std::cout << var_in.x << ' ' << var_in.y << ' ' << var_in.z << ' ' << var_in.P_hperp << ' ' << var_in.had_0.M << std::endl;

      //cout << "g " << var_in.XS << ' ' << var_in.XS_int << ' ' << thrown_pdf_val << endl;

      var_in.weight = var_in.XS / thrown_pdf_val;
      var_in.weight_int = var_in.XS_int / thrown_pdf_val;
   };


};
