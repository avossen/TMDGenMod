
#include "TMDGen.d/SIDIS_1had.h"

#include "Common.d/ParseInput.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/Consts.h"
#include "Common.d/LundPID.h"
#include "Common.d/Var.h"

#include "Thrower.d/BasicThrower_1had.h"
//#include "Thrower.d/PWC_Thrower_1had.h"

#include "XSec.d/XSec.h"
#include "XSec.d/SIDIS_1had.h"

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

#include <iostream>
using std::endl;
using std::cerr;

namespace TMDGen {

   int SIDIS_1had_t::Init_Root( const std::string& filename ){
#ifndef USE_ROOT
      cerr << "\t\tERROR: Root not available.  Was not compiled with '-DUSE_ROOT' option" << endl;
      return 127;
#else
      int ierr = 1;
      rootfile = new TFile (filename.data(), "RECREATE", "TMDGen Output");

      if( !rootfile->IsZombie() ){
         ierr = 0;
         tree = new TTree("tree", "event tree");

         int iflavor;

         tree->Branch("XS",            &var.XS, "XS/D" );
         tree->Branch("Thrown_Weight", &var.thrown_weight, "Thrown_Weight/D" );
         tree->Branch("Weight",        &var.weight, "Weight/D" );

         tree->Branch("XTrue", &var.x, "XTrue/D" );
         tree->Branch("YTrue", &var.y, "YTrue/D" );
         tree->Branch("ZTrue", &var.z, "ZTrue/D" );
         tree->Branch("pT",    &var.pT, "pT/D" );
         tree->Branch("kT",    &var.kT, "kT/D" );
         tree->Branch("P_hperp_True",   &var.P_hperp,   "P_hperp_True/D" );

         tree->Branch("Q2True", &var.Q2, "Q2True/D" );
         tree->Branch("W2True", &var.W2, "W2True/D" );
         tree->Branch("quark_flavor", &iflavor, "quark_flavor/I" );

         tree->Branch("phi_h_True", &var.phi_h, "phi_h_True/D" );
         tree->Branch("phi_S_True", &var.phi_S, "phi_S_True/D" );

         tree->Branch("e2_rP",     &var.lep_3.E, "e2_rP/D" );
         tree->Branch("e2_rTheta", &var.lep_3.theta, "e2_rTheta/D" );
         tree->Branch("e2_rPhi",   &var.lep_3.phi, "e2_rPhi/D" );

         tree->Branch("h1_rP",     &var.had.P, "h1_rP/D" );
         tree->Branch("h1_rTheta", &var.had.theta, "h1_rTheta/D" );
         tree->Branch("h1_rPhi",   &var.had.phi, "h1_rPhi/D" );

      };

      return ierr;
#endif
   };

   void SIDIS_1had_t::Construct_Child( const sgUtil::ParseInputReturn_t& parsed_input ) {

      // check process
      sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find("Process");

      if( it == parsed_input.end() ){
         throw Error::Constructing( "SIDIS_1had_t", "no 'Process' directive given" );
      };

      if( it->second != "SIDIS" ){
         throw Error::Constructing( "SIDIS_1had_t", std::string("invalid process '") +
                                    it->second + "' for given class.  Expected 'SIDIS'." );
      };

      // check final state

      it = parsed_input.find("Final_State");
      if( it == parsed_input.end() )
         throw Error::Constructing( "XSec_t", "no 'Final_State' directive given" );
      if( it->second != "Single Hadron" ){
         throw Error::Constructing( "SIDIS_1had_t", std::string("invalid final state '") +
                                    it->second + "' for given class.  Expected 'Single Hadron'" );
      };

      // allocate cross section
      try{
         xsec = new XSec::SIDIS_1had_t( parsed_input, min, max );
      }
      catch( std::exception& e ){
         throw Error::Constructing( "TMDGen::SIDIS_1had_t", std::string("error in required class XSec::SIDIS_1had_t:\n\t\t") + e.what() );
      };

      // construct parent
      Construct( parsed_input );

      // check range
      if( max.P_hperp <= 0 ){
         delete xsec;
         delete r;
         throw Error::Constructing( "TMDGen_t", "non positive maximum value of P_hperp" );
      };

      // allocate variable thrower
      it = parsed_input.find("VarThrower");
      if( it == parsed_input.end() ){
         throw Error::Constructing("TMDGen_t", "No 'VarThrower' directive");
      } else if ( it == parsed_input.end() || it->second == "Basic" ){
         int N_warmup =    10000;
         int N_calls =     1000000;
         int N_max_calls = 10000000;

         it = parsed_input.find("Basic_VarThrower_N_warmup");
         if( it != parsed_input.end() )
            N_warmup = atoi( it->second.data() );

         it = parsed_input.find("Basic_VarThrower_N_calls");
         if( it != parsed_input.end() )
            N_calls = atoi( it->second.data() );

         it = parsed_input.find("Basic_VarThrower_N_max_calls");
         if( it != parsed_input.end() )
            N_max_calls = atoi( it->second.data() );

         var_thrower = new Thrower::BasicThrower_1had_t( r, N_warmup, N_calls, N_max_calls, min, max );

         /*
      } else if ( it->second == "PWC" ){
         int N_evals_per_bin = 100000;

         Var_t nper;
         nper.x = 5;
         nper.y = 1;
         nper.z = 5;
         nper.P_hperp = 1;
         nper.phi_h = 1;
         nper.pT = 1;
         nper.phi_pT = 1;

         it = parsed_input.find("PWC_VarThrower_N_evals_per_bin");
         if( it != parsed_input.end() )
            N_evals_per_bin = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_x");
         if( it != parsed_input.end() )
            nper.x = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_y");
         if( it != parsed_input.end() )
            nper.y = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_z");
         if( it != parsed_input.end() )
            nper.z = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_P_hperp");
         if( it != parsed_input.end() )
            nper.P_hperp = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_phi_h");
         if( it != parsed_input.end() )
            nper.phi_h = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_pT");
         if( it != parsed_input.end() )
            nper.pT = atoi( it->second.data() );

         it = parsed_input.find("PWC_VarThrower_N_phi_pT");
         if( it != parsed_input.end() )
            nper.phi_pT = atoi( it->second.data() );

         var_thrower = new Thrower::PWC_Thrower_1had_t( r, min, max, N_evals_per_bin, nper );
         */
      } else {
         throw Error::Constructing( "TMDGen::SIDIS_1had", std::string("Invalid 'VarThrower' directive: '") +
                                    it->second + "'" );
      };

      cerr << endl << "\tVariable Ranges" << endl;
      cerr << "\t\tQ2 in " << min.Q2 << " -> " << max.Q2 << endl;
      cerr << "\t\tW2 in " << min.W2 << " -> " << max.W2 << endl;
      cerr << "\t\tx in " << min.x << " -> " << max.x << endl;
      cerr << "\t\ty in " << min.y << " -> " << max.y << endl;
      cerr << "\t\tz in " << min.z << " -> " << max.z << endl;
      cerr << "\t\tpT in " << min.pT << " -> " << max.pT << endl;
      cerr << "\t\tkT in " << min.kT << " -> " << max.kT << endl;
      cerr << "\t\tP_hperp in " << min.P_hperp << " -> " << max.P_hperp << endl;
   };


   // NOW PUBLIC MEMBER FUNCTIONS


   SIDIS_1had_t::SIDIS_1had_t( const char* parameter_file_name ) : TMDGen_t() {
      sgUtil::ParseInputReturn_t parsed_input;

      if( sgUtil::ParseInput( parameter_file_name, parsed_input ) ){
         throw Error::Constructing( "TMDGen_t", "problem with input file" );
      };

      Construct_Child( parsed_input );
   };

   SIDIS_1had_t::SIDIS_1had_t( const sgUtil::ParseInputReturn_t& parsed_input ) {
      Construct_Child( parsed_input );
   };

   SIDIS_1had_t::~SIDIS_1had_t(){
      if( xsec )
         delete xsec;
      if( var_thrower )
         delete var_thrower;
   };

   // fill in all event variables, as the other
   // procedure fills in just minimal set of values
   int SIDIS_1had_t::ReconstructEvent_details(){
      int ierr = var.ReconstructEvent_SIDIS_1had( xsec->GetTargetPolState() );
      if( check_Box_acceptance ){
         if( !ierr )
            ierr = !var.lep_2.InBoxAcc();
         if( !ierr )
            ierr = !var.had.InBoxAcc();
      };
      return ierr;
   };

};
