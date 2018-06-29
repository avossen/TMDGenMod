#include "TMDGen.d/SIDIS_2had.h"

#include "Common.d/Consts.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/LundPID.h"
#include "Common.d/ParseInput.h"
#include "Common.d/Var.h"

#include "Decay.d/pi0_ROOT.h"

#include "Thrower.d/BasicThrower_2had.h"

#include "XSec.d/XSec.h"
#include "XSec.d/SIDIS_2had.h"

#include <iostream>
using std::endl;
using std::cerr;
using std::cout;

#ifdef USE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

namespace TMDGen {

   int SIDIS_2had_t::Init_Root( const std::string& filename ){
#ifndef USE_ROOT
      cerr << "\t\tERROR: Root not available.  Was not compiled with '-DUSE_ROOT' option" << endl;
      return 127;
#else
      int ierr = 1;
      rootfile = new TFile (filename.data(), "RECREATE", "TMDGen Output");

      if( !rootfile->IsZombie() ){
         ierr = 0;
         tree = new TTree("tree", "SIDIS 2had event tree");

         int iflavor;

         tree->Branch("XS", &var.XS, "XS/D" );
         tree->Branch("ISR", &var.ISR, "ISR/I" );
         tree->Branch("FSR", &var.FSR, "FSR/I" );
         tree->Branch("Thrown_Weight", &var.thrown_weight, "Thrown_Weight/D" );
         tree->Branch("Weight", &var.weight, "Weight/D" );
         tree->Branch("P_T", &var.P_T, "P_T/D" );
         tree->Branch("P_L", &var.P_L, "P_L/D" );

         tree->Branch("XTrue", &var.x, "XTrue/D" );
         tree->Branch("YTrue", &var.y, "YTrue/D" );
         tree->Branch("ZTrue", &var.z, "ZTrue/D" );
         tree->Branch("pT",    &var.pT, "pT/D" );
         tree->Branch("phi_pT",&var.phi_pT, "phi_pT/D" );
         tree->Branch("kT",    &var.kT, "kT/D" );
         tree->Branch("phi_kT",&var.phi_kT, "phi_kT/D" );
         tree->Branch("P_hperp_True", &var.P_hperp, "P_hperp_True/D" );
         tree->Branch("M_hh_True",    &var.had_0.M, "M_hh_True/D" );

         tree->Branch("Q2True", &var.Q2, "Q2True/D" );
         tree->Branch("W2True", &var.W2, "W2True/D" );
         tree->Branch("quark_flavor", &iflavor, "quark_flavor/I" );

         tree->Branch("phi_h_True", &var.phi_h, "phi_h_True/D" );
         tree->Branch("phi_S_True", &var.phi_S, "phi_S_True/D" );
         tree->Branch("phi_R_True", &var.phi_R, "phi_R_True/D" );
         tree->Branch("cos_vartheta_True", &var.cos_vartheta, "cos_vartheta_True/D" );
         tree->Branch("psi_True",   &var.psi,    "psi_True/D" );

         tree->Branch("e2_iLType", &e2_iLType,       "e2_iLType/I" );
         tree->Branch("e2_rP",     &var.lep_3.E,     "e2_rP/D" );
         tree->Branch("e2_rTheta", &var.lep_3.theta, "e2_rTheta/D" );
         tree->Branch("e2_rPhi",   &var.lep_3.phi,   "e2_rPhi/D" );

         // parent hadron
         h1_P = 0;
         h2_P = 0;

         tree->Branch("h0_rP",     &var.had_0.P,     "h0_rP/D" );
         tree->Branch("h0_rTheta", &var.had_0.theta, "h0_rTheta/D" );
         tree->Branch("h0_rPhi",   &var.had_0.phi,   "h0_rPhi/D" );

         tree->Branch("h1_iLType", &h1_iLType,       "h1_iLType/I" );
         tree->Branch("h1_rP",     &h1_P,            "h1_rP/D" );
         tree->Branch("h1_rTheta", &var.had_1.theta, "h1_rTheta/D" );
         tree->Branch("h1_rPhi",   &var.had_1.phi,   "h1_rPhi/D" );

         tree->Branch("h2_iLType", &h2_iLType,       "h2_iLType/I" );
         tree->Branch("h2_rE",     &var.had_2.E,     "h2_E/D" );
         tree->Branch("h2_rP",     &h2_P,            "h2_rP/D" );
         tree->Branch("h2_rTheta", &var.had_2.theta, "h2_rTheta/D" );
         tree->Branch("h2_rPhi",   &var.had_2.phi,   "h2_rPhi/D" );

         if( xsec->Get_Had1_PID() == PI_ZERO || xsec->Get_Had2_PID() == PI_ZERO ){
            tree->Branch("c1_rP",     &var.gamma_1.P,     "c1_rP/D" );
            tree->Branch("c1_rTheta", &var.gamma_1.theta, "c1_rTheta/D" );
            tree->Branch("c1_rPhi",   &var.gamma_1.phi,   "c1_rPhi/D" );

            tree->Branch("c2_rP",     &var.gamma_2.P,     "c2_rP/D" );
            tree->Branch("c2_rTheta", &var.gamma_2.theta, "c2_rTheta/D" );
            tree->Branch("c2_rPhi",   &var.gamma_2.phi,   "c2_rPhi/D" );
         };
      };

      return ierr;
#endif
   };

   int SIDIS_2had_t::Output_to_Root(){
#ifndef USE_ROOT
      cerr << "\t\tERROR: Root not available.  Was not compiled with '-DUSE_ROOT' option" << endl;
      return 127;
#else
      int ierr = 1;
      if( rootfile ){
         // convert enum to int
         int iflavor = var.flavor;
         tree->SetBranchAddress("quark_flavor", &iflavor );

         e2_iLType = -11 * beam_charge;
         h1_iLType = xsec->Get_Had1_PID();
         h2_iLType = xsec->Get_Had2_PID();

         h1_P = var.had_1.P;
         if( xsec->Get_Had1_PID() < 0 )
            h1_P *= -1.;

         h2_P = var.had_2.P;
         if( xsec->Get_Had2_PID() < 0 )
            h2_P *= -1.;

         // fill tree
         tree->Fill();
         ierr = 0;
      };
      return ierr;
#endif
   };

   void SIDIS_2had_t::Construct_Child( const sgUtil::ParseInputReturn_t& parsed_input ) {

      // check process
      sgUtil::ParseInputReturn_const_iterator_t it = parsed_input.find("Process");

      if( it == parsed_input.end() ){
         throw Error::Constructing( "SIDIS_2had_t", "no 'Process' directive given" );
      };

      if( it->second != "SIDIS" ){
         throw Error::Constructing( "SIDIS_2had_t", std::string("invalid process '") +
                                    it->second + "' for given class.  Expected 'SIDIS'." );
      };

      // check final state

      it = parsed_input.find("Final_State");
      if( it == parsed_input.end() )
         throw Error::Constructing( "XSec_t", "no 'Final_State' directive given" );
      if( it->second != "Dihadron" ){
         throw Error::Constructing( "SIDIS_2had_t", std::string("invalid final state '") +
                                    it->second + "' for given class.  Expected 'Dihadron'" );
      };

      // allocate cross section
      try {
         xsec = new XSec::SIDIS_2had_t( parsed_input, min, max );
      }
      catch( std::exception& e ){
         cerr << "Caught error: " << e.what() << endl;
         throw Error::Constructing( "TMDGen::SIDIS_2had_t", std::string("error in required class XSec::SIDIS_2had_t:\n\t\t") + e.what() );
      }
      catch(...){
         cerr << "Error Constructing 'TMDGen::SIDIS_2had_t': unknown error in required class XSec::SIDIS_2had_t" << endl << endl;
         throw;
      };

      // construct parent TMDGen class
      Construct( parsed_input );

      // Var_t only has two gammas, not four, so cannot have two pi0's
      // would also have to update box acceptance
      if( xsec->Get_Had1_PID() == PI_ZERO && xsec->Get_Had2_PID() == PI_ZERO )
         throw Error::Constructing( "TMDGen::SIDIS_2had_t", "Not yet programmed: pi0 pi0 pair");

      // check if need decay
      if( xsec->Get_Had1_PID() == PI_ZERO || xsec->Get_Had2_PID() == PI_ZERO ){
#ifdef USE_ROOT
         decay_pi0 = new Decay::pi0_ROOT();
#else
         throw Error::Constructing( "TMDGen::SIDIS_2had_t", "Not yet programmed pi0 decay without ROOT");
#endif
      };

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

         var_thrower = new Thrower::BasicThrower_2had_t( r, N_warmup, N_calls, N_max_calls, min, max );
      } else {
         throw Error::Constructing( "TMDGen::SIDIS_2had", std::string("Invalid 'VarThrower' directive: '") +
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
      cerr << "\t\tM_hh in " << min.had_0.M << " -> " << max.had_0.M << endl;
   };


   // NOW PUBLIC MEMBER FUNCTIONS


   SIDIS_2had_t::SIDIS_2had_t( const char* parameter_file_name ) : TMDGen_t() {
      sgUtil::ParseInputReturn_t parsed_input;

      if( sgUtil::ParseInput( parameter_file_name, parsed_input ) ){
         throw Error::Constructing( "TMDGen_t", "problem with input file" );
      };

      Construct_Child( parsed_input );
   };

   SIDIS_2had_t::SIDIS_2had_t( const sgUtil::ParseInputReturn_t& parsed_input ) : TMDGen_t() {
      Construct_Child( parsed_input );
   };

   SIDIS_2had_t::~SIDIS_2had_t(){
      if( xsec )
         delete xsec;
      if( var_thrower )
         delete var_thrower;
      if( decay_pi0 )
         delete decay_pi0;
   };

   // fill in all event variables, as the other
   // procedure fills in just minimal set of values
   int SIDIS_2had_t::ReconstructEvent_details(){
      int ierr = var.ReconstructEvent_SIDIS_2had( xsec->GetTargetPolState() );

      if( xsec->Get_Had1_PID() == PI_ZERO )
         (*decay_pi0)( var.had_1, var.gamma_1, var.gamma_2, r );
      else if( xsec->Get_Had2_PID() == PI_ZERO )
         (*decay_pi0)( var.had_2, var.gamma_1, var.gamma_2, r );

      if( check_Box_acceptance ){

         //cout << "LEPTON" << endl;
         int l2 = var.lep_2.InBoxAcc();

         //cout << "HADRON 1" << endl;
         int h1 = ( xsec->Get_Had1_PID() == PI_ZERO ? var.gamma_1.InBoxAcc() && var.gamma_2.InBoxAcc() : var.had_1.InBoxAcc() );

         //cout << "HADRON 2" << endl;
         int h2 = ( xsec->Get_Had2_PID() == PI_ZERO ? var.gamma_1.InBoxAcc() && var.gamma_2.InBoxAcc() : var.had_2.InBoxAcc() );

         //cout << l2 << ' ' << h1 << ' ' << h2 << endl;

         if( !ierr )
            ierr = !( l2 && h1 && h2 );

//          if( !ierr )
//             ierr = !var.lep_2.InBoxAcc();
//          if( !ierr )
//             ierr = !var.had_1.InBoxAcc();
//          if( !ierr )
//             ierr = !var.had_2.InBoxAcc();
      };
      return ierr;
   };


};
