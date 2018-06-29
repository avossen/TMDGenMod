#include "TMDGen_Wrapper.h"

#include "TMDGen.d/SIDIS_2had.h"
#include "Common.d/ParseInput.h"
#include "Common.d/Consts.h"

#include <iostream>
using std::cerr;
using std::endl;

// global variables

TMDGen::SIDIS_2had_t *gmc_trans = 0;
TMDGen::Var_t var;

// the functions

int init_tmdgen_( const char* TMDGen_instr_filename, int length ){

   cerr << "Length = " << length << endl;

   char filename[80];

   int i = 0;
   for( ; i<length && i<80 && TMDGen_instr_filename[i] != ' '; i++ ){
      filename[i] = TMDGen_instr_filename[i];
   };
   filename[i] = '\0';

   cerr << "Filename is --'" << filename << "'--" << endl;

   sgUtil::ParseInputReturn_t TMDGen_input;
   if( sgUtil::ParseInput( filename, TMDGen_input ) ){
      cerr << "Error loading instruction file." << endl;
      return 127;
   };

   // make TMDGen class
   TMDGen_input["Generate_with_Weights"] = "Yes";  // generate with weights
   TMDGen_input["MaxWeight"] = "Skip";             // do not estimate the maximum weight
   TMDGen_input["Output_To_Root_File"] = "";       // do not output

   if( gmc_trans )
      delete gmc_trans;

   gmc_trans = new TMDGen::SIDIS_2had_t ( TMDGen_input );


   // initialize TMDGen
   cerr << "\tInitializing TMDGen/TMDGen" << endl;
   if( gmc_trans->Initialize() ){
      cerr << "ERROR" << endl;
      return 127;
   };

   double E_beam_lab = 27.57;
   var.lep_0.E = E_beam_lab;
   var.lep_0.P = E_beam_lab;
   var.lep_0.M = ELECTRON_MASS;
   var.lep_1.E = E_beam_lab;
   var.lep_1.P = E_beam_lab;
   var.lep_1.M = ELECTRON_MASS;
   var.P_L = 0;

   return 0;
};


void eval_tmdgen_xsec_(){

   var.x = event_.x;
   var.y = event_.y;
   var.z = event_.z;
   var.P_hperp = event_.P_hperp;
   var.had_0.M = event_.M_h;
   var.cos_vartheta = event_.cos_vartheta;
   var.phi_h = event_.phi_h;
   var.phi_R = event_.phi_R;
   var.phi_S = event_.phi_S;
   var.psi = event_.e2_rPhi;
   var.P_T = event_.P_T;

   gmc_trans->WeightEvent( var );

   event_.w_int = var.weight_int;
   event_.w_full = var.weight;

   //cerr << "C++: out weights are " << event_.w_int << " and " << event_.w_full << endl;

};

void clear_tmdgen_(){
   if( gmc_trans )
      delete gmc_trans;
};
