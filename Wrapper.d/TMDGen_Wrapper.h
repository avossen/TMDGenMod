// Simple wrapper to allow FORTRAN to call evaluate the cross section from C++ TMDGen 

extern "C" {

   // structures to include variables
   extern struct Event_t {
      double x, y, z, P_hperp, M_h, cos_vartheta, phi_h, phi_R, e2_rPhi, phi_S, w_int, w_full;
      int P_T;
   } event_;

/*    extern struct Event_Out_t { */
/*       double w_int, w_full; */
/*    } event_out_; */


   // function to initialize the generator
   // requires filename of input file
   // returns error code (0 = all is OK)
   int init_tmdgen_( const char* TMDGen_instr_filename, int length );

   // Evaluates the cross section, returning both the angular integrated term and the full cross section
   void eval_tmdgen_xsec_();

   // frees memory allocated by the initialization routine
   void clear_tmdgen_();


};
