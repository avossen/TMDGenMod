/*

   Struct to include current value of variables common to various processes
   Initializes values to zero
   Also contains member functions to compute certain variables based on other variables

*/

#ifndef _VAR_H_
#define _VAR_H_

#include "Common.d/Enums.h"
#include "Common.d/Flavor.h" 
#include "Common.d/FourMom.h"

#include <string>

namespace TMDGen {

   struct Var_t {
      flavor_t flavor;

      double x, y, Q2, W2, nu, phi_h, psi, phi_S;    // kinematic and angular variables

      double P_L, P_T, S_L, S_T;                     // polarization values

      FourMom_t lep_0,                               // lep_0 : beam lepton
         lep_1,                                      // lep_1 : incident lepton (beam after possible radiative effects)
         lep_2,                                      // lep_2 : scattered lepton after primary nuclean interaction
         lep_3,                                      // lep_3 : lep_2 after possible radiative effects
         had;                                        // had   : the produced hadron (or hadron system)

      double cos_theta_gamma, sin_theta_gamma;       // other temporary quantities
      double cos_theta_gamma_hadron, sin_theta_gamma_hadron;
      double q, x_F, gamma_sq;                       // q is mag. of virtual photon 3-mom
                                                     // gamma is 2Mx/Q
      double P_h_IIa_x, P_h_IIa_y;                   // P_h in IIa
      double P_h_Ia[3];                              // P_h in Ia
      double lep_p_Ia[3];                            // scattered lepton momenta in Ia
      double lep_2_theta_Ib;                         // angle between incident and scattered leptons


      double z, P_hperp, pT, kT, phi_pT, phi_kT;     // SIDIS 1 had

      double vartheta, cos_vartheta, sin_vartheta, phi_R;  // SIDIS 2 had
      double R_3mag, R_4mag_sq, R_0_IIIa;
      FourMom_t had_1, had_2, &had_0, gamma_1, gamma_2;

      double t, mtprime;                             // exclusive

      double vertex[3];                              // location (in lab coordinants) of primary vertex

      double &E_lep_in;                              // energy of incident lepton entering primary physical cross section
                                                     // actually a reference to lep_1.E

      double weight, thrown_weight, weight_int,
         XS, XS_int;                                 // Event level quantities
      int ISR, FSR;                                  // whether exists initial or final state radiation (lepton bremsstrahlung)


      // extra status-type things
      int integrating;

      static int N_INSTANTIATIONS;


      // member functions

      // initialize everything
      Var_t( std::string name = "default" );
      ~Var_t();

      // Based on the values of ( E_beam_lab, x, y, z, P_hperp, had.M )
      // this function computes ( Q2, nu, W2, had.E, had.P, kT )
      // also checks ranges on Q2 and W2 and ensures P_hperp <= had.P
      int ComputeOtherVar_SIDIS_1had( const Var_t& min, const Var_t& max );

      // calls ComputeOtherVar_SIDIS_1had, and then computes R_3mag, R_4mag, P_h_LCminus
      int ComputeOtherVar_SIDIS_2had( const Var_t& min, const Var_t& max );

      // compute scattered lepton and hadron system momenta from SIDIS variables (i.e. SIDIS 1 hadron case)
      int ReconstructEvent_SIDIS_1had( Enum::pol_state_t target_pol_state );

      // determine momenta of two hadrons and everything in Reconstruct_SIDIS_1had()
      int ReconstructEvent_SIDIS_2had( Enum::pol_state_t target_pol_state );

      // transform 3 momenta from System Ic to Ia
      void Transform_Ic_to_Ia( const double p_Ic[3], double p_Ia[3] );

      // transform 3 momenta from System Ib to Ia
      void Transform_Ib_to_Ia( const double p_Ib[3], double p_Ia[3] );

      // transform cartesion to spherical
      void Transform_Cart_to_Sph( const double p[3], FourMom_t& mom );

      void SmallCopy( const Var_t& var_in );

   };
};

#endif
