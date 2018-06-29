
//double dot;

/* 

   Note: initial state radiation assumed already included in lep_1
         final state radiation can be added to lep_2 (to yeild lep_3) post this class

*/

//TODO:
//    phi_S -> psi for independent

#include "Common.d/Consts.h"
#include "Common.d/Enums.h"
#include "Common.d/Exceptions.h"
#include "Common.d/FourMom.h"
#include "Common.d/Var.h"

#include <iostream>
#include <cmath>
#include <string>

namespace TMDGen {
   // initialize everything
   Var_t::Var_t( std::string name ) : flavor(ANTI_BOT_FLAV), x(0), y(0), Q2(0), W2(0), nu(0), phi_h(0), psi(0), phi_S(0), 
                                      P_L(0), P_T(0), S_L(0), S_T(0), 
                                      cos_theta_gamma(0), sin_theta_gamma(0),
                                      cos_theta_gamma_hadron(0), x_F(0),
                                      z(0), P_hperp(0), pT(0), kT(0), phi_pT(0), phi_kT(0),
                                      vartheta(0), cos_vartheta(0), sin_vartheta(0), phi_R(0), had_0( had ),
                                      t(0), mtprime(0),
                                      E_lep_in( lep_1.E ),
                                      weight(-999), thrown_weight(-999), weight_int(-999), XS(-999), XS_int(-999), integrating(0) {
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "*** Constructing '" << name << "', #" << ++N_INSTANTIATIONS << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
   };

   Var_t::~Var_t(){
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "*** Deconstructing #" << N_INSTANTIATIONS-- << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
//       std::cerr << "***" << std::endl;
   };

   int Var_t::N_INSTANTIATIONS = 0;

   // input: E_lep_in, x, y, z, P_hperp, had.M, P_T, P_L
   // output: Q2, nu, W2, had.E, had.P, kT, gamma, phi_S, S_T, S_L
   int Var_t::ComputeOtherVar_SIDIS_1had( const Var_t& min, const Var_t& max ) {
      int ierr = 1;

      Q2 = x*y*TWICE_PROTON_MASS * E_lep_in;
      //std::cout << "Q2 = " << Q2 << ' ' << x << ' ' << y << ' ' << E_lep_in << std::endl;

//       std::cerr << "E_lep_in = " << E_lep_in << std::endl;
//       std::cerr << &E_lep_in << ' ' << &lep_1.E << ' ' << this << std::endl;
//      std::cerr << "Q2: " << Q2 << ' ' << min.Q2 << ' ' << max.Q2 << std::endl;

      gamma_sq = 0;
      W2 = 0;

      if( Q2 > min.Q2 && Q2 < max.Q2 ){

         nu = E_lep_in * y;
         W2 = PROTON_MASS_SQUARED + TWICE_PROTON_MASS * nu - Q2;

         gamma_sq = 4. * PROTON_MASS_SQUARED * x * x / Q2;

         //std::cerr << "I AM HERE" << std::endl;
         //std::cerr << "W2: " << W2 << ' ' << min.W2 << ' ' << max.W2 << std::endl;

         // std::cerr << "W2 = " << W2 << std::endl;

         if( W2 > min.W2 && W2 < max.W2 ){

            had.E = z * y * E_lep_in;
            had.P = sqrt( had.E*had.E - had.M*had.M );

            if( P_hperp < had.P ){

               double kT_x = pT * cos( phi_pT ) - P_hperp / z * cos( phi_h );   // 03 Jun 2010--noticed missing 1/z
               double kT_y = pT * sin( phi_pT ) - P_hperp / z * sin( phi_h );

               kT = sqrt( kT_x*kT_x + kT_y*kT_y );
               phi_kT = atan2( kT_y, kT_x );

               if( kT > min.kT && kT < max.kT ){
                  ierr = 0;
               };

            };
         };
      };

      // zero in case of unpolarized
      phi_S = 0;
      S_T = 0;
      S_L = 0;

      //std::cerr << "ierr = " << ierr << std::endl;
      //throw ierr;

      // if no error
      if( !ierr ){
         ierr = 1;

         // working towards polarization values and phi_S

         double costheta = 1 - 0.5*Q2 / (1-y) / E_lep_in / E_lep_in;
         if( costheta >= 0 && costheta <= 1 ){

            // scattered electron theta in Ib
            // will later transform lep_2 to Ia
            lep_2_theta_Ib = acos( costheta );

            // scattered beam energy
            lep_2.E = E_lep_in - nu;

            // virtual photon 3 momentum
            double q = sqrt( nu*nu + Q2 );

            // the sin and cos of theta_gamma
            sin_theta_gamma = -lep_2.E / q * sin( lep_2_theta_Ib );

            if( sin_theta_gamma >= -1 && sin_theta_gamma <= 1 ){
               cos_theta_gamma = sqrt( 1 - sin_theta_gamma*sin_theta_gamma );

               // precompute some trig functions
               double cos_psi = cos( psi );
               double sin_psi = sin( psi );
               double cos_theta_e1 = cos( lep_1.theta );
               double sin_theta_e1 = sin( lep_1.theta );
               double cos_phi_e1 = cos( lep_1.phi );
               double sin_phi_e1 = sin( lep_1.phi );

               double S_X = 0, S_Y = 0;
               S_L = 0;

               if( P_L ){
                  S_X = cos_theta_gamma * cos_psi * sin_theta_e1 + sin_theta_gamma * cos_theta_e1;
                  S_X *= P_L;

                  S_Y = P_L * sin_psi * sin_theta_e1;

                  S_L = cos_theta_gamma * cos_theta_e1 - sin_theta_gamma * cos_psi * sin_theta_e1;
                  S_L *= P_L;
               };
               if( P_T ){
                  double temp = sin_theta_gamma * sin_theta_e1 * sin_phi_e1
                     - cos_theta_gamma * cos_psi * cos_theta_e1 * sin_phi_e1
                     + cos_theta_gamma * sin_psi * cos_phi_e1;
                  temp *= P_T;
                  S_X += temp;

                  temp = sin_psi * cos_theta_e1 * sin_phi_e1 - cos_psi * cos_phi_e1;
                  temp *= P_T;
                  S_Y += temp;

                  temp = sin_theta_gamma * cos_psi * cos_theta_e1 * sin_phi_e1
                     + cos_theta_gamma * sin_theta_e1 * sin_phi_e1
                     - sin_theta_gamma * sin_psi * cos_phi_e1;
                  temp *= P_T;

                  S_L -= temp;
               };

               if( P_L || P_T ){
                  phi_S = atan2( S_Y, S_X );
                  if( phi_S < 0 )
                     phi_S += TWO_PI;
                  S_T = sqrt( S_X*S_X + S_Y*S_Y );
                  //if( P_T < 0 )
                  //   S_T = -S_T;
                  //std::cout << "P_T, S_T " << P_T << ' ' << S_T << std::endl;
               };

               ierr = 0;
            };
         };
      };

      return ierr;
   };

   // input: E_lep_in, x, y, z, P_hperp, had.M, had_1.M, had_2.M
   // output: all of ComputeOtherVar_SIDIS_1had plus  R_3mag, R_4mag
   int Var_t::ComputeOtherVar_SIDIS_2had( const Var_t& min, const Var_t& max ) {

      R_3mag = 0;
      R_4mag_sq = 0;

      int ierr = ComputeOtherVar_SIDIS_1had( min, max );

      if( !ierr ){
         vartheta = acos( cos_vartheta );
         //cos_vartheta = cos(vartheta);
         sin_vartheta = sin(vartheta);

//          std::cerr << "had_0.M = " << had_0.M << std::endl;
//          std::cerr << "had_1.M = " << had_1.M << std::endl;
//          std::cerr << "had_2.M = " << had_2.M << std::endl;

         double m0_sq = had_0.M*had_0.M;
         double m1_sq = had_1.M*had_1.M;
         double m2_sq = (had_2.M == had_1.M) ? m1_sq : had_2.M*had_2.M;

         double R_3mag_sq = 0.25*m0_sq - 0.5*(m1_sq+m2_sq);   // is negative

         R_0_IIIa = (m1_sq - m2_sq) * 0.5 / had_0.M;

         R_4mag_sq = R_0_IIIa*R_0_IIIa - R_3mag_sq;

         if( R_3mag_sq > 0 )
            R_3mag = sqrt( R_3mag_sq );
      };

      //std::cerr << "ierr == " << ierr << std::endl;
      return ierr;
   };



   // reconstruct event, based on SIDIS variables, single hadron
   int Var_t::ReconstructEvent_SIDIS_1had( Enum::pol_state_t target_pol_state ){
      int ierr = 1;

      // scattered lepton
      {
         //double lep_p_Ia[3] = { 0, 0, 0 };

         double lep_p_Ib[3] = {
            lep_2.E*sin( lep_2_theta_Ib ),
            0,
            lep_2.E*cos( lep_2_theta_Ib )
         };

         Transform_Ib_to_Ia( lep_p_Ib, lep_p_Ia );
         Transform_Cart_to_Sph( lep_p_Ia, lep_2 );

      };


      // direction of hadron momentum

      // precompute some things

      // note negative because...
      sin_theta_gamma_hadron = P_hperp/had.P;   // already have checks to ensure this quantity is in (-1,1)
      cos_theta_gamma_hadron = sqrt( 1 - sin_theta_gamma_hadron*sin_theta_gamma_hadron );  // assume cos_theta_gamma_hadron < pi/2


      // direction of 3-momenta in Diehl/Sapeta unprimed coordinants (x, y, z)
      // same as Gliske's Ic frame
      double P_h_Ic[3] = { had.P*cos(phi_h)*sin_theta_gamma_hadron,
                           had.P*sin(phi_h)*sin_theta_gamma_hadron,
                           had.P*cos_theta_gamma_hadron };

      // save quantities
      // note x & y unchanged by boost between I & II
      P_h_IIa_x = P_h_Ic[0];
      P_h_IIa_y = P_h_Ic[1];

      // rotate to System Ia (detector coordinant system)
      Transform_Ic_to_Ia( P_h_Ic, P_h_Ia );

//       std::cout << "P_h_Ic = " << P_h_Ic[0] << ' ' << P_h_Ic[1] << ' ' << P_h_Ic[2] << std::endl;
//       std::cout << "P_h_Ia = " << P_h_Ia[0] << ' ' << P_h_Ia[1] << ' ' << P_h_Ia[2] << std::endl;

      // convert to spherical coordinants
      Transform_Cart_to_Sph( P_h_Ia, had );

      ierr = 0;

      // xF
      {
         q = ( nu*nu + Q2 ); // really |q|^2 at this point

         ierr = 1;
         if( q > 0 ){
            q = sqrt( q );

            double E_tot = nu + PROTON_MASS;
            double beta2 = q*q/E_tot/E_tot;
            if( beta2 < 1.0 ){
               ierr = 0;

               double eta = 1.0 / sqrt( 1.0 - beta2 );

               double had_dot_gamma = cos_theta_gamma_hadron * q * had.P;
               double had_dot_v = had_dot_gamma/E_tot;
               double factor = (((eta-1.0)/beta2)*(had_dot_v) - eta*had.E);

               double had_long = had_dot_gamma / q + factor/E_tot*q;
               x_F = 2*had_long/sqrt(W2);
            };
         };
      };


//       {
//          std::cout << "** " << sin_theta_gamma*q << " = " << lep_2.E*sin( lep_2_theta_Ib ) << std::endl;

//          std::cout << "\t\t\tP_h_Ia: " << P_h_Ia[0] << ' ' << P_h_Ia[1] << ' ' << P_h_Ia[2] << std::endl;
//          std::cout << "\t\t\tgamma: " << -lep_p_Ia[0] << ' ' << -lep_p_Ia[1] << ' ' << lep_1.E-lep_p_Ia[2] << std::endl;
//          std::cout << "\t\t\tgamma: " << q*sin_theta_gamma << ' ' << 0 << ' ' << q*cos_theta_gamma << std::endl;
//          dot = P_h_Ia[0]*(-lep_p_Ia[0]) + P_h_Ia[1]*(-lep_p_Ia[1]) + P_h_Ia[2]*(lep_1.E-lep_p_Ia[2]);
//       };

      return ierr;
   };

   // reconstruct event, based on SIDIS variables, dihadron
   int Var_t::ReconstructEvent_SIDIS_2had( Enum::pol_state_t target_pol_state ){

      int ierr = ReconstructEvent_SIDIS_1had( target_pol_state );

      // need to compute measured hadron momenta from parent hadron momenta
      if( !ierr ){
         // start with R in System IIIa

         //double sin_vartheta = sqrt( 1 - cos_vartheta*cos_vartheta );

         //std::cout << "a " << thrown_weight << ' ' << weight << ' ' << vartheta << ' ' << cos_vartheta << std::endl;
         //std::cout << "b " << thrown_weight/(1. + 0.5*cos_vartheta) << std::endl;

         double R_IIIa[3] = { R_3mag*cos(phi_R)*sin_vartheta,
                              R_3mag*sin(phi_R)*sin_vartheta,
                              R_3mag*cos_vartheta };


         // need to transform P_h to IIb 

         double W = sqrt( W2 );
         double gamma_I_II = (nu + PROTON_MASS) / W;
         double beta_gamma_I_II =  - q / W;

         double P_h_IIa_z = gamma_I_II * had_0.P * cos_theta_gamma_hadron + beta_gamma_I_II * had_0.E;
         double E_h_II = gamma_I_II * had_0.E                          + beta_gamma_I_II * had_0.P * cos_theta_gamma_hadron;

         // z in IIb is the magnitude of the 3 mom of P_h in IIa
         double P_h_IIb_z = sqrt( P_hperp*P_hperp + P_h_IIa_z*P_h_IIa_z );


         // are now prepared to transform R from IIIa to IIb

         double gamma_II_III = E_h_II/had_0.M;
         double gamma_beta_II_III = P_h_IIb_z/had_0.M;

         double R_0_II = gamma_II_III*R_0_IIIa + gamma_beta_II_III*R_IIIa[2];
         double R_IIb[3] = { R_IIIa[0],
                             R_IIIa[1],
                             gamma_II_III*R_IIIa[2] + gamma_beta_II_III*R_0_IIIa };

         // transform IIb to IIa
         double P_h_IIa_xz = sqrt( P_h_IIa_x*P_h_IIa_x + P_h_IIa_z*P_h_IIa_z );

         double Rot_IIb_IIa[3][3] =
            {{  P_h_IIa_z/P_h_IIa_xz, -P_h_IIa_x * P_h_IIa_y / P_h_IIa_xz / P_h_IIb_z, P_h_IIa_x/P_h_IIb_z },
             {                     0,                          P_h_IIa_xz / P_h_IIb_z, P_h_IIa_y/P_h_IIb_z },
             { -P_h_IIa_x/P_h_IIa_xz, -P_h_IIa_y * P_h_IIa_z / P_h_IIa_xz / P_h_IIb_z, P_h_IIa_z/P_h_IIb_z } };

         double R_IIa[3] = { Rot_IIb_IIa[0][0]*R_IIb[0] + Rot_IIb_IIa[0][1]*R_IIb[1] + Rot_IIb_IIa[0][2]*R_IIb[2],
                             Rot_IIb_IIa[1][0]*R_IIb[0] + Rot_IIb_IIa[1][1]*R_IIb[1] + Rot_IIb_IIa[1][2]*R_IIb[2],
                             Rot_IIb_IIa[2][0]*R_IIb[0] + Rot_IIb_IIa[2][1]*R_IIb[1] + Rot_IIb_IIa[2][2]*R_IIb[2] };

         // transform IIa to Ic
         double R_0_I = gamma_I_II * R_0_II - beta_gamma_I_II * R_IIa[2];

         double R_Ic[3] = {
            R_IIa[0],
            R_IIa[1],
            gamma_I_II * R_IIa[2] - beta_gamma_I_II * R_0_II };


         // transform from Ic to Ia
         double R_Ia[3] = {0,0,0};

         Transform_Ic_to_Ia( R_Ic, R_Ia );

         double p1[3] = { 0.5*P_h_Ia[0] + R_Ia[0],
                          0.5*P_h_Ia[1] + R_Ia[1],
                          0.5*P_h_Ia[2] + R_Ia[2] };

         double p2[3] = { 0.5*P_h_Ia[0] - R_Ia[0],
                          0.5*P_h_Ia[1] - R_Ia[1],
                          0.5*P_h_Ia[2] - R_Ia[2] };

         // convert to spherical
         Transform_Cart_to_Sph( p1, had_1 );
         Transform_Cart_to_Sph( p2, had_2 );

         // SANITY CHECKS

         // double check energy (equivelant to checking mass)
         had_1.E = sqrt( had_1.P * had_1.P + had_1.M * had_1.M );
         had_2.E = sqrt( had_2.P * had_2.P + had_2.M * had_2.M );

         if( fabs((float)(had_1.E + had_2.E) - (float)(had_0.E)) > 1e-5 ){
            std::cerr << had_1.E << " + " << had_2.E << " != " << had_0.E << std::endl;
            std::cerr << "Delta = " << fabs((float)(had_1.E + had_2.E) - (float)(had_0.E)) << std::endl;

            throw Error::SanityCheckFailure("Var_t::ReconstructEvent_SIDIS_2had(...)", "Inconsistant result when computing hadron system energy.");
         };


      };


      //       std::cout << "\t\t\t\t" << dot << ' ' << dot/q/had_0.P << ' ' << cos_theta_gamma_hadron << std::endl;

//       dot = P_h_Ia[0]*(-lep_p_Ia[0]) + P_h_Ia[1]*(-lep_p_Ia[1]) + P_h_Ia[2]*(lep_1.E-lep_p_Ia[2]);

//       std::cout << "\t\t\t\t" << dot << ' ' << dot/q/had_0.P << ' ' << cos_theta_gamma_hadron << std::endl;

//        // check in 1b gamma is as expected 

//        double lep_p_Ib[3] = {
//           lep_2.E*sin( lep_2_theta_Ib ),
//           0,
//           lep_2.E*cos( lep_2_theta_Ib )
//        };
//        double gamma_Ib_a[3] = {





      return ierr;
   };


   // transform 3 momentum from System Ic to Ia ( just rotations )
   void Var_t::Transform_Ic_to_Ia( const double p_Ic[3], double p_Ia[3] ){

      // rotate to Ib system
      double p_Ib[3] = { cos_theta_gamma*p_Ic[0]+sin_theta_gamma*p_Ic[2], 
                        p_Ic[1],
                        -sin_theta_gamma*p_Ic[0]+cos_theta_gamma*p_Ic[2]  };

//       // check cos_theta_gamma_had
//       {
//          double p = sqrt( p_Ib[0]*p_Ib[0] + p_Ib[1]*p_Ib[1] + p_Ib[2]*p_Ib[2] );
//          double gamma[3] = { q*sin_theta_gamma, 0, q*cos_theta_gamma };
//          double dot = p_Ib[0]*gamma[0] + p_Ib[1]*gamma[1] + p_Ib[2]*gamma[2];

//          if( p > 3)
//             std::cout << p << ' ' << dot/q/p << ' ' << cos_theta_gamma_hadron << std::endl;
//       };

      double p = sqrt( p_Ib[0]*p_Ib[0] + p_Ib[1]*p_Ib[1] + p_Ib[2]*p_Ib[2] );

      // rotate to A system
      Transform_Ib_to_Ia( p_Ib, p_Ia );

//       if( p > 3 ){
//          std::cout << "P_h_Ic = " << p_Ic[0] << ' ' << p_Ic[1] << ' ' << p_Ic[2] << std::endl;
//          std::cout << "P_h_Ib = " << p_Ib[0] << ' ' << p_Ib[1] << ' ' << p_Ib[2] << std::endl;
//          std::cout << "P_h_Ia = " << p_Ia[0] << ' ' << p_Ia[1] << ' ' << p_Ia[2] << std::endl;
//       };


   };

   // transform 3 momentum from System Ib to Ia ( just rotations )
   void Var_t::Transform_Ib_to_Ia( const double p_Ib[3], double p_Ia[3] ){

      // rotate about z
      // 23 August, 2010: changed sign of psi
      // now lep_2.phi = psi when lep_1 == lep_0
      // 03 March, 2011: changed sign of psi back
      // should have lep_2.phi = -psi when lep_1 == lep_0
      double p_temp_1[3] = {
         cos(psi)*p_Ib[0] + sin(psi)*p_Ib[1],
        -sin(psi)*p_Ib[0] + cos(psi)*p_Ib[1],
         p_Ib[2]
      };

      // rotate about y
      double p_temp_2[3] = {
          cos(lep_1.theta)*p_temp_1[0] + sin(lep_1.theta)*p_temp_1[2],
         p_temp_1[1],
         -sin(lep_1.theta)*p_temp_1[0] + cos(lep_1.theta)*p_temp_1[2]
      };

      // rotate about z
      p_Ia[0] = cos(lep_1.phi)*p_temp_2[0] - sin(lep_1.phi)*p_temp_2[1];
      p_Ia[1] = sin(lep_1.phi)*p_temp_2[0] + cos(lep_1.phi)*p_temp_2[1];
      p_Ia[2] = p_temp_2[2];

   };

   // transform Cartesian to Spherical coordinants
   void Var_t::Transform_Cart_to_Sph( const double p[3], FourMom_t& mom ){
      mom.P = sqrt( p[0]*p[0] + p[1]*p[1] + p[2]*p[2] );
      mom.theta = 0;
      mom.phi = 0;

      if( mom.P ){
         mom.theta = acos( p[2]/mom.P );
         mom.phi = atan2( p[1], p[0] );
         if( mom.phi < 0 )
            mom.phi += TWO_PI;
      };

      //std::cout << p[0] << ' ' << p[1] << ' ' << p[2] << ' ' << mom.P << ' ' << mom.theta << ' ' << mom.phi << std::endl;
   };

   void Var_t::SmallCopy( const Var_t& var_in ){
      x = var_in.x;
      y = var_in.y;
      z = var_in.z;
      Q2 = var_in.Q2;
      W2 = var_in.W2;
      P_hperp = var_in.P_hperp;
      had_0.M = var_in.had_0.M;
      kT = var_in.kT;

   };


};
