
#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include "Common.d/FourMom.h"
#include "RNG.d/GSL_RNG.h"
#include "Decay.d/pi0_ROOT.h"


namespace TMDGen {
   namespace Decay {

      // assume cartesian coordinants are valid for pi0
      int pi0_ROOT::operator() ( const FourMom_t& p_pi0, FourMom_t& gamma_1, FourMom_t& gamma_2, TMDGen::GSL_RNG_t *r ){

         double p_sin = p_pi0.P*sin( p_pi0.theta );

         TLorentzVector p_pi0_ROOT( p_sin*cos( p_pi0.phi ),
                                    p_sin*sin( p_pi0.phi ),
                                    p_pi0.P*cos( p_pi0.theta ),
                                    p_pi0.E );

         double mass[2] = { 0, 0 };
         int ierr = !genps.SetDecay( p_pi0_ROOT, 2, mass );

         if( !ierr ){
            // get a valid event via acceptance/rejection
            while ( genps.Generate() < genps.GetWtMax() * r->EvalUnif() ) { /* */ };

            TLorentzVector *g1 = genps.GetDecay(0);
            TLorentzVector *g2 = genps.GetDecay(1);

            gamma_1.E = g1->E();
            gamma_1.M = 0;
            gamma_1.P = g1->Rho();
            gamma_1.theta = g1->Theta();
            gamma_1.phi = g1->Phi();

            gamma_2.E = g2->E();
            gamma_2.M = 0;
            gamma_2.P = g2->Rho();
            gamma_2.theta = g2->Theta();
            gamma_2.phi = g2->Phi();
         };
 
         return ierr;
      };

   };
};

