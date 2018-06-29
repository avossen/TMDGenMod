
/*
  Struct to keep track of parameters
*/

#include "Common.d/Flavor.h"
#include "DF.d/h1_Torino_Params.h"
#include "Common.d/Exceptions.h"

#include <string>


namespace TMDGen {
   namespace DF {

      h1_Torino_Params_t::h1_Torino_Params_t( std::string init_code ){

         for( double *ptr = &N[0]; ptr != &N[ GMC_TRANS_N_FLAVORS ]; ++ptr )
            (*ptr) = 0;

         if( init_code == "Ringberg" ){
            // from arXiv:0812.4366v1
            // Invited talk delivered by U. D'Alesio at the Ringberg Workshop
            // "New Trends in HERA Physics 2008", Ringberg Castle, Tegernsee, Germany,
            // October 5-10, 2008.

            N[ UP_FLAV ] = 0.64;
            N[ DOWN_FLAV ] = -1.00;
            alpha = 0.73;
            beta = 0.84;
            ave_pT2 = 0.25;

          } else if( init_code == "Ringberg-b" ){
             // from arXiv:0812.4366v1
             // Invited talk delivered by U. D'Alesio at the Ringberg Workshop
             // "New Trends in HERA Physics 2008", Ringberg Castle, Tegernsee, Germany,
             // October 5-10, 2008.
             N[ UP_FLAV ] = 0.64;
             N[ DOWN_FLAV ] = -1.00;
             alpha = 0.73;
             beta = 0.84;
             ave_pT2 = 0.18;   // changed to closer ally with ave pT for dihadrons

         } else if ( init_code == "PRD_07_A12" ){
            // from 
            // Journal reference: Phys.Rev.D75:054032,2007
            // DOI: 10.1103/PhysRevD.75.054032
            // Cite as: arXiv:hep-ph/0701006v3

            N[ UP_FLAV ] = 0.48;
            N[ DOWN_FLAV ] = -0.62;
            alpha = 1.14;
            beta = 4.74;
            ave_pT2 = 0.25;

         } else if ( init_code == "PRD_07_A0" ){
            // also from 
            // Journal reference: Phys.Rev.D75:054032,2007
            // DOI: 10.1103/PhysRevD.75.054032
            // Cite as: arXiv:hep-ph/0701006v3

            N[ UP_FLAV ] = 0.42;
            N[ DOWN_FLAV ] = -0.53;
            alpha = 1.2;
            beta = 5.09;
            ave_pT2 = 0.25;

         } else {
            throw Error::Constructing("h1_Torino_Params_t", std::string("Unknown paramter set code: '") + init_code + "'" );
         };
      };

   };
};
