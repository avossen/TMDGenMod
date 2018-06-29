
/*
  Struct to keep track of parameters
*/

#include "Common.d/Flavor.h"
#include "DF.d/Torino_nT_odd_Params.h"
#include "Common.d/Exceptions.h"

#include <string>


namespace TMDGen {
   namespace DF {

      Torino_nT_odd_Params_t::Torino_nT_odd_Params_t( std::string init_code ){

         for( int i=0; i<GMC_TRANS_N_FLAVORS; ++i ){
            A[i] = 0;
            alpha[i] = 0;
            lambda[i] = 1;
         };

         // basic setup
         A[ UP_FLAV ] = -0.35;
         A[ DOWN_FLAV ] = 0.90;
         A[ STR_FLAV ] = 0.24;
         A[ ANTI_UP_FLAV ] = -0.04;
         A[ ANTI_DOWN_FLAV ] = 0.40;
         A[ ANTI_STR_FLAV ] = -1;

         alpha[ UP_FLAV ] = 0.73;
         alpha[ DOWN_FLAV ] = 1.08;
         alpha[ STR_FLAV ] = 0.79;
         alpha[ ANTI_UP_FLAV ] = 0.79;
         alpha[ ANTI_DOWN_FLAV ] = 0.79;
         alpha[ ANTI_STR_FLAV ] = 0.79;

         beta = 3.46;
         ave_pT_sq = 0.25;

         M1_sq = 0.34;

         if( init_code.substr(0, 8) == "h_1^perp" ){
            lambda[ UP_FLAV ] = 2.0;
            lambda[ DOWN_FLAV ] = -1.111;
            lambda[ STR_FLAV ] = -99999;        // flag to set to negative abs. value
            lambda[ ANTI_UP_FLAV ] = -99999;
            lambda[ ANTI_DOWN_FLAV ] = -99999;
            lambda[ ANTI_STR_FLAV ] = -99999;
         };


         if( init_code == "f_1T^perp_08" ){
            // nothing different than the basic

         } else if( init_code == "f_1T^perp_08-Fit-2" ){
            // to match h_1^perp_08-Fit-2
            ave_pT_sq = 0.18;

         } else if( init_code == "h_1^perp_08-Fit-1" ){
            // also nothing different from the basic

         } else if( init_code == "h_1^perp_08-Fit-2" ){

            lambda[ UP_FLAV ] = 2.1;
            ave_pT_sq = 0.18;


         } else {
            throw Error::Constructing("Torino_nT_odd_Params_t", std::string("Unknown paramter set code: '") + init_code + "'" );
         };
      };

   };
};
