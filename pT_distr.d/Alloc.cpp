
// allocates a pT or kT distribution


#include "pT_distr.d/pT_distr.h"
#include "pT_distr.d/Const_Exp.h"
#include "pT_distr.d/NonConst_Exp.h"
#include "pT_distr.d/Null_pT.h"
#include "pT_distr.d/Torino_Exp.h"

#include <iostream>
#include <string>
#include <sstream>

namespace TMDGen {
   namespace pT_distr {

      pT_distr_t* Alloc( const std::string& input_string ){

         std::stringstream ss( input_string );

         std::string word;
         ss >> word;

         bool is_kT = 0;

         if( word == "pT" ){
            is_kT = 0;
         } else if ( word == "kT" ){
            is_kT = 1;
         } else {
            std::cerr << "\t\tError in TMDGen::pT_distr::Alloc(...): string must start with 'pT ' or 'kT '" << std::endl;
            return 0;
         };

         ss >> word;

         pT_distr_t* pT_func = 0;

         if( word == "Exp" ){
            ss >> word;

            if( word == "Const" ){
               double alpha = 0.4;
               ss >> alpha;
               pT_func = new Const_Exp_t( alpha, is_kT );
               //std::cerr << "Allocating Const_Exp_t( " << alpha << ", " << is_kT << ")\tString = " << input_string << std::endl;
            } else if( word == "Torino" ){
               double alpha = 0.2;
               ss >> alpha;
               pT_func = new Torino_Exp_t( alpha, is_kT );
            } else if( word == "NonConst" ){
               int model = 0;
               double C = 0.4;
               ss >> model >> C;
               pT_func = new NonConst_Exp_t( model, C, is_kT );
               //std::cerr << "Allocating NonConst_Exp_t( " << model << ", " << C << ", " << is_kT << ")\tString = " << input_string << std::endl;
            } else {
               std::cerr << "\t\tUnknown type of Exp. pT (kT) distribution: '" << word << "'" << std::endl;
            };
         } else if (word == "NULL"){
            pT_func = new Null_pT_t();
         } else {
            std::cerr << "\t\tUnknown type of pT (kT) distribution: '" << word << "'" << std::endl;
         };

         return pT_func;
      };
   };
};
