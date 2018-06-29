/*
  Main class for constructing SIDIS Cross Sections
*/

#ifndef _XSEC_SIDIS_2HAD_H_
#define _XSEC_SIDIS_2HAD_H_

#include "Common.d/ParseInput.h"
#include "Common.d/LundPID.h"
#include "Common.d/Var.h"

#include "XSec.d/XSec.h"

#include "DF.d/Full_DF_Set.h"
#include "FF.d/Full_FF_Set.h"

#include <vector>
#include <cmath>

namespace TMDGen {
   namespace XSec {
      class SIDIS_2had_t : public XSec_t {

         double had_1_M, had_2_M;
         LundPID_t had_1_PID, had_2_PID;

         NO_COPY_CONSTR_W_CONST( SIDIS_2had_t ) : XSec_t( other.min, other.max ) { /* */ };
         NO_EQ_OP( SIDIS_2had_t );

         // Determine valid domain
         virtual void DetermineDomain( const sgUtil::ParseInputReturn_t& parsed_input );

         // Computing dependent variables
         virtual int ComputeOtherVar( Var_t& );

         // for overall kinematic and constant factors for the cross section
         virtual double Overall_Kinematic_Factor( const Var_t& var ) const;

         // make needed terms
         virtual void MakeTerms( const sgUtil::ParseInputReturn_t& parsed_input );

         // to precomputed Legendre polynomials
         double cos_array[6];


      public:
         SIDIS_2had_t( const sgUtil::ParseInputReturn_t& input, Var_t& min_in, Var_t& max_in );
         virtual ~SIDIS_2had_t();

         LundPID_t Get_Had1_PID() const;
         LundPID_t Get_Had2_PID() const;

      };
   };
};


#endif
