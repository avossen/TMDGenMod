/*
   Abstract base class for random number generation
*/

#ifndef _RNG_H_
#define _RNG_H_


namespace TMDGen {

   class RNG_t {
   public:
      virtual ~RNG_t(){ /* */ };
      virtual double EvalUnif() = 0;
   };

};

#endif
