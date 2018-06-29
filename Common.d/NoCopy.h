
// include these as private members to prevent copy const. and equals op.
#define NO_COPY_CONSTR( class ) class( const class& other ){  };
#define NO_COPY_CONSTR_W_CONST( class ) class( const class& other )
#define NO_EQ_OP( class )       class& operator=( const class& ){ return *this; };
