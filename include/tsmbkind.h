!*    kindef: define default KIND macros
! --------------------------------------


USE PARKIND1, ONLY :&
 &JPIT, JPIS, JPIM, JPIB,&
 &JPRT, JPRS, JPRM, JPRB


#ifndef INTEGER_T
#define INTEGER_T INTEGER(KIND=JPIT)
#define INTEGER_S INTEGER(KIND=JPIS)
#define INTEGER_M INTEGER(KIND=JPIM)
#define INTEGER_B INTEGER(KIND=JPIB)
#endif

#ifndef REAL_T
#define REAL_T REAL(KIND=JPRT)
#define REAL_S REAL(KIND=JPRS)
#define REAL_M REAL(KIND=JPRM)
#define REAL_B REAL(KIND=JPRB)
#endif

#ifndef _0T
#define _0T  0.0_JPRT
#define _0S  0.0_JPRS
#define _0M  0.0_JPRM
#define _0B  0.0_JPRB
#endif

#ifndef _05T
#define _05T 0.5_JPRT
#define _05S 0.5_JPRS
#define _05M 0.5_JPRM
#define _05B 0.5_JPRB
#endif

#ifndef _1T
#define _1T  1.0_JPRT
#define _1S  1.0_JPRS
#define _1M  1.0_JPRM
#define _1B  1.0_JPRB
#endif

#ifndef _2T
#define _2T  2.0_JPRT
#define _2S  2.0_JPRS
#define _2M  2.0_JPRM
#define _2B  2.0_JPRB
#endif

#define _ZERO_   _0B
#define _HALF_   _05B
#define _ONE_    _1B
#define _TWO_    _2B

! --------------------------------------
