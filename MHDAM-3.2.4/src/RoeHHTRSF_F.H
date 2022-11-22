#ifndef _ROEHHTRSF_F_H_
#define _ROEHHTRSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDROEHHTF 
#define GUARDROEHHTF 
// Prototype for Fortran procedure ROEHHTF ...
//
void FORTRAN_NAME( ROEHHTF ,roehhtf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_ROEHHTF FORTRAN_NAME( inlineROEHHTF, inlineROEHHTF)
#define FORTNT_ROEHHTF FORTRAN_NAME( ROEHHTF, roehhtf)

inline void FORTRAN_NAME(inlineROEHHTF, inlineROEHHTF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ROEHHTF");
 FORTRAN_NAME( ROEHHTF ,roehhtf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDROEHHTF 

}

#endif
