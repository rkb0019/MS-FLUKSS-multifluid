#ifndef _RUSANOVF_F_H_
#define _RUSANOVF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRUSANOVF 
#define GUARDRUSANOVF 
// Prototype for Fortran procedure RUSANOVF ...
//
void FORTRAN_NAME( RUSANOVF ,rusanovf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_RUSANOVF FORTRAN_NAME( inlineRUSANOVF, inlineRUSANOVF)
#define FORTNT_RUSANOVF FORTRAN_NAME( RUSANOVF, rusanovf)

inline void FORTRAN_NAME(inlineRUSANOVF, inlineRUSANOVF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RUSANOVF");
 FORTRAN_NAME( RUSANOVF ,rusanovf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRUSANOVF 

}

#endif
