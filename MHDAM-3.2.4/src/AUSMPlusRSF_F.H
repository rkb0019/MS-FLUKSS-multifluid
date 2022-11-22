#ifndef _AUSMPLUSRSF_F_H_
#define _AUSMPLUSRSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAUSMPLUSF 
#define GUARDAUSMPLUSF 
// Prototype for Fortran procedure AUSMPLUSF ...
//
void FORTRAN_NAME( AUSMPLUSF ,ausmplusf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_AUSMPLUSF FORTRAN_NAME( inlineAUSMPLUSF, inlineAUSMPLUSF)
#define FORTNT_AUSMPLUSF FORTRAN_NAME( AUSMPLUSF, ausmplusf)

inline void FORTRAN_NAME(inlineAUSMPLUSF, inlineAUSMPLUSF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_AUSMPLUSF");
 FORTRAN_NAME( AUSMPLUSF ,ausmplusf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAUSMPLUSF 

}

#endif
