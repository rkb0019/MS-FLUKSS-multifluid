#ifndef _ROE7WAVEF_F_H_
#define _ROE7WAVEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRIEMANN7F 
#define GUARDRIEMANN7F 
// Prototype for Fortran procedure RIEMANN7F ...
//
void FORTRAN_NAME( RIEMANN7F ,riemann7f )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_RIEMANN7F FORTRAN_NAME( inlineRIEMANN7F, inlineRIEMANN7F)
#define FORTNT_RIEMANN7F FORTRAN_NAME( RIEMANN7F, riemann7f)

inline void FORTRAN_NAME(inlineRIEMANN7F, inlineRIEMANN7F)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANN7F");
 FORTRAN_NAME( RIEMANN7F ,riemann7f )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANN7F 

}

#endif
