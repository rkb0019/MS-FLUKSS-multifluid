#ifndef _VANLEERRSF_F_H_
#define _VANLEERRSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVANLEERF 
#define GUARDVANLEERF 
// Prototype for Fortran procedure VANLEERF ...
//
void FORTRAN_NAME( VANLEERF ,vanleerf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_VANLEERF FORTRAN_NAME( inlineVANLEERF, inlineVANLEERF)
#define FORTNT_VANLEERF FORTRAN_NAME( VANLEERF, vanleerf)

inline void FORTRAN_NAME(inlineVANLEERF, inlineVANLEERF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANLEERF");
 FORTRAN_NAME( VANLEERF ,vanleerf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANLEERF 

#ifndef GUARDVANLEERRS 
#define GUARDVANLEERRS 
// Prototype for Fortran procedure VANLEERRS ...
//
void FORTRAN_NAME( VANLEERRS ,vanleerrs )(
      CHFp_CONST_REAL(PGasL)
      ,CHFp_CONST_REAL(RhoL)
      ,CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(VL)
      ,CHFp_CONST_REAL(WL)
      ,CHFp_CONST_REAL(PGasR)
      ,CHFp_CONST_REAL(RhoR)
      ,CHFp_CONST_REAL(UR)
      ,CHFp_CONST_REAL(VR)
      ,CHFp_CONST_REAL(WR)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUX)
      ,CHFp_REAL(FLUXRUY)
      ,CHFp_REAL(FLUXRUZ)
      ,CHFp_REAL(FLUXE) );

#define FORT_VANLEERRS FORTRAN_NAME( inlineVANLEERRS, inlineVANLEERRS)
#define FORTNT_VANLEERRS FORTRAN_NAME( VANLEERRS, vanleerrs)

inline void FORTRAN_NAME(inlineVANLEERRS, inlineVANLEERRS)(
      CHFp_CONST_REAL(PGasL)
      ,CHFp_CONST_REAL(RhoL)
      ,CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(VL)
      ,CHFp_CONST_REAL(WL)
      ,CHFp_CONST_REAL(PGasR)
      ,CHFp_CONST_REAL(RhoR)
      ,CHFp_CONST_REAL(UR)
      ,CHFp_CONST_REAL(VR)
      ,CHFp_CONST_REAL(WR)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUX)
      ,CHFp_REAL(FLUXRUY)
      ,CHFp_REAL(FLUXRUZ)
      ,CHFp_REAL(FLUXE) )
{
 CH_TIMELEAF("FORT_VANLEERRS");
 FORTRAN_NAME( VANLEERRS ,vanleerrs )(
      CHFt_CONST_REAL(PGasL)
      ,CHFt_CONST_REAL(RhoL)
      ,CHFt_CONST_REAL(UL)
      ,CHFt_CONST_REAL(VL)
      ,CHFt_CONST_REAL(WL)
      ,CHFt_CONST_REAL(PGasR)
      ,CHFt_CONST_REAL(RhoR)
      ,CHFt_CONST_REAL(UR)
      ,CHFt_CONST_REAL(VR)
      ,CHFt_CONST_REAL(WR)
      ,CHFt_REAL(FLUXRHO)
      ,CHFt_REAL(FLUXRUX)
      ,CHFt_REAL(FLUXRUY)
      ,CHFt_REAL(FLUXRUZ)
      ,CHFt_REAL(FLUXE) );
}
#endif  // GUARDVANLEERRS 

}

#endif
