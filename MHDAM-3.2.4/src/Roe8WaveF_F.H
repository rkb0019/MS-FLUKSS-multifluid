#ifndef _ROE8WAVEF_F_H_
#define _ROE8WAVEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRIEMANNF 
#define GUARDRIEMANNF 
// Prototype for Fortran procedure RIEMANNF ...
//
void FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNF FORTRAN_NAME( inlineRIEMANNF, inlineRIEMANNF)
#define FORTNT_RIEMANNF FORTRAN_NAME( RIEMANNF, riemannf)

inline void FORTRAN_NAME(inlineRIEMANNF, inlineRIEMANNF)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNF");
 FORTRAN_NAME( RIEMANNF ,riemannf )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNF 

#ifndef GUARDROE8WAVEF 
#define GUARDROE8WAVEF 
// Prototype for Fortran procedure ROE8WAVEF ...
//
void FORTRAN_NAME( ROE8WAVEF ,roe8wavef )(
      CHFp_CONST_REAL(PGasL)
      ,CHFp_CONST_REAL(RhoL)
      ,CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(VL)
      ,CHFp_CONST_REAL(WL)
      ,CHFp_CONST_REAL(BXL)
      ,CHFp_CONST_REAL(BYL)
      ,CHFp_CONST_REAL(BZL)
      ,CHFp_CONST_REAL(PGasR)
      ,CHFp_CONST_REAL(RhoR)
      ,CHFp_CONST_REAL(UR)
      ,CHFp_CONST_REAL(VR)
      ,CHFp_CONST_REAL(WR)
      ,CHFp_CONST_REAL(BXR)
      ,CHFp_CONST_REAL(BYR)
      ,CHFp_CONST_REAL(BZR)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUX)
      ,CHFp_REAL(FLUXRUY)
      ,CHFp_REAL(FLUXRUZ)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBX)
      ,CHFp_REAL(FLUXBY)
      ,CHFp_REAL(FLUXBZ) );

#define FORT_ROE8WAVEF FORTRAN_NAME( inlineROE8WAVEF, inlineROE8WAVEF)
#define FORTNT_ROE8WAVEF FORTRAN_NAME( ROE8WAVEF, roe8wavef)

inline void FORTRAN_NAME(inlineROE8WAVEF, inlineROE8WAVEF)(
      CHFp_CONST_REAL(PGasL)
      ,CHFp_CONST_REAL(RhoL)
      ,CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(VL)
      ,CHFp_CONST_REAL(WL)
      ,CHFp_CONST_REAL(BXL)
      ,CHFp_CONST_REAL(BYL)
      ,CHFp_CONST_REAL(BZL)
      ,CHFp_CONST_REAL(PGasR)
      ,CHFp_CONST_REAL(RhoR)
      ,CHFp_CONST_REAL(UR)
      ,CHFp_CONST_REAL(VR)
      ,CHFp_CONST_REAL(WR)
      ,CHFp_CONST_REAL(BXR)
      ,CHFp_CONST_REAL(BYR)
      ,CHFp_CONST_REAL(BZR)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUX)
      ,CHFp_REAL(FLUXRUY)
      ,CHFp_REAL(FLUXRUZ)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBX)
      ,CHFp_REAL(FLUXBY)
      ,CHFp_REAL(FLUXBZ) )
{
 CH_TIMELEAF("FORT_ROE8WAVEF");
 FORTRAN_NAME( ROE8WAVEF ,roe8wavef )(
      CHFt_CONST_REAL(PGasL)
      ,CHFt_CONST_REAL(RhoL)
      ,CHFt_CONST_REAL(UL)
      ,CHFt_CONST_REAL(VL)
      ,CHFt_CONST_REAL(WL)
      ,CHFt_CONST_REAL(BXL)
      ,CHFt_CONST_REAL(BYL)
      ,CHFt_CONST_REAL(BZL)
      ,CHFt_CONST_REAL(PGasR)
      ,CHFt_CONST_REAL(RhoR)
      ,CHFt_CONST_REAL(UR)
      ,CHFt_CONST_REAL(VR)
      ,CHFt_CONST_REAL(WR)
      ,CHFt_CONST_REAL(BXR)
      ,CHFt_CONST_REAL(BYR)
      ,CHFt_CONST_REAL(BZR)
      ,CHFt_REAL(FLUXRHO)
      ,CHFt_REAL(FLUXRUX)
      ,CHFt_REAL(FLUXRUY)
      ,CHFt_REAL(FLUXRUZ)
      ,CHFt_REAL(FLUXE)
      ,CHFt_REAL(FLUXBX)
      ,CHFt_REAL(FLUXBY)
      ,CHFt_REAL(FLUXBZ) );
}
#endif  // GUARDROE8WAVEF 

}

#endif
