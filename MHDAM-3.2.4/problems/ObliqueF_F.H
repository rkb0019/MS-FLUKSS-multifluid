#ifndef _OBLIQUEF_F_H_
#define _OBLIQUEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETOBLIQUE 
#define GUARDSETOBLIQUE 
// Prototype for Fortran procedure SETOBLIQUE ...
//
void FORTRAN_NAME( SETOBLIQUE ,setoblique )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auAngleSW)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auObliqueAe)
      ,CHFp_CONST_REAL(auObliqueAv)
      ,CHFp_CONST_REAL(auObliquePsi)
      ,CHFp_CONST_REAL(auObliqueK) );

#define FORT_SETOBLIQUE FORTRAN_NAME( inlineSETOBLIQUE, inlineSETOBLIQUE)
#define FORTNT_SETOBLIQUE FORTRAN_NAME( SETOBLIQUE, setoblique)

inline void FORTRAN_NAME(inlineSETOBLIQUE, inlineSETOBLIQUE)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auAngleSW)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auObliqueAe)
      ,CHFp_CONST_REAL(auObliqueAv)
      ,CHFp_CONST_REAL(auObliquePsi)
      ,CHFp_CONST_REAL(auObliqueK) )
{
 CH_TIMELEAF("FORT_SETOBLIQUE");
 FORTRAN_NAME( SETOBLIQUE ,setoblique )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auDensityL)
      ,CHFt_CONST_REAL(auPressureL)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auVelyL)
      ,CHFt_CONST_REAL(auAngleSW)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_REAL(auObliqueAe)
      ,CHFt_CONST_REAL(auObliqueAv)
      ,CHFt_CONST_REAL(auObliquePsi)
      ,CHFt_CONST_REAL(auObliqueK) );
}
#endif  // GUARDSETOBLIQUE 

#ifndef GUARDOBLIQUEINIT 
#define GUARDOBLIQUEINIT 
// Prototype for Fortran procedure OBLIQUEINIT ...
//
void FORTRAN_NAME( OBLIQUEINIT ,obliqueinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_OBLIQUEINIT FORTRAN_NAME( inlineOBLIQUEINIT, inlineOBLIQUEINIT)
#define FORTNT_OBLIQUEINIT FORTRAN_NAME( OBLIQUEINIT, obliqueinit)

inline void FORTRAN_NAME(inlineOBLIQUEINIT, inlineOBLIQUEINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OBLIQUEINIT");
 FORTRAN_NAME( OBLIQUEINIT ,obliqueinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOBLIQUEINIT 

#ifndef GUARDFLUXBCOS 
#define GUARDFLUXBCOS 
// Prototype for Fortran procedure FLUXBCOS ...
//
void FORTRAN_NAME( FLUXBCOS ,fluxbcos )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_FLUXBCOS FORTRAN_NAME( inlineFLUXBCOS, inlineFLUXBCOS)
#define FORTNT_FLUXBCOS FORTRAN_NAME( FLUXBCOS, fluxbcos)

inline void FORTRAN_NAME(inlineFLUXBCOS, inlineFLUXBCOS)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXBCOS");
 FORTRAN_NAME( FLUXBCOS ,fluxbcos )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXBCOS 

#ifndef GUARDOBLIQUEGS 
#define GUARDOBLIQUEGS 
// Prototype for Fortran procedure OBLIQUEGS ...
//
void FORTRAN_NAME( OBLIQUEGS ,obliquegs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_OBLIQUEGS FORTRAN_NAME( inlineOBLIQUEGS, inlineOBLIQUEGS)
#define FORTNT_OBLIQUEGS FORTRAN_NAME( OBLIQUEGS, obliquegs)

inline void FORTRAN_NAME(inlineOBLIQUEGS, inlineOBLIQUEGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OBLIQUEGS");
 FORTRAN_NAME( OBLIQUEGS ,obliquegs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOBLIQUEGS 

#ifndef GUARDOBLIQUESHOCK 
#define GUARDOBLIQUESHOCK 
// Prototype for Fortran procedure ObliqueShock ...
//
void FORTRAN_NAME( OBLIQUESHOCK ,obliqueshock )(
      CHFp_CONST_REAL(Rho1)
      ,CHFp_CONST_REAL(U1)
      ,CHFp_CONST_REAL(V1)
      ,CHFp_CONST_REAL(P1)
      ,CHFp_CONST_REAL(Gam)
      ,CHFp_CONST_REAL(ANGLE)
      ,CHFp_REAL(Rho2)
      ,CHFp_REAL(U2)
      ,CHFp_REAL(V2)
      ,CHFp_REAL(P2) );

#define FORT_OBLIQUESHOCK FORTRAN_NAME( inlineOBLIQUESHOCK, inlineOBLIQUESHOCK)
#define FORTNT_OBLIQUESHOCK FORTRAN_NAME( OBLIQUESHOCK, obliqueshock)

inline void FORTRAN_NAME(inlineOBLIQUESHOCK, inlineOBLIQUESHOCK)(
      CHFp_CONST_REAL(Rho1)
      ,CHFp_CONST_REAL(U1)
      ,CHFp_CONST_REAL(V1)
      ,CHFp_CONST_REAL(P1)
      ,CHFp_CONST_REAL(Gam)
      ,CHFp_CONST_REAL(ANGLE)
      ,CHFp_REAL(Rho2)
      ,CHFp_REAL(U2)
      ,CHFp_REAL(V2)
      ,CHFp_REAL(P2) )
{
 CH_TIMELEAF("FORT_OBLIQUESHOCK");
 FORTRAN_NAME( OBLIQUESHOCK ,obliqueshock )(
      CHFt_CONST_REAL(Rho1)
      ,CHFt_CONST_REAL(U1)
      ,CHFt_CONST_REAL(V1)
      ,CHFt_CONST_REAL(P1)
      ,CHFt_CONST_REAL(Gam)
      ,CHFt_CONST_REAL(ANGLE)
      ,CHFt_REAL(Rho2)
      ,CHFt_REAL(U2)
      ,CHFt_REAL(V2)
      ,CHFt_REAL(P2) );
}
#endif  // GUARDOBLIQUESHOCK 

}

#endif
