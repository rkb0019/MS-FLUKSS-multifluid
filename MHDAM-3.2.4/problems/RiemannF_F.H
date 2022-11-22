#ifndef _RIEMANNF_F_H_
#define _RIEMANNF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETRIEMANN 
#define GUARDSETRIEMANN 
// Prototype for Fortran procedure SETRIEMANN ...
//
void FORTRAN_NAME( SETRIEMANN ,setriemann )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auPressureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auXs) );

#define FORT_SETRIEMANN FORTRAN_NAME( inlineSETRIEMANN, inlineSETRIEMANN)
#define FORTNT_SETRIEMANN FORTRAN_NAME( SETRIEMANN, setriemann)

inline void FORTRAN_NAME(inlineSETRIEMANN, inlineSETRIEMANN)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auPressureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auXs) )
{
 CH_TIMELEAF("FORT_SETRIEMANN");
 FORTRAN_NAME( SETRIEMANN ,setriemann )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auDensityL)
      ,CHFt_CONST_REAL(auDensityR)
      ,CHFt_CONST_REAL(auPressureL)
      ,CHFt_CONST_REAL(auPressureR)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auVelxR)
      ,CHFt_CONST_REAL(auVelyL)
      ,CHFt_CONST_REAL(auVelyR)
      ,CHFt_CONST_REAL(auVelzL)
      ,CHFt_CONST_REAL(auVelzR)
      ,CHFt_CONST_REAL(auBxL)
      ,CHFt_CONST_REAL(auBxR)
      ,CHFt_CONST_REAL(auByL)
      ,CHFt_CONST_REAL(auByR)
      ,CHFt_CONST_REAL(auBzL)
      ,CHFt_CONST_REAL(auBzR)
      ,CHFt_CONST_REAL(auXs) );
}
#endif  // GUARDSETRIEMANN 

#ifndef GUARDRIEMANNINIT 
#define GUARDRIEMANNINIT 
// Prototype for Fortran procedure RIEMANNINIT ...
//
void FORTRAN_NAME( RIEMANNINIT ,riemanninit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNINIT FORTRAN_NAME( inlineRIEMANNINIT, inlineRIEMANNINIT)
#define FORTNT_RIEMANNINIT FORTRAN_NAME( RIEMANNINIT, riemanninit)

inline void FORTRAN_NAME(inlineRIEMANNINIT, inlineRIEMANNINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNINIT");
 FORTRAN_NAME( RIEMANNINIT ,riemanninit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNINIT 

#ifndef GUARDFLUXBC 
#define GUARDFLUXBC 
// Prototype for Fortran procedure FLUXBC ...
//
void FORTRAN_NAME( FLUXBC ,fluxbc )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_FLUXBC FORTRAN_NAME( inlineFLUXBC, inlineFLUXBC)
#define FORTNT_FLUXBC FORTRAN_NAME( FLUXBC, fluxbc)

inline void FORTRAN_NAME(inlineFLUXBC, inlineFLUXBC)(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXBC");
 FORTRAN_NAME( FLUXBC ,fluxbc )(
      CHFt_FRA(F)
      ,CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXBC 

#ifndef GUARDRIEMANNGS 
#define GUARDRIEMANNGS 
// Prototype for Fortran procedure RIEMANNGS ...
//
void FORTRAN_NAME( RIEMANNGS ,riemanngs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNGS FORTRAN_NAME( inlineRIEMANNGS, inlineRIEMANNGS)
#define FORTNT_RIEMANNGS FORTRAN_NAME( RIEMANNGS, riemanngs)

inline void FORTRAN_NAME(inlineRIEMANNGS, inlineRIEMANNGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNGS");
 FORTRAN_NAME( RIEMANNGS ,riemanngs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNGS 

}

#endif
