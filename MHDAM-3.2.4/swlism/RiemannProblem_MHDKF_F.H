#ifndef _RIEMANNPROBLEM_MHDKF_F_H_
#define _RIEMANNPROBLEM_MHDKF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETRIEMANN_MHDK 
#define GUARDSETRIEMANN_MHDK 
// Prototype for Fortran procedure SETRIEMANN_MHDK ...
//
void FORTRAN_NAME( SETRIEMANN_MHDK ,setriemann_mhdk )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auNumdenL)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auNumdenR)
      ,CHFp_CONST_REAL(auTemperatureR)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByR) );

#define FORT_SETRIEMANN_MHDK FORTRAN_NAME( inlineSETRIEMANN_MHDK, inlineSETRIEMANN_MHDK)
#define FORTNT_SETRIEMANN_MHDK FORTRAN_NAME( SETRIEMANN_MHDK, setriemann_mhdk)

inline void FORTRAN_NAME(inlineSETRIEMANN_MHDK, inlineSETRIEMANN_MHDK)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auNumdenL)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auNumdenR)
      ,CHFp_CONST_REAL(auTemperatureR)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByR) )
{
 CH_TIMELEAF("FORT_SETRIEMANN_MHDK");
 FORTRAN_NAME( SETRIEMANN_MHDK ,setriemann_mhdk )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_REAL(auNetnumL)
      ,CHFt_CONST_REAL(auNumdenL)
      ,CHFt_CONST_REAL(auTemperatureL)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auVelyL)
      ,CHFt_CONST_REAL(auBxL)
      ,CHFt_CONST_REAL(auByL)
      ,CHFt_CONST_REAL(auNumdenR)
      ,CHFt_CONST_REAL(auTemperatureR)
      ,CHFt_CONST_REAL(auVelxR)
      ,CHFt_CONST_REAL(auVelyR)
      ,CHFt_CONST_REAL(auBxR)
      ,CHFt_CONST_REAL(auByR) );
}
#endif  // GUARDSETRIEMANN_MHDK 

#ifndef GUARDRIEMANNINIT_MHDK 
#define GUARDRIEMANNINIT_MHDK 
// Prototype for Fortran procedure RIEMANNINIT_MHDK ...
//
void FORTRAN_NAME( RIEMANNINIT_MHDK ,riemanninit_mhdk )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNINIT_MHDK FORTRAN_NAME( inlineRIEMANNINIT_MHDK, inlineRIEMANNINIT_MHDK)
#define FORTNT_RIEMANNINIT_MHDK FORTRAN_NAME( RIEMANNINIT_MHDK, riemanninit_mhdk)

inline void FORTRAN_NAME(inlineRIEMANNINIT_MHDK, inlineRIEMANNINIT_MHDK)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNINIT_MHDK");
 FORTRAN_NAME( RIEMANNINIT_MHDK ,riemanninit_mhdk )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNINIT_MHDK 

#ifndef GUARDFLUXBC_MHDK 
#define GUARDFLUXBC_MHDK 
// Prototype for Fortran procedure FLUXBC_MHDK ...
//
void FORTRAN_NAME( FLUXBC_MHDK ,fluxbc_mhdk )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_FLUXBC_MHDK FORTRAN_NAME( inlineFLUXBC_MHDK, inlineFLUXBC_MHDK)
#define FORTNT_FLUXBC_MHDK FORTRAN_NAME( FLUXBC_MHDK, fluxbc_mhdk)

inline void FORTRAN_NAME(inlineFLUXBC_MHDK, inlineFLUXBC_MHDK)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXBC_MHDK");
 FORTRAN_NAME( FLUXBC_MHDK ,fluxbc_mhdk )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXBC_MHDK 

}

#endif
