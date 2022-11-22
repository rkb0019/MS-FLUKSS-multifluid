#ifndef _HELIOSPHERICTURBF_F_H_
#define _HELIOSPHERICTURBF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETHELIOS_TM 
#define GUARDSETHELIOS_TM 
// Prototype for Fortran procedure SETHELIOS_TM ...
//
void FORTRAN_NAME( SETHELIOS_TM ,sethelios_tm )(
      CHFp_CONST_REAL(aulismZ2)
      ,CHFp_CONST_REAL(aulismLambda)
      ,CHFp_CONST_REAL(aulismSigmaC)
      ,CHFp_CONST_REAL(ausunZ2)
      ,CHFp_CONST_REAL(ausunLambda)
      ,CHFp_CONST_REAL(ausunSigmaC) );

#define FORT_SETHELIOS_TM FORTRAN_NAME( inlineSETHELIOS_TM, inlineSETHELIOS_TM)
#define FORTNT_SETHELIOS_TM FORTRAN_NAME( SETHELIOS_TM, sethelios_tm)

inline void FORTRAN_NAME(inlineSETHELIOS_TM, inlineSETHELIOS_TM)(
      CHFp_CONST_REAL(aulismZ2)
      ,CHFp_CONST_REAL(aulismLambda)
      ,CHFp_CONST_REAL(aulismSigmaC)
      ,CHFp_CONST_REAL(ausunZ2)
      ,CHFp_CONST_REAL(ausunLambda)
      ,CHFp_CONST_REAL(ausunSigmaC) )
{
 CH_TIMELEAF("FORT_SETHELIOS_TM");
 FORTRAN_NAME( SETHELIOS_TM ,sethelios_tm )(
      CHFt_CONST_REAL(aulismZ2)
      ,CHFt_CONST_REAL(aulismLambda)
      ,CHFt_CONST_REAL(aulismSigmaC)
      ,CHFt_CONST_REAL(ausunZ2)
      ,CHFt_CONST_REAL(ausunLambda)
      ,CHFt_CONST_REAL(ausunSigmaC) );
}
#endif  // GUARDSETHELIOS_TM 

#ifndef GUARDHELIOINIT_TM 
#define GUARDHELIOINIT_TM 
// Prototype for Fortran procedure HELIOINIT_TM ...
//
void FORTRAN_NAME( HELIOINIT_TM ,helioinit_tm )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) );

#define FORT_HELIOINIT_TM FORTRAN_NAME( inlineHELIOINIT_TM, inlineHELIOINIT_TM)
#define FORTNT_HELIOINIT_TM FORTRAN_NAME( HELIOINIT_TM, helioinit_tm)

inline void FORTRAN_NAME(inlineHELIOINIT_TM, inlineHELIOINIT_TM)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINIT_TM");
 FORTRAN_NAME( HELIOINIT_TM ,helioinit_tm )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINIT_TM 

#ifndef GUARDHELIOINITSPHERICAL_TM 
#define GUARDHELIOINITSPHERICAL_TM 
// Prototype for Fortran procedure HELIOINITSPHERICAL_TM ...
//
void FORTRAN_NAME( HELIOINITSPHERICAL_TM ,helioinitspherical_tm )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOINITSPHERICAL_TM FORTRAN_NAME( inlineHELIOINITSPHERICAL_TM, inlineHELIOINITSPHERICAL_TM)
#define FORTNT_HELIOINITSPHERICAL_TM FORTRAN_NAME( HELIOINITSPHERICAL_TM, helioinitspherical_tm)

inline void FORTRAN_NAME(inlineHELIOINITSPHERICAL_TM, inlineHELIOINITSPHERICAL_TM)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINITSPHERICAL_TM");
 FORTRAN_NAME( HELIOINITSPHERICAL_TM ,helioinitspherical_tm )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINITSPHERICAL_TM 

#ifndef GUARDHELIOREINIT_TM 
#define GUARDHELIOREINIT_TM 
// Prototype for Fortran procedure HELIOREINIT_TM ...
//
void FORTRAN_NAME( HELIOREINIT_TM ,helioreinit_tm )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINIT_TM FORTRAN_NAME( inlineHELIOREINIT_TM, inlineHELIOREINIT_TM)
#define FORTNT_HELIOREINIT_TM FORTRAN_NAME( HELIOREINIT_TM, helioreinit_tm)

inline void FORTRAN_NAME(inlineHELIOREINIT_TM, inlineHELIOREINIT_TM)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINIT_TM");
 FORTRAN_NAME( HELIOREINIT_TM ,helioreinit_tm )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINIT_TM 

#ifndef GUARDHELIOBC_TM 
#define GUARDHELIOBC_TM 
// Prototype for Fortran procedure HELIOBC_TM ...
//
void FORTRAN_NAME( HELIOBC_TM ,heliobc_tm )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) );

#define FORT_HELIOBC_TM FORTRAN_NAME( inlineHELIOBC_TM, inlineHELIOBC_TM)
#define FORTNT_HELIOBC_TM FORTRAN_NAME( HELIOBC_TM, heliobc_tm)

inline void FORTRAN_NAME(inlineHELIOBC_TM, inlineHELIOBC_TM)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOBC_TM");
 FORTRAN_NAME( HELIOBC_TM ,heliobc_tm )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOBC_TM 

#ifndef GUARDHELIOGS_TM 
#define GUARDHELIOGS_TM 
// Prototype for Fortran procedure HELIOGS_TM ...
//
void FORTRAN_NAME( HELIOGS_TM ,heliogs_tm )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) );

#define FORT_HELIOGS_TM FORTRAN_NAME( inlineHELIOGS_TM, inlineHELIOGS_TM)
#define FORTNT_HELIOGS_TM FORTRAN_NAME( HELIOGS_TM, heliogs_tm)

inline void FORTRAN_NAME(inlineHELIOGS_TM, inlineHELIOGS_TM)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGS_TM");
 FORTRAN_NAME( HELIOGS_TM ,heliogs_tm )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGS_TM 

#ifndef GUARDHELIOBCSPHERICAL_TM 
#define GUARDHELIOBCSPHERICAL_TM 
// Prototype for Fortran procedure HELIOBCSPHERICAL_TM ...
//
void FORTRAN_NAME( HELIOBCSPHERICAL_TM ,heliobcspherical_tm )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) );

#define FORT_HELIOBCSPHERICAL_TM FORTRAN_NAME( inlineHELIOBCSPHERICAL_TM, inlineHELIOBCSPHERICAL_TM)
#define FORTNT_HELIOBCSPHERICAL_TM FORTRAN_NAME( HELIOBCSPHERICAL_TM, heliobcspherical_tm)

inline void FORTRAN_NAME(inlineHELIOBCSPHERICAL_TM, inlineHELIOBCSPHERICAL_TM)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOBCSPHERICAL_TM");
 FORTRAN_NAME( HELIOBCSPHERICAL_TM ,heliobcspherical_tm )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOBCSPHERICAL_TM 

#ifndef GUARDHELIOGSSPHERICAL_TM 
#define GUARDHELIOGSSPHERICAL_TM 
// Prototype for Fortran procedure HELIOGSSPHERICAL_TM ...
//
void FORTRAN_NAME( HELIOGSSPHERICAL_TM ,heliogsspherical_tm )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSSPHERICAL_TM FORTRAN_NAME( inlineHELIOGSSPHERICAL_TM, inlineHELIOGSSPHERICAL_TM)
#define FORTNT_HELIOGSSPHERICAL_TM FORTRAN_NAME( HELIOGSSPHERICAL_TM, heliogsspherical_tm)

inline void FORTRAN_NAME(inlineHELIOGSSPHERICAL_TM, inlineHELIOGSSPHERICAL_TM)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSSPHERICAL_TM");
 FORTRAN_NAME( HELIOGSSPHERICAL_TM ,heliogsspherical_tm )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSSPHERICAL_TM 

#ifndef GUARDOMNIGSSPHERICAL_TM 
#define GUARDOMNIGSSPHERICAL_TM 
// Prototype for Fortran procedure OMNIGSSPHERICAL_TM ...
//
void FORTRAN_NAME( OMNIGSSPHERICAL_TM ,omnigsspherical_tm )(
      CHFp_FRA(W)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_Z2,ia_Z2hi0)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_BOX(box) );

#define FORT_OMNIGSSPHERICAL_TM FORTRAN_NAME( inlineOMNIGSSPHERICAL_TM, inlineOMNIGSSPHERICAL_TM)
#define FORTNT_OMNIGSSPHERICAL_TM FORTRAN_NAME( OMNIGSSPHERICAL_TM, omnigsspherical_tm)

inline void FORTRAN_NAME(inlineOMNIGSSPHERICAL_TM, inlineOMNIGSSPHERICAL_TM)(
      CHFp_FRA(W)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_Z2,ia_Z2hi0)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OMNIGSSPHERICAL_TM");
 FORTRAN_NAME( OMNIGSSPHERICAL_TM ,omnigsspherical_tm )(
      CHFt_FRA(W)
      ,CHFt_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFt_CONST_R1D(a_Z2,ia_Z2hi0)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_REAL(ltime)
      ,CHFt_CONST_REAL(startBC)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOMNIGSSPHERICAL_TM 

}

#endif
