#ifndef _HELIOSPHERICPIF_F_H_
#define _HELIOSPHERICPIF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETHELIOS_PI 
#define GUARDSETHELIOS_PI 
// Prototype for Fortran procedure SETHELIOS_PI ...
//
void FORTRAN_NAME( SETHELIOS_PI ,sethelios_pi )(
      CHFp_CONST_REAL(aulismNPI)
      ,CHFp_CONST_REAL(aulismTPI)
      ,CHFp_CONST_REAL(ausunNPI)
      ,CHFp_CONST_REAL(ausunTPI) );

#define FORT_SETHELIOS_PI FORTRAN_NAME( inlineSETHELIOS_PI, inlineSETHELIOS_PI)
#define FORTNT_SETHELIOS_PI FORTRAN_NAME( SETHELIOS_PI, sethelios_pi)

inline void FORTRAN_NAME(inlineSETHELIOS_PI, inlineSETHELIOS_PI)(
      CHFp_CONST_REAL(aulismNPI)
      ,CHFp_CONST_REAL(aulismTPI)
      ,CHFp_CONST_REAL(ausunNPI)
      ,CHFp_CONST_REAL(ausunTPI) )
{
 CH_TIMELEAF("FORT_SETHELIOS_PI");
 FORTRAN_NAME( SETHELIOS_PI ,sethelios_pi )(
      CHFt_CONST_REAL(aulismNPI)
      ,CHFt_CONST_REAL(aulismTPI)
      ,CHFt_CONST_REAL(ausunNPI)
      ,CHFt_CONST_REAL(ausunTPI) );
}
#endif  // GUARDSETHELIOS_PI 

#ifndef GUARDHELIOINIT_PI 
#define GUARDHELIOINIT_PI 
// Prototype for Fortran procedure HELIOINIT_PI ...
//
void FORTRAN_NAME( HELIOINIT_PI ,helioinit_pi )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_HELIOINIT_PI FORTRAN_NAME( inlineHELIOINIT_PI, inlineHELIOINIT_PI)
#define FORTNT_HELIOINIT_PI FORTRAN_NAME( HELIOINIT_PI, helioinit_pi)

inline void FORTRAN_NAME(inlineHELIOINIT_PI, inlineHELIOINIT_PI)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINIT_PI");
 FORTRAN_NAME( HELIOINIT_PI ,helioinit_pi )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINIT_PI 

#ifndef GUARDHELIOINITSPHERICAL_PI 
#define GUARDHELIOINITSPHERICAL_PI 
// Prototype for Fortran procedure HELIOINITSPHERICAL_PI ...
//
void FORTRAN_NAME( HELIOINITSPHERICAL_PI ,helioinitspherical_pi )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOINITSPHERICAL_PI FORTRAN_NAME( inlineHELIOINITSPHERICAL_PI, inlineHELIOINITSPHERICAL_PI)
#define FORTNT_HELIOINITSPHERICAL_PI FORTRAN_NAME( HELIOINITSPHERICAL_PI, helioinitspherical_pi)

inline void FORTRAN_NAME(inlineHELIOINITSPHERICAL_PI, inlineHELIOINITSPHERICAL_PI)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINITSPHERICAL_PI");
 FORTRAN_NAME( HELIOINITSPHERICAL_PI ,helioinitspherical_pi )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINITSPHERICAL_PI 

#ifndef GUARDHELIOREINIT_PI 
#define GUARDHELIOREINIT_PI 
// Prototype for Fortran procedure HELIOREINIT_PI ...
//
void FORTRAN_NAME( HELIOREINIT_PI ,helioreinit_pi )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINIT_PI FORTRAN_NAME( inlineHELIOREINIT_PI, inlineHELIOREINIT_PI)
#define FORTNT_HELIOREINIT_PI FORTRAN_NAME( HELIOREINIT_PI, helioreinit_pi)

inline void FORTRAN_NAME(inlineHELIOREINIT_PI, inlineHELIOREINIT_PI)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINIT_PI");
 FORTRAN_NAME( HELIOREINIT_PI ,helioreinit_pi )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINIT_PI 


#ifndef GUARDHELIOREINITSPHERICAL_PI 
#define GUARDHELIOREINITSPHERICAL_PI
// Prototype for Fortran procedure HELIOREINITSPHERICAL_PI ...
//
void FORTRAN_NAME( HELIOREINITSPHERICAL_PI ,helioreinitspherical_pi )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINITSPHERICAL_PI FORTRAN_NAME( inlineHELIOREINITSPHERICAL_PI, inlineHELIOREINITSPHERICAL_PI)
#define FORTNT_HELIOREINITSPHERICAL_PI FORTRAN_NAME( HELIOREINITSPHERICAL_PI, helioreinitspherical_pi)

inline void FORTRAN_NAME(inlineHELIOREINITSPHERICAL_PI, inlineHELIOREINITSPHERICAL_PI)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINITSPHERICAL_PI");
 FORTRAN_NAME( HELIOREINITSPHERICAL_PI ,helioreinitspherical_pi )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_REAL(Rinit)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINITSPHERICAL_PI 

#ifndef GUARDHELIOBC_PI 
#define GUARDHELIOBC_PI 
// Prototype for Fortran procedure HELIOBC_PI ...
//
void FORTRAN_NAME( HELIOBC_PI ,heliobc_pi )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_HELIOBC_PI FORTRAN_NAME( inlineHELIOBC_PI, inlineHELIOBC_PI)
#define FORTNT_HELIOBC_PI FORTRAN_NAME( HELIOBC_PI, heliobc_pi)

inline void FORTRAN_NAME(inlineHELIOBC_PI, inlineHELIOBC_PI)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOBC_PI");
 FORTRAN_NAME( HELIOBC_PI ,heliobc_pi )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOBC_PI 

#ifndef GUARDHELIOGS_PI 
#define GUARDHELIOGS_PI 
// Prototype for Fortran procedure HELIOGS_PI ...
//
void FORTRAN_NAME( HELIOGS_PI ,heliogs_pi )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_HELIOGS_PI FORTRAN_NAME( inlineHELIOGS_PI, inlineHELIOGS_PI)
#define FORTNT_HELIOGS_PI FORTRAN_NAME( HELIOGS_PI, heliogs_pi)

inline void FORTRAN_NAME(inlineHELIOGS_PI, inlineHELIOGS_PI)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGS_PI");
 FORTRAN_NAME( HELIOGS_PI ,heliogs_pi )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGS_PI 

#ifndef GUARDHELIOBCSPHERICAL_PI 
#define GUARDHELIOBCSPHERICAL_PI 
// Prototype for Fortran procedure HELIOBCSPHERICAL_PI ...
//
void FORTRAN_NAME( HELIOBCSPHERICAL_PI ,heliobcspherical_pi )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_HELIOBCSPHERICAL_PI FORTRAN_NAME( inlineHELIOBCSPHERICAL_PI, inlineHELIOBCSPHERICAL_PI)
#define FORTNT_HELIOBCSPHERICAL_PI FORTRAN_NAME( HELIOBCSPHERICAL_PI, heliobcspherical_pi)

inline void FORTRAN_NAME(inlineHELIOBCSPHERICAL_PI, inlineHELIOBCSPHERICAL_PI)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOBCSPHERICAL_PI");
 FORTRAN_NAME( HELIOBCSPHERICAL_PI ,heliobcspherical_pi )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOBCSPHERICAL_PI 

#ifndef GUARDHELIOGSSPHERICAL_PI 
#define GUARDHELIOGSSPHERICAL_PI 
// Prototype for Fortran procedure HELIOGSSPHERICAL_PI ...
//
void FORTRAN_NAME( HELIOGSSPHERICAL_PI ,heliogsspherical_pi )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSSPHERICAL_PI FORTRAN_NAME( inlineHELIOGSSPHERICAL_PI, inlineHELIOGSSPHERICAL_PI)
#define FORTNT_HELIOGSSPHERICAL_PI FORTRAN_NAME( HELIOGSSPHERICAL_PI, heliogsspherical_pi)

inline void FORTRAN_NAME(inlineHELIOGSSPHERICAL_PI, inlineHELIOGSSPHERICAL_PI)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSSPHERICAL_PI");
 FORTRAN_NAME( HELIOGSSPHERICAL_PI ,heliogsspherical_pi )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSSPHERICAL_PI 

#ifndef GUARDOMNIGSSPHERICAL_PI 
#define GUARDOMNIGSSPHERICAL_PI 
// Prototype for Fortran procedure OMNIGSSPHERICAL_PI ...
//
void FORTRAN_NAME( OMNIGSSPHERICAL_PI ,omnigsspherical_pi )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_OMNIGSSPHERICAL_PI FORTRAN_NAME( inlineOMNIGSSPHERICAL_PI, inlineOMNIGSSPHERICAL_PI)
#define FORTNT_OMNIGSSPHERICAL_PI FORTRAN_NAME( OMNIGSSPHERICAL_PI, omnigsspherical_pi)

inline void FORTRAN_NAME(inlineOMNIGSSPHERICAL_PI, inlineOMNIGSSPHERICAL_PI)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OMNIGSSPHERICAL_PI");
 FORTRAN_NAME( OMNIGSSPHERICAL_PI ,omnigsspherical_pi )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOMNIGSSPHERICAL_PI 

}

#endif
