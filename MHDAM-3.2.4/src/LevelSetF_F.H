#ifndef _LEVELSETF_F_H_
#define _LEVELSETF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDFLUXESHANCOCK_LS 
#define GUARDFLUXESHANCOCK_LS 
// Prototype for Fortran procedure FLUXESHANCOCK_LS ...
//
void FORTRAN_NAME( FLUXESHANCOCK_LS ,fluxeshancock_ls )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) );

#define FORT_FLUXESHANCOCK_LS FORTRAN_NAME( inlineFLUXESHANCOCK_LS, inlineFLUXESHANCOCK_LS)
#define FORTNT_FLUXESHANCOCK_LS FORTRAN_NAME( FLUXESHANCOCK_LS, fluxeshancock_ls)

inline void FORTRAN_NAME(inlineFLUXESHANCOCK_LS, inlineFLUXESHANCOCK_LS)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXESHANCOCK_LS");
 FORTRAN_NAME( FLUXESHANCOCK_LS ,fluxeshancock_ls )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXESHANCOCK_LS 

#ifndef GUARDLEVELSETFLUX 
#define GUARDLEVELSETFLUX 
// Prototype for Fortran procedure LEVELSETFLUX ...
//
void FORTRAN_NAME( LEVELSETFLUX ,levelsetflux )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) );

#define FORT_LEVELSETFLUX FORTRAN_NAME( inlineLEVELSETFLUX, inlineLEVELSETFLUX)
#define FORTNT_LEVELSETFLUX FORTRAN_NAME( LEVELSETFLUX, levelsetflux)

inline void FORTRAN_NAME(inlineLEVELSETFLUX, inlineLEVELSETFLUX)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LEVELSETFLUX");
 FORTRAN_NAME( LEVELSETFLUX ,levelsetflux )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLEVELSETFLUX 

#ifndef GUARDMINDT_LSCARTESIAN 
#define GUARDMINDT_LSCARTESIAN 
// Prototype for Fortran procedure MINDT_LSCARTESIAN ...
//
void FORTRAN_NAME( MINDT_LSCARTESIAN ,mindt_lscartesian )(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_MINDT_LSCARTESIAN FORTRAN_NAME( inlineMINDT_LSCARTESIAN, inlineMINDT_LSCARTESIAN)
#define FORTNT_MINDT_LSCARTESIAN FORTRAN_NAME( MINDT_LSCARTESIAN, mindt_lscartesian)

inline void FORTRAN_NAME(inlineMINDT_LSCARTESIAN, inlineMINDT_LSCARTESIAN)(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINDT_LSCARTESIAN");
 FORTRAN_NAME( MINDT_LSCARTESIAN ,mindt_lscartesian )(
      CHFt_FRA1(dt)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINDT_LSCARTESIAN 

#ifndef GUARDMINDT_LSSPHERICAL 
#define GUARDMINDT_LSSPHERICAL 
// Prototype for Fortran procedure MINDT_LSSPHERICAL ...
//
void FORTRAN_NAME( MINDT_LSSPHERICAL ,mindt_lsspherical )(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_MINDT_LSSPHERICAL FORTRAN_NAME( inlineMINDT_LSSPHERICAL, inlineMINDT_LSSPHERICAL)
#define FORTNT_MINDT_LSSPHERICAL FORTRAN_NAME( MINDT_LSSPHERICAL, mindt_lsspherical)

inline void FORTRAN_NAME(inlineMINDT_LSSPHERICAL, inlineMINDT_LSSPHERICAL)(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINDT_LSSPHERICAL");
 FORTRAN_NAME( MINDT_LSSPHERICAL ,mindt_lsspherical )(
      CHFt_FRA1(dt)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINDT_LSSPHERICAL 

#ifndef GUARDUPDATELEVELSET 
#define GUARDUPDATELEVELSET 
// Prototype for Fortran procedure UPDATELEVELSET ...
//
void FORTRAN_NAME( UPDATELEVELSET ,updatelevelset )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) );

#define FORT_UPDATELEVELSET FORTRAN_NAME( inlineUPDATELEVELSET, inlineUPDATELEVELSET)
#define FORTNT_UPDATELEVELSET FORTRAN_NAME( UPDATELEVELSET, updatelevelset)

inline void FORTRAN_NAME(inlineUPDATELEVELSET, inlineUPDATELEVELSET)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATELEVELSET");
 FORTRAN_NAME( UPDATELEVELSET ,updatelevelset )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA1(dx)
      ,CHFt_CONST_FRA1(g)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATELEVELSET 

#ifndef GUARDUPDATELEVELSET2 
#define GUARDUPDATELEVELSET2 
// Prototype for Fortran procedure UPDATELEVELSET2 ...
//
void FORTRAN_NAME( UPDATELEVELSET2 ,updatelevelset2 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) );

#define FORT_UPDATELEVELSET2 FORTRAN_NAME( inlineUPDATELEVELSET2, inlineUPDATELEVELSET2)
#define FORTNT_UPDATELEVELSET2 FORTRAN_NAME( UPDATELEVELSET2, updatelevelset2)

inline void FORTRAN_NAME(inlineUPDATELEVELSET2, inlineUPDATELEVELSET2)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATELEVELSET2");
 FORTRAN_NAME( UPDATELEVELSET2 ,updatelevelset2 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA1(dx)
      ,CHFt_CONST_FRA1(g)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATELEVELSET2 

#ifndef GUARDUPDATELEVELSET3 
#define GUARDUPDATELEVELSET3 
// Prototype for Fortran procedure UPDATELEVELSET3 ...
//
void FORTRAN_NAME( UPDATELEVELSET3 ,updatelevelset3 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(WMinus)
      ,CHFp_CONST_FRA(WPlus)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) );

#define FORT_UPDATELEVELSET3 FORTRAN_NAME( inlineUPDATELEVELSET3, inlineUPDATELEVELSET3)
#define FORTNT_UPDATELEVELSET3 FORTRAN_NAME( UPDATELEVELSET3, updatelevelset3)

inline void FORTRAN_NAME(inlineUPDATELEVELSET3, inlineUPDATELEVELSET3)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(WMinus)
      ,CHFp_CONST_FRA(WPlus)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATELEVELSET3");
 FORTRAN_NAME( UPDATELEVELSET3 ,updatelevelset3 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA(WMinus)
      ,CHFt_CONST_FRA(WPlus)
      ,CHFt_CONST_FRA1(dx)
      ,CHFt_CONST_FRA1(g)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATELEVELSET3 

#ifndef GUARDUPDATELEVELSETHANCOCK 
#define GUARDUPDATELEVELSETHANCOCK 
// Prototype for Fortran procedure UPDATELEVELSETHANCOCK ...
//
void FORTRAN_NAME( UPDATELEVELSETHANCOCK ,updatelevelsethancock )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) );

#define FORT_UPDATELEVELSETHANCOCK FORTRAN_NAME( inlineUPDATELEVELSETHANCOCK, inlineUPDATELEVELSETHANCOCK)
#define FORTNT_UPDATELEVELSETHANCOCK FORTRAN_NAME( UPDATELEVELSETHANCOCK, updatelevelsethancock)

inline void FORTRAN_NAME(inlineUPDATELEVELSETHANCOCK, inlineUPDATELEVELSETHANCOCK)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_FRA1(dx)
      ,CHFp_CONST_FRA1(g)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(ivf)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATELEVELSETHANCOCK");
 FORTRAN_NAME( UPDATELEVELSETHANCOCK ,updatelevelsethancock )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(FMinus)
      ,CHFt_CONST_FRA(FPlus)
      ,CHFt_CONST_FRA1(dx)
      ,CHFt_CONST_FRA1(g)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(ivf)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATELEVELSETHANCOCK 

#ifndef GUARDUPDATELEVELSETHANCOCK2 
#define GUARDUPDATELEVELSETHANCOCK2 
// Prototype for Fortran procedure UPDATELEVELSETHANCOCK2 ...
//
void FORTRAN_NAME( UPDATELEVELSETHANCOCK2 ,updatelevelsethancock2 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_UPDATELEVELSETHANCOCK2 FORTRAN_NAME( inlineUPDATELEVELSETHANCOCK2, inlineUPDATELEVELSETHANCOCK2)
#define FORTNT_UPDATELEVELSETHANCOCK2 FORTRAN_NAME( UPDATELEVELSETHANCOCK2, updatelevelsethancock2)

inline void FORTRAN_NAME(inlineUPDATELEVELSETHANCOCK2, inlineUPDATELEVELSETHANCOCK2)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATELEVELSETHANCOCK2");
 FORTRAN_NAME( UPDATELEVELSETHANCOCK2 ,updatelevelsethancock2 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(FMinus)
      ,CHFt_CONST_FRA(FPlus)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_FRA1(invvol)
      ,CHFt_CONST_FRA1(areas)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATELEVELSETHANCOCK2 

#ifndef GUARDPOSTPROCESSING_LS 
#define GUARDPOSTPROCESSING_LS 
// Prototype for Fortran procedure POSTPROCESSING_LS ...
//
void FORTRAN_NAME( POSTPROCESSING_LS ,postprocessing_ls )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(nls)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING_LS FORTRAN_NAME( inlinePOSTPROCESSING_LS, inlinePOSTPROCESSING_LS)
#define FORTNT_POSTPROCESSING_LS FORTRAN_NAME( POSTPROCESSING_LS, postprocessing_ls)

inline void FORTRAN_NAME(inlinePOSTPROCESSING_LS, inlinePOSTPROCESSING_LS)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_INT(ils)
      ,CHFp_CONST_INT(nls)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING_LS");
 FORTRAN_NAME( POSTPROCESSING_LS ,postprocessing_ls )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Uold)
      ,CHFt_CONST_INT(ils)
      ,CHFt_CONST_INT(nls)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING_LS 

#ifndef GUARDLEVELSETBCOR 
#define GUARDLEVELSETBCOR 
// Prototype for Fortran procedure LEVELSETBCOR ...
//
void FORTRAN_NAME( LEVELSETBCOR ,levelsetbcor )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) );

#define FORT_LEVELSETBCOR FORTRAN_NAME( inlineLEVELSETBCOR, inlineLEVELSETBCOR)
#define FORTNT_LEVELSETBCOR FORTRAN_NAME( LEVELSETBCOR, levelsetbcor)

inline void FORTRAN_NAME(inlineLEVELSETBCOR, inlineLEVELSETBCOR)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LEVELSETBCOR");
 FORTRAN_NAME( LEVELSETBCOR ,levelsetbcor )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(ils)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLEVELSETBCOR 

#ifndef GUARDLEVELSETBCOR2 
#define GUARDLEVELSETBCOR2 
// Prototype for Fortran procedure LEVELSETBCOR2 ...
//
void FORTRAN_NAME( LEVELSETBCOR2 ,levelsetbcor2 )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA1(HCS)
      ,CHFp_BOX(box) );

#define FORT_LEVELSETBCOR2 FORTRAN_NAME( inlineLEVELSETBCOR2, inlineLEVELSETBCOR2)
#define FORTNT_LEVELSETBCOR2 FORTRAN_NAME( LEVELSETBCOR2, levelsetbcor2)

inline void FORTRAN_NAME(inlineLEVELSETBCOR2, inlineLEVELSETBCOR2)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA1(HCS)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LEVELSETBCOR2");
 FORTRAN_NAME( LEVELSETBCOR2 ,levelsetbcor2 )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA1(HCS)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLEVELSETBCOR2 

#ifndef GUARDLSREINITIALIZE 
#define GUARDLSREINITIALIZE 
// Prototype for Fortran procedure LSREINITIALIZE ...
//
void FORTRAN_NAME( LSREINITIALIZE ,lsreinitialize )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) );

#define FORT_LSREINITIALIZE FORTRAN_NAME( inlineLSREINITIALIZE, inlineLSREINITIALIZE)
#define FORTNT_LSREINITIALIZE FORTRAN_NAME( LSREINITIALIZE, lsreinitialize)

inline void FORTRAN_NAME(inlineLSREINITIALIZE, inlineLSREINITIALIZE)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(ils)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LSREINITIALIZE");
 FORTRAN_NAME( LSREINITIALIZE ,lsreinitialize )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(ils)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLSREINITIALIZE 

}

#endif
