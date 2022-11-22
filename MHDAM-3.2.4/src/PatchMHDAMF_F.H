#ifndef _PATCHMHDAMF_F_H_
#define _PATCHMHDAMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDGETGRAD 
#define GUARDGETGRAD 
// Prototype for Fortran procedure GETGRAD ...
//
void FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETGRAD FORTRAN_NAME( inlineGETGRAD, inlineGETGRAD)
#define FORTNT_GETGRAD FORTRAN_NAME( GETGRAD, getgrad)

inline void FORTRAN_NAME(inlineGETGRAD, inlineGETGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETGRAD");
 FORTRAN_NAME( GETGRAD ,getgrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETGRAD 

#ifndef GUARDMAGNITUDE 
#define GUARDMAGNITUDE 
// Prototype for Fortran procedure MAGNITUDE ...
//
void FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) );

#define FORT_MAGNITUDE FORTRAN_NAME( inlineMAGNITUDE, inlineMAGNITUDE)
#define FORTNT_MAGNITUDE FORTRAN_NAME( MAGNITUDE, magnitude)

inline void FORTRAN_NAME(inlineMAGNITUDE, inlineMAGNITUDE)(
      CHFp_FRA1(magdata)
      ,CHFp_CONST_FRA(data)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAGNITUDE");
 FORTRAN_NAME( MAGNITUDE ,magnitude )(
      CHFt_FRA1(magdata)
      ,CHFt_CONST_FRA(data)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAGNITUDE 

#ifndef GUARDSECONDSLOPEDIFFS 
#define GUARDSECONDSLOPEDIFFS 
// Prototype for Fortran procedure SECONDSLOPEDIFFS ...
//
void FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_SECONDSLOPEDIFFS FORTRAN_NAME( inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)
#define FORTNT_SECONDSLOPEDIFFS FORTRAN_NAME( SECONDSLOPEDIFFS, secondslopediffs)

inline void FORTRAN_NAME(inlineSECONDSLOPEDIFFS, inlineSECONDSLOPEDIFFS)(
      CHFp_FRA(deltaWC)
      ,CHFp_FRA(deltaWL)
      ,CHFp_FRA(deltaWR)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_SECONDSLOPEDIFFS");
 FORTRAN_NAME( SECONDSLOPEDIFFS ,secondslopediffs )(
      CHFt_FRA(deltaWC)
      ,CHFt_FRA(deltaWL)
      ,CHFt_FRA(deltaWR)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDSECONDSLOPEDIFFS 

#ifndef GUARDFOURTHSLOPEDIFFS 
#define GUARDFOURTHSLOPEDIFFS 
// Prototype for Fortran procedure FOURTHSLOPEDIFFS ...
//
void FORTRAN_NAME( FOURTHSLOPEDIFFS ,fourthslopediffs )(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_FOURTHSLOPEDIFFS FORTRAN_NAME( inlineFOURTHSLOPEDIFFS, inlineFOURTHSLOPEDIFFS)
#define FORTNT_FOURTHSLOPEDIFFS FORTRAN_NAME( FOURTHSLOPEDIFFS, fourthslopediffs)

inline void FORTRAN_NAME(inlineFOURTHSLOPEDIFFS, inlineFOURTHSLOPEDIFFS)(
      CHFp_FRA(delta4WC)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(delta2W)
      ,CHFp_CONST_INT(numSlopes)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_FOURTHSLOPEDIFFS");
 FORTRAN_NAME( FOURTHSLOPEDIFFS ,fourthslopediffs )(
      CHFt_FRA(delta4WC)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(delta2W)
      ,CHFt_CONST_INT(numSlopes)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDFOURTHSLOPEDIFFS 

#ifndef GUARDDIVUEDGE 
#define GUARDDIVUEDGE 
// Prototype for Fortran procedure DIVUEDGE ...
//
void FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(uNorm)
      ,CHFp_CONST_FRA(duTan)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_DIVUEDGE FORTRAN_NAME( inlineDIVUEDGE, inlineDIVUEDGE)
#define FORTNT_DIVUEDGE FORTRAN_NAME( DIVUEDGE, divuedge)

inline void FORTRAN_NAME(inlineDIVUEDGE, inlineDIVUEDGE)(
      CHFp_FRA1(divu)
      ,CHFp_CONST_FRA1(uNorm)
      ,CHFp_CONST_FRA(duTan)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_DIVUEDGE");
 FORTRAN_NAME( DIVUEDGE ,divuedge )(
      CHFt_FRA1(divu)
      ,CHFt_CONST_FRA1(uNorm)
      ,CHFt_CONST_FRA(duTan)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDDIVUEDGE 

#ifndef GUARDARTVISC 
#define GUARDARTVISC 
// Prototype for Fortran procedure ARTVISC ...
//
void FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_ARTVISC FORTRAN_NAME( inlineARTVISC, inlineARTVISC)
#define FORTNT_ARTVISC FORTRAN_NAME( ARTVISC, artvisc)

inline void FORTRAN_NAME(inlineARTVISC, inlineARTVISC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_FRA1(divu)
      ,CHFp_CONST_REAL(coeff)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ARTVISC");
 FORTRAN_NAME( ARTVISC ,artvisc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_FRA1(divu)
      ,CHFt_CONST_REAL(coeff)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDARTVISC 

#ifndef GUARDAVERAGESTATE 
#define GUARDAVERAGESTATE 
// Prototype for Fortran procedure AVERAGESTATE ...
//
void FORTRAN_NAME( AVERAGESTATE ,averagestate )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UOLD)
      ,CHFp_CONST_FRA(UPr)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) );

#define FORT_AVERAGESTATE FORTRAN_NAME( inlineAVERAGESTATE, inlineAVERAGESTATE)
#define FORTNT_AVERAGESTATE FORTRAN_NAME( AVERAGESTATE, averagestate)

inline void FORTRAN_NAME(inlineAVERAGESTATE, inlineAVERAGESTATE)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UOLD)
      ,CHFp_CONST_FRA(UPr)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_AVERAGESTATE");
 FORTRAN_NAME( AVERAGESTATE ,averagestate )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(UOLD)
      ,CHFt_CONST_FRA(UPr)
      ,CHFt_CONST_INT(iBGN)
      ,CHFt_CONST_INT(iEND)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAVERAGESTATE 

#ifndef GUARDGETRELGRAD 
#define GUARDGETRELGRAD 
// Prototype for Fortran procedure GETRELGRAD ...
//
void FORTRAN_NAME( GETRELGRAD ,getrelgrad )(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_GETRELGRAD FORTRAN_NAME( inlineGETRELGRAD, inlineGETRELGRAD)
#define FORTNT_GETRELGRAD FORTRAN_NAME( GETRELGRAD, getrelgrad)

inline void FORTRAN_NAME(inlineGETRELGRAD, inlineGETRELGRAD)(
      CHFp_FRA1(du)
      ,CHFp_CONST_FRA1(u)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_GETRELGRAD");
 FORTRAN_NAME( GETRELGRAD ,getrelgrad )(
      CHFt_FRA1(du)
      ,CHFt_CONST_FRA1(u)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDGETRELGRAD 

#ifndef GUARDGETVECTMAGNITUDE 
#define GUARDGETVECTMAGNITUDE 
// Prototype for Fortran procedure GETVECTMAGNITUDE ...
//
void FORTRAN_NAME( GETVECTMAGNITUDE ,getvectmagnitude )(
      CHFp_FRA1(VMag)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) );

#define FORT_GETVECTMAGNITUDE FORTRAN_NAME( inlineGETVECTMAGNITUDE, inlineGETVECTMAGNITUDE)
#define FORTNT_GETVECTMAGNITUDE FORTRAN_NAME( GETVECTMAGNITUDE, getvectmagnitude)

inline void FORTRAN_NAME(inlineGETVECTMAGNITUDE, inlineGETVECTMAGNITUDE)(
      CHFp_FRA1(VMag)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GETVECTMAGNITUDE");
 FORTRAN_NAME( GETVECTMAGNITUDE ,getvectmagnitude )(
      CHFt_FRA1(VMag)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iBGN)
      ,CHFt_CONST_INT(iEND)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGETVECTMAGNITUDE 

#ifndef GUARDGETRELGRAD_SPHERICAL 
#define GUARDGETRELGRAD_SPHERICAL 
// Prototype for Fortran procedure GETRELGRAD_SPHERICAL ...
//
void FORTRAN_NAME( GETRELGRAD_SPHERICAL ,getrelgrad_spherical )(
      CHFp_FRA1(modgrad)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_GETRELGRAD_SPHERICAL FORTRAN_NAME( inlineGETRELGRAD_SPHERICAL, inlineGETRELGRAD_SPHERICAL)
#define FORTNT_GETRELGRAD_SPHERICAL FORTRAN_NAME( GETRELGRAD_SPHERICAL, getrelgrad_spherical)

inline void FORTRAN_NAME(inlineGETRELGRAD_SPHERICAL, inlineGETRELGRAD_SPHERICAL)(
      CHFp_FRA1(modgrad)
      ,CHFp_CONST_FRA1(var)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_GETRELGRAD_SPHERICAL");
 FORTRAN_NAME( GETRELGRAD_SPHERICAL ,getrelgrad_spherical )(
      CHFt_FRA1(modgrad)
      ,CHFt_CONST_FRA1(var)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDGETRELGRAD_SPHERICAL 

#ifndef GUARDTRANSFORMFLUX_CYL 
#define GUARDTRANSFORMFLUX_CYL 
// Prototype for Fortran procedure TRANSFORMFLUX_CYL ...
//
void FORTRAN_NAME( TRANSFORMFLUX_CYL ,transformflux_cyl )(
      CHFp_FRA(F)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_TRANSFORMFLUX_CYL FORTRAN_NAME( inlineTRANSFORMFLUX_CYL, inlineTRANSFORMFLUX_CYL)
#define FORTNT_TRANSFORMFLUX_CYL FORTRAN_NAME( TRANSFORMFLUX_CYL, transformflux_cyl)

inline void FORTRAN_NAME(inlineTRANSFORMFLUX_CYL, inlineTRANSFORMFLUX_CYL)(
      CHFp_FRA(F)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_TRANSFORMFLUX_CYL");
 FORTRAN_NAME( TRANSFORMFLUX_CYL ,transformflux_cyl )(
      CHFt_FRA(F)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTRANSFORMFLUX_CYL 

#ifndef GUARDTRANSFORMUCYL_FORWARD 
#define GUARDTRANSFORMUCYL_FORWARD 
// Prototype for Fortran procedure TRANSFORMUCYL_FORWARD ...
//
void FORTRAN_NAME( TRANSFORMUCYL_FORWARD ,transformucyl_forward )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_TRANSFORMUCYL_FORWARD FORTRAN_NAME( inlineTRANSFORMUCYL_FORWARD, inlineTRANSFORMUCYL_FORWARD)
#define FORTNT_TRANSFORMUCYL_FORWARD FORTRAN_NAME( TRANSFORMUCYL_FORWARD, transformucyl_forward)

inline void FORTRAN_NAME(inlineTRANSFORMUCYL_FORWARD, inlineTRANSFORMUCYL_FORWARD)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_TRANSFORMUCYL_FORWARD");
 FORTRAN_NAME( TRANSFORMUCYL_FORWARD ,transformucyl_forward )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTRANSFORMUCYL_FORWARD 

#ifndef GUARDTRANSFORMUCYL_BACKWARD 
#define GUARDTRANSFORMUCYL_BACKWARD 
// Prototype for Fortran procedure TRANSFORMUCYL_BACKWARD ...
//
void FORTRAN_NAME( TRANSFORMUCYL_BACKWARD ,transformucyl_backward )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_TRANSFORMUCYL_BACKWARD FORTRAN_NAME( inlineTRANSFORMUCYL_BACKWARD, inlineTRANSFORMUCYL_BACKWARD)
#define FORTNT_TRANSFORMUCYL_BACKWARD FORTRAN_NAME( TRANSFORMUCYL_BACKWARD, transformucyl_backward)

inline void FORTRAN_NAME(inlineTRANSFORMUCYL_BACKWARD, inlineTRANSFORMUCYL_BACKWARD)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_TRANSFORMUCYL_BACKWARD");
 FORTRAN_NAME( TRANSFORMUCYL_BACKWARD ,transformucyl_backward )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTRANSFORMUCYL_BACKWARD 

#ifndef GUARDUPDATECONS 
#define GUARDUPDATECONS 
// Prototype for Fortran procedure UPDATECONS ...
//
void FORTRAN_NAME( UPDATECONS ,updatecons )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_UPDATECONS FORTRAN_NAME( inlineUPDATECONS, inlineUPDATECONS)
#define FORTNT_UPDATECONS FORTRAN_NAME( UPDATECONS, updatecons)

inline void FORTRAN_NAME(inlineUPDATECONS, inlineUPDATECONS)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATECONS");
 FORTRAN_NAME( UPDATECONS ,updatecons )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_INT(iBgn)
      ,CHFt_CONST_INT(iEnd)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_FRA1(invvol)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATECONS 

#ifndef GUARDUPDATECONSHANCOCK 
#define GUARDUPDATECONSHANCOCK 
// Prototype for Fortran procedure UPDATECONSHANCOCK ...
//
void FORTRAN_NAME( UPDATECONSHANCOCK ,updateconshancock )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_UPDATECONSHANCOCK FORTRAN_NAME( inlineUPDATECONSHANCOCK, inlineUPDATECONSHANCOCK)
#define FORTNT_UPDATECONSHANCOCK FORTRAN_NAME( UPDATECONSHANCOCK, updateconshancock)

inline void FORTRAN_NAME(inlineUPDATECONSHANCOCK, inlineUPDATECONSHANCOCK)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(FMinus)
      ,CHFp_CONST_FRA(FPlus)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATECONSHANCOCK");
 FORTRAN_NAME( UPDATECONSHANCOCK ,updateconshancock )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(FMinus)
      ,CHFt_CONST_FRA(FPlus)
      ,CHFt_CONST_INT(iBgn)
      ,CHFt_CONST_INT(iEnd)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_FRA1(invvol)
      ,CHFt_CONST_FRA1(areas)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATECONSHANCOCK 

#ifndef GUARDUPDATEVALUESONFACES 
#define GUARDUPDATEVALUESONFACES 
// Prototype for Fortran procedure UPDATEVALUESONFACES ...
//
void FORTRAN_NAME( UPDATEVALUESONFACES ,updatevaluesonfaces )(
      CHFp_FRA(WMinus)
      ,CHFp_FRA(WPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(WOld)
      ,CHFp_BOX(box) );

#define FORT_UPDATEVALUESONFACES FORTRAN_NAME( inlineUPDATEVALUESONFACES, inlineUPDATEVALUESONFACES)
#define FORTNT_UPDATEVALUESONFACES FORTRAN_NAME( UPDATEVALUESONFACES, updatevaluesonfaces)

inline void FORTRAN_NAME(inlineUPDATEVALUESONFACES, inlineUPDATEVALUESONFACES)(
      CHFp_FRA(WMinus)
      ,CHFp_FRA(WPlus)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(WOld)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPDATEVALUESONFACES");
 FORTRAN_NAME( UPDATEVALUESONFACES ,updatevaluesonfaces )(
      CHFt_FRA(WMinus)
      ,CHFt_FRA(WPlus)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(WOld)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPDATEVALUESONFACES 

#ifndef GUARDADDSOURCES 
#define GUARDADDSOURCES 
// Prototype for Fortran procedure ADDSOURCES ...
//
void FORTRAN_NAME( ADDSOURCES ,addsources )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) );

#define FORT_ADDSOURCES FORTRAN_NAME( inlineADDSOURCES, inlineADDSOURCES)
#define FORTNT_ADDSOURCES FORTRAN_NAME( ADDSOURCES, addsources)

inline void FORTRAN_NAME(inlineADDSOURCES, inlineADDSOURCES)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADDSOURCES");
 FORTRAN_NAME( ADDSOURCES ,addsources )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(S)
      ,CHFt_CONST_INT(iBGN)
      ,CHFt_CONST_INT(iEND)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADDSOURCES 

}

#endif
