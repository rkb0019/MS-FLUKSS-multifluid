#ifndef _PATCHIDEALMHDF_F_H_
#define _PATCHIDEALMHDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAXWAVESPEED 
#define GUARDMAXWAVESPEED 
// Prototype for Fortran procedure MAXWAVESPEED ...
//
void FORTRAN_NAME( MAXWAVESPEED ,maxwavespeed )(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_BOX(box) );

#define FORT_MAXWAVESPEED FORTRAN_NAME( inlineMAXWAVESPEED, inlineMAXWAVESPEED)
#define FORTNT_MAXWAVESPEED FORTRAN_NAME( MAXWAVESPEED, maxwavespeed)

inline void FORTRAN_NAME(inlineMAXWAVESPEED, inlineMAXWAVESPEED)(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXWAVESPEED");
 FORTRAN_NAME( MAXWAVESPEED ,maxwavespeed )(
      CHFt_FRA1(cmax)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXWAVESPEED 

#ifndef GUARDMINDT_SPHERICAL 
#define GUARDMINDT_SPHERICAL 
// Prototype for Fortran procedure MINDT_SPHERICAL ...
//
void FORTRAN_NAME( MINDT_SPHERICAL ,mindt_spherical )(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_MINDT_SPHERICAL FORTRAN_NAME( inlineMINDT_SPHERICAL, inlineMINDT_SPHERICAL)
#define FORTNT_MINDT_SPHERICAL FORTRAN_NAME( MINDT_SPHERICAL, mindt_spherical)

inline void FORTRAN_NAME(inlineMINDT_SPHERICAL, inlineMINDT_SPHERICAL)(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINDT_SPHERICAL");
 FORTRAN_NAME( MINDT_SPHERICAL ,mindt_spherical )(
      CHFt_FRA1(dt)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINDT_SPHERICAL 

#ifndef GUARDCONSTOPRIM 
#define GUARDCONSTOPRIM 
// Prototype for Fortran procedure CONSTOPRIM ...
//
void FORTRAN_NAME( CONSTOPRIM ,constoprim )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) );

#define FORT_CONSTOPRIM FORTRAN_NAME( inlineCONSTOPRIM, inlineCONSTOPRIM)
#define FORTNT_CONSTOPRIM FORTRAN_NAME( CONSTOPRIM, constoprim)

inline void FORTRAN_NAME(inlineCONSTOPRIM, inlineCONSTOPRIM)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONSTOPRIM");
 FORTRAN_NAME( CONSTOPRIM ,constoprim )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONSTOPRIM 

#ifndef GUARDPRIMTOCONS 
#define GUARDPRIMTOCONS 
// Prototype for Fortran procedure PRIMTOCONS ...
//
void FORTRAN_NAME( PRIMTOCONS ,primtocons )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_BOX(box) );

#define FORT_PRIMTOCONS FORTRAN_NAME( inlinePRIMTOCONS, inlinePRIMTOCONS)
#define FORTNT_PRIMTOCONS FORTRAN_NAME( PRIMTOCONS, primtocons)

inline void FORTRAN_NAME(inlinePRIMTOCONS, inlinePRIMTOCONS)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRIMTOCONS");
 FORTRAN_NAME( PRIMTOCONS ,primtocons )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRIMTOCONS 

#ifndef GUARDFLUXESHANCOCK 
#define GUARDFLUXESHANCOCK 
// Prototype for Fortran procedure FLUXESHANCOCK ...
//
void FORTRAN_NAME( FLUXESHANCOCK ,fluxeshancock )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_FLUXESHANCOCK FORTRAN_NAME( inlineFLUXESHANCOCK, inlineFLUXESHANCOCK)
#define FORTNT_FLUXESHANCOCK FORTRAN_NAME( FLUXESHANCOCK, fluxeshancock)

inline void FORTRAN_NAME(inlineFLUXESHANCOCK, inlineFLUXESHANCOCK)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXESHANCOCK");
 FORTRAN_NAME( FLUXESHANCOCK ,fluxeshancock )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXESHANCOCK 

#ifndef GUARDGETDIVBSTEP2 
#define GUARDGETDIVBSTEP2 
// Prototype for Fortran procedure GETDIVBSTEP2 ...
//
void FORTRAN_NAME( GETDIVBSTEP2 ,getdivbstep2 )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA(duB)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(dx) );

#define FORT_GETDIVBSTEP2 FORTRAN_NAME( inlineGETDIVBSTEP2, inlineGETDIVBSTEP2)
#define FORTNT_GETDIVBSTEP2 FORTRAN_NAME( GETDIVBSTEP2, getdivbstep2)

inline void FORTRAN_NAME(inlineGETDIVBSTEP2, inlineGETDIVBSTEP2)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA(duB)
      ,CHFp_BOX(box)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_GETDIVBSTEP2");
 FORTRAN_NAME( GETDIVBSTEP2 ,getdivbstep2 )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA(duB)
      ,CHFt_BOX(box)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDGETDIVBSTEP2 

#ifndef GUARDCORRECT_B 
#define GUARDCORRECT_B 
// Prototype for Fortran procedure CORRECT_B ...
//
void FORTRAN_NAME( CORRECT_B ,correct_b )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_CORRECT_B FORTRAN_NAME( inlineCORRECT_B, inlineCORRECT_B)
#define FORTNT_CORRECT_B FORTRAN_NAME( CORRECT_B, correct_b)

inline void FORTRAN_NAME(inlineCORRECT_B, inlineCORRECT_B)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CORRECT_B");
 FORTRAN_NAME( CORRECT_B ,correct_b )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCORRECT_B 

#ifndef GUARDCORRECTBN 
#define GUARDCORRECTBN 
// Prototype for Fortran procedure CORRECTBN ...
//
void FORTRAN_NAME( CORRECTBN ,correctbn )(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_CORRECTBN FORTRAN_NAME( inlineCORRECTBN, inlineCORRECTBN)
#define FORTNT_CORRECTBN FORTRAN_NAME( CORRECTBN, correctbn)

inline void FORTRAN_NAME(inlineCORRECTBN, inlineCORRECTBN)(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA1(phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CORRECTBN");
 FORTRAN_NAME( CORRECTBN ,correctbn )(
      CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA1(phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCORRECTBN 

#ifndef GUARDCOMPUTEBN 
#define GUARDCOMPUTEBN 
// Prototype for Fortran procedure COMPUTEBN ...
//
void FORTRAN_NAME( COMPUTEBN ,computebn )(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_COMPUTEBN FORTRAN_NAME( inlineCOMPUTEBN, inlineCOMPUTEBN)
#define FORTNT_COMPUTEBN FORTRAN_NAME( COMPUTEBN, computebn)

inline void FORTRAN_NAME(inlineCOMPUTEBN, inlineCOMPUTEBN)(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_COMPUTEBN");
 FORTRAN_NAME( COMPUTEBN ,computebn )(
      CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDCOMPUTEBN 

#ifndef GUARDCOMPUTEBNCD 
#define GUARDCOMPUTEBNCD 
// Prototype for Fortran procedure COMPUTEBNCD ...
//
void FORTRAN_NAME( COMPUTEBNCD ,computebncd )(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_COMPUTEBNCD FORTRAN_NAME( inlineCOMPUTEBNCD, inlineCOMPUTEBNCD)
#define FORTNT_COMPUTEBNCD FORTRAN_NAME( COMPUTEBNCD, computebncd)

inline void FORTRAN_NAME(inlineCOMPUTEBNCD, inlineCOMPUTEBNCD)(
      CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_COMPUTEBNCD");
 FORTRAN_NAME( COMPUTEBNCD ,computebncd )(
      CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDCOMPUTEBNCD 

#ifndef GUARDCOMPUTEUN 
#define GUARDCOMPUTEUN 
// Prototype for Fortran procedure COMPUTEUN ...
//
void FORTRAN_NAME( COMPUTEUN ,computeun )(
      CHFp_FRA1(Un)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_COMPUTEUN FORTRAN_NAME( inlineCOMPUTEUN, inlineCOMPUTEUN)
#define FORTNT_COMPUTEUN FORTRAN_NAME( COMPUTEUN, computeun)

inline void FORTRAN_NAME(inlineCOMPUTEUN, inlineCOMPUTEUN)(
      CHFp_FRA1(Un)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_COMPUTEUN");
 FORTRAN_NAME( COMPUTEUN ,computeun )(
      CHFt_FRA1(Un)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDCOMPUTEUN 

#ifndef GUARDCOMPUTEUNCD 
#define GUARDCOMPUTEUNCD 
// Prototype for Fortran procedure COMPUTEUNCD ...
//
void FORTRAN_NAME( COMPUTEUNCD ,computeuncd )(
      CHFp_FRA1(Un)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) );

#define FORT_COMPUTEUNCD FORTRAN_NAME( inlineCOMPUTEUNCD, inlineCOMPUTEUNCD)
#define FORTNT_COMPUTEUNCD FORTRAN_NAME( COMPUTEUNCD, computeuncd)

inline void FORTRAN_NAME(inlineCOMPUTEUNCD, inlineCOMPUTEUNCD)(
      CHFp_FRA1(Un)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_COMPUTEUNCD");
 FORTRAN_NAME( COMPUTEUNCD ,computeuncd )(
      CHFt_FRA1(Un)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDCOMPUTEUNCD 

#ifndef GUARDCOMPUTEDIVB 
#define GUARDCOMPUTEDIVB 
// Prototype for Fortran procedure COMPUTEDIVB ...
//
void FORTRAN_NAME( COMPUTEDIVB ,computedivb )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) );

#define FORT_COMPUTEDIVB FORTRAN_NAME( inlineCOMPUTEDIVB, inlineCOMPUTEDIVB)
#define FORTNT_COMPUTEDIVB FORTRAN_NAME( COMPUTEDIVB, computedivb)

inline void FORTRAN_NAME(inlineCOMPUTEDIVB, inlineCOMPUTEDIVB)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVB");
 FORTRAN_NAME( COMPUTEDIVB ,computedivb )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA1(Bn)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(Box) );
}
#endif  // GUARDCOMPUTEDIVB 

#ifndef GUARDCOMPUTEDIVB_CYL 
#define GUARDCOMPUTEDIVB_CYL 
// Prototype for Fortran procedure COMPUTEDIVB_CYL ...
//
void FORTRAN_NAME( COMPUTEDIVB_CYL ,computedivb_cyl )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) );

#define FORT_COMPUTEDIVB_CYL FORTRAN_NAME( inlineCOMPUTEDIVB_CYL, inlineCOMPUTEDIVB_CYL)
#define FORTNT_COMPUTEDIVB_CYL FORTRAN_NAME( COMPUTEDIVB_CYL, computedivb_cyl)

inline void FORTRAN_NAME(inlineCOMPUTEDIVB_CYL, inlineCOMPUTEDIVB_CYL)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVB_CYL");
 FORTRAN_NAME( COMPUTEDIVB_CYL ,computedivb_cyl )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA1(Bn)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(Box) );
}
#endif  // GUARDCOMPUTEDIVB_CYL 

#ifndef GUARDCOMPUTEDIVB_POLAR 
#define GUARDCOMPUTEDIVB_POLAR 
// Prototype for Fortran procedure COMPUTEDIVB_POLAR ...
//
void FORTRAN_NAME( COMPUTEDIVB_POLAR ,computedivb_polar )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) );

#define FORT_COMPUTEDIVB_POLAR FORTRAN_NAME( inlineCOMPUTEDIVB_POLAR, inlineCOMPUTEDIVB_POLAR)
#define FORTNT_COMPUTEDIVB_POLAR FORTRAN_NAME( COMPUTEDIVB_POLAR, computedivb_polar)

inline void FORTRAN_NAME(inlineCOMPUTEDIVB_POLAR, inlineCOMPUTEDIVB_POLAR)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVB_POLAR");
 FORTRAN_NAME( COMPUTEDIVB_POLAR ,computedivb_polar )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA1(Bn)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(Box) );
}
#endif  // GUARDCOMPUTEDIVB_POLAR 

#ifndef GUARDCOMPUTEDIVB_SPHERICAL 
#define GUARDCOMPUTEDIVB_SPHERICAL 
// Prototype for Fortran procedure COMPUTEDIVB_SPHERICAL ...
//
void FORTRAN_NAME( COMPUTEDIVB_SPHERICAL ,computedivb_spherical )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) );

#define FORT_COMPUTEDIVB_SPHERICAL FORTRAN_NAME( inlineCOMPUTEDIVB_SPHERICAL, inlineCOMPUTEDIVB_SPHERICAL)
#define FORTNT_COMPUTEDIVB_SPHERICAL FORTRAN_NAME( COMPUTEDIVB_SPHERICAL, computedivb_spherical)

inline void FORTRAN_NAME(inlineCOMPUTEDIVB_SPHERICAL, inlineCOMPUTEDIVB_SPHERICAL)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA1(Bn)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVB_SPHERICAL");
 FORTRAN_NAME( COMPUTEDIVB_SPHERICAL ,computedivb_spherical )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA1(Bn)
      ,CHFt_CONST_FRA1(areas)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(Box) );
}
#endif  // GUARDCOMPUTEDIVB_SPHERICAL 


#ifndef GUARDCOMPUTEDIVU 
#define GUARDCOMPUTEDIVU
// Prototype for Fortran procedure COMPUTEDIVU ...
//
void FORTRAN_NAME( COMPUTEDIVU ,computedivu )(
      CHFp_FRA1(divU)
      ,CHFp_CONST_FRA1(Un)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) );

#define FORT_COMPUTEDIVU FORTRAN_NAME( inlineCOMPUTEDIVU, inlineCOMPUTEDIVU)
#define FORTNT_COMPUTEDIVU FORTRAN_NAME( COMPUTEDIVU, computedivu)

inline void FORTRAN_NAME(inlineCOMPUTEDIVU, inlineCOMPUTEDIVU)(
      CHFp_FRA1(divU)
      ,CHFp_CONST_FRA1(Un)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(Box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVU");
 FORTRAN_NAME( COMPUTEDIVU ,computedivu )(
      CHFt_FRA1(divU)
      ,CHFt_CONST_FRA1(Un)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(Box) );
}
#endif  // GUARDCOMPUTEDIVU


#ifndef GUARDCOMPUTEDIVU_SPHERICAL 
#define GUARDCOMPUTEDIVU_SPHERICAL 
// Prototype for Fortran procedure COMPUTEDIVU_SPHERICAL ...
//
void FORTRAN_NAME( COMPUTEDIVU_SPHERICAL ,computedivu_spherical )(
      CHFp_FRA1(divU)
      ,CHFp_CONST_FRA1(Un)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_COMPUTEDIVU_SPHERICAL FORTRAN_NAME( inlineCOMPUTEDIVU_SPHERICAL, inlineCOMPUTEDIVU_SPHERICAL)
#define FORTNT_COMPUTEDIVU_SPHERICAL FORTRAN_NAME( COMPUTEDIVU_SPHERICAL, computedivu_spherical)

inline void FORTRAN_NAME(inlineCOMPUTEDIVU_SPHERICAL, inlineCOMPUTEDIVU_SPHERICAL)(
      CHFp_FRA1(divU)
      ,CHFp_CONST_FRA1(Un)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPUTEDIVU_SPHERICAL");
 FORTRAN_NAME( COMPUTEDIVU_SPHERICAL ,computedivu_spherical )(
      CHFt_FRA1(divU)
      ,CHFt_CONST_FRA1(Un)
      ,CHFt_CONST_FRA1(areas)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPUTEDIVU_SPHERICAL 

#ifndef GUARDAXISYMMETRICDBZBYDZ 
#define GUARDAXISYMMETRICDBZBYDZ 
// Prototype for Fortran procedure AXISYMMETRICDBZBYDZ ...
//
void FORTRAN_NAME( AXISYMMETRICDBZBYDZ ,axisymmetricdbzbydz )(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_AXISYMMETRICDBZBYDZ FORTRAN_NAME( inlineAXISYMMETRICDBZBYDZ, inlineAXISYMMETRICDBZBYDZ)
#define FORTNT_AXISYMMETRICDBZBYDZ FORTRAN_NAME( AXISYMMETRICDBZBYDZ, axisymmetricdbzbydz)

inline void FORTRAN_NAME(inlineAXISYMMETRICDBZBYDZ, inlineAXISYMMETRICDBZBYDZ)(
      CHFp_FRA1(divB)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_AXISYMMETRICDBZBYDZ");
 FORTRAN_NAME( AXISYMMETRICDBZBYDZ ,axisymmetricdbzbydz )(
      CHFt_FRA1(divB)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDAXISYMMETRICDBZBYDZ 

#ifndef GUARDSOURCE8WAVES 
#define GUARDSOURCE8WAVES 
// Prototype for Fortran procedure SOURCE8WAVES ...
//
void FORTRAN_NAME( SOURCE8WAVES ,source8waves )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_BOX(box) );

#define FORT_SOURCE8WAVES FORTRAN_NAME( inlineSOURCE8WAVES, inlineSOURCE8WAVES)
#define FORTNT_SOURCE8WAVES FORTRAN_NAME( SOURCE8WAVES, source8waves)

inline void FORTRAN_NAME(inlineSOURCE8WAVES, inlineSOURCE8WAVES)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dtbydx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCE8WAVES");
 FORTRAN_NAME( SOURCE8WAVES ,source8waves )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(divBFra)
      ,CHFt_CONST_REAL(dtbydx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCE8WAVES 

#ifndef GUARDSOURCE8WAVES_POLAR 
#define GUARDSOURCE8WAVES_POLAR 
// Prototype for Fortran procedure SOURCE8WAVES_POLAR ...
//
void FORTRAN_NAME( SOURCE8WAVES_POLAR ,source8waves_polar )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCE8WAVES_POLAR FORTRAN_NAME( inlineSOURCE8WAVES_POLAR, inlineSOURCE8WAVES_POLAR)
#define FORTNT_SOURCE8WAVES_POLAR FORTRAN_NAME( SOURCE8WAVES_POLAR, source8waves_polar)

inline void FORTRAN_NAME(inlineSOURCE8WAVES_POLAR, inlineSOURCE8WAVES_POLAR)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCE8WAVES_POLAR");
 FORTRAN_NAME( SOURCE8WAVES_POLAR ,source8waves_polar )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(divBFra)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCE8WAVES_POLAR 

#ifndef GUARDSOURCE8WAVES_SPHERICAL 
#define GUARDSOURCE8WAVES_SPHERICAL 
// Prototype for Fortran procedure SOURCE8WAVES_SPHERICAL ...
//
void FORTRAN_NAME( SOURCE8WAVES_SPHERICAL ,source8waves_spherical )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCE8WAVES_SPHERICAL FORTRAN_NAME( inlineSOURCE8WAVES_SPHERICAL, inlineSOURCE8WAVES_SPHERICAL)
#define FORTNT_SOURCE8WAVES_SPHERICAL FORTRAN_NAME( SOURCE8WAVES_SPHERICAL, source8waves_spherical)

inline void FORTRAN_NAME(inlineSOURCE8WAVES_SPHERICAL, inlineSOURCE8WAVES_SPHERICAL)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_FRA1(invvol)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCE8WAVES_SPHERICAL");
 FORTRAN_NAME( SOURCE8WAVES_SPHERICAL ,source8waves_spherical )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(divBFra)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_FRA1(invvol)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCE8WAVES_SPHERICAL 

#ifndef GUARDSOURCEAXISYMMETRIC 
#define GUARDSOURCEAXISYMMETRIC 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC ,sourceaxisymmetric )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC FORTRAN_NAME( inlineSOURCEAXISYMMETRIC, inlineSOURCEAXISYMMETRIC)
#define FORTNT_SOURCEAXISYMMETRIC FORTRAN_NAME( SOURCEAXISYMMETRIC, sourceaxisymmetric)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC, inlineSOURCEAXISYMMETRIC)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC");
 FORTRAN_NAME( SOURCEAXISYMMETRIC ,sourceaxisymmetric )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC 

#ifndef GUARDSOURCEAXISYMMETRIC_CYL 
#define GUARDSOURCEAXISYMMETRIC_CYL 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_CYL ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL ,sourceaxisymmetric_cyl )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_CYL FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_CYL, inlineSOURCEAXISYMMETRIC_CYL)
#define FORTNT_SOURCEAXISYMMETRIC_CYL FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL, sourceaxisymmetric_cyl)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_CYL, inlineSOURCEAXISYMMETRIC_CYL)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_CYL");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL ,sourceaxisymmetric_cyl )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_CYL 

#ifndef GUARDSOURCEAXISYMMETRIC_POLAR 
#define GUARDSOURCEAXISYMMETRIC_POLAR 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_POLAR ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR ,sourceaxisymmetric_polar )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_POLAR FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_POLAR, inlineSOURCEAXISYMMETRIC_POLAR)
#define FORTNT_SOURCEAXISYMMETRIC_POLAR FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR, sourceaxisymmetric_polar)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_POLAR, inlineSOURCEAXISYMMETRIC_POLAR)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_POLAR");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR ,sourceaxisymmetric_polar )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_POLAR 

#ifndef GUARDPOSTPROCESSING 
#define GUARDPOSTPROCESSING 
// Prototype for Fortran procedure POSTPROCESSING ...
//
void FORTRAN_NAME( POSTPROCESSING ,postprocessing )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING FORTRAN_NAME( inlinePOSTPROCESSING, inlinePOSTPROCESSING)
#define FORTNT_POSTPROCESSING FORTRAN_NAME( POSTPROCESSING, postprocessing)

inline void FORTRAN_NAME(inlinePOSTPROCESSING, inlinePOSTPROCESSING)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING");
 FORTRAN_NAME( POSTPROCESSING ,postprocessing )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Uold)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING 

#ifndef GUARDPOSTPROCESSINGDT 
#define GUARDPOSTPROCESSINGDT 
// Prototype for Fortran procedure POSTPROCESSINGDT ...
//
void FORTRAN_NAME( POSTPROCESSINGDT ,postprocessingdt )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSINGDT FORTRAN_NAME( inlinePOSTPROCESSINGDT, inlinePOSTPROCESSINGDT)
#define FORTNT_POSTPROCESSINGDT FORTRAN_NAME( POSTPROCESSINGDT, postprocessingdt)

inline void FORTRAN_NAME(inlinePOSTPROCESSINGDT, inlinePOSTPROCESSINGDT)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSINGDT");
 FORTRAN_NAME( POSTPROCESSINGDT ,postprocessingdt )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Uold)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSINGDT 

#ifndef GUARDCHARANALYSISF 
#define GUARDCHARANALYSISF 
// Prototype for Fortran procedure CHARANALYSISF ...
//
void FORTRAN_NAME( CHARANALYSISF ,charanalysisf )(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_CHARANALYSISF FORTRAN_NAME( inlineCHARANALYSISF, inlineCHARANALYSISF)
#define FORTNT_CHARANALYSISF FORTRAN_NAME( CHARANALYSISF, charanalysisf)

inline void FORTRAN_NAME(inlineCHARANALYSISF, inlineCHARANALYSISF)(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARANALYSISF");
 FORTRAN_NAME( CHARANALYSISF ,charanalysisf )(
      CHFt_FRA(dWLeft)
      ,CHFt_FRA(dWRight)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARANALYSISF 

#ifndef GUARDCHARSYNTHESISF 
#define GUARDCHARSYNTHESISF 
// Prototype for Fortran procedure CHARSYNTHESISF ...
//
void FORTRAN_NAME( CHARSYNTHESISF ,charsynthesisf )(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_CHARSYNTHESISF FORTRAN_NAME( inlineCHARSYNTHESISF, inlineCHARSYNTHESISF)
#define FORTNT_CHARSYNTHESISF FORTRAN_NAME( CHARSYNTHESISF, charsynthesisf)

inline void FORTRAN_NAME(inlineCHARSYNTHESISF, inlineCHARSYNTHESISF)(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARSYNTHESISF");
 FORTRAN_NAME( CHARSYNTHESISF ,charsynthesisf )(
      CHFt_FRA(dWLeft)
      ,CHFt_FRA(dWRight)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARSYNTHESISF 

#ifndef GUARDREADFORT55 
#define GUARDREADFORT55 
// Prototype for Fortran procedure READFORT55 ...
//
void FORTRAN_NAME( READFORT55 ,readfort55 )(
      CHFp_FRA(Wall) );

#define FORT_READFORT55 FORTRAN_NAME( inlineREADFORT55, inlineREADFORT55)
#define FORTNT_READFORT55 FORTRAN_NAME( READFORT55, readfort55)

inline void FORTRAN_NAME(inlineREADFORT55, inlineREADFORT55)(
      CHFp_FRA(Wall) )
{
 CH_TIMELEAF("FORT_READFORT55");
 FORTRAN_NAME( READFORT55 ,readfort55 )(
      CHFt_FRA(Wall) );
}
#endif  // GUARDREADFORT55 

}

#endif
