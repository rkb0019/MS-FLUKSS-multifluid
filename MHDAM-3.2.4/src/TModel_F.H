#ifndef _TMODEL_F_H_
#define _TMODEL_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCONST_TM 
#define GUARDSETCONST_TM 
// Prototype for Fortran procedure SETCONST_TM ...
//
void FORTRAN_NAME( SETCONST_TM ,setconst_tm )(
      CHFp_CONST_REAL(auAlpha)
      ,CHFp_CONST_REAL(auBeta)
      ,CHFp_CONST_REAL(auSmallTM)
      ,CHFp_CONST_REAL(auSigmaD)
      ,CHFp_CONST_REAL(auFD)
      ,CHFp_CONST_REAL(auNH)
      ,CHFp_CONST_REAL(auTion)
      ,CHFp_CONST_REAL(auLcav)
      ,CHFp_CONST_REAL(auUr1AU)
      ,CHFp_CONST_REAL(auVa1AU)
      ,CHFp_CONST_REAL(auNSW1AU) );

#define FORT_SETCONST_TM FORTRAN_NAME( inlineSETCONST_TM, inlineSETCONST_TM)
#define FORTNT_SETCONST_TM FORTRAN_NAME( SETCONST_TM, setconst_tm)

inline void FORTRAN_NAME(inlineSETCONST_TM, inlineSETCONST_TM)(
      CHFp_CONST_REAL(auAlpha)
      ,CHFp_CONST_REAL(auBeta)
      ,CHFp_CONST_REAL(auSmallTM)
      ,CHFp_CONST_REAL(auSigmaD)
      ,CHFp_CONST_REAL(auFD)
      ,CHFp_CONST_REAL(auNH)
      ,CHFp_CONST_REAL(auTion)
      ,CHFp_CONST_REAL(auLcav)
      ,CHFp_CONST_REAL(auUr1AU)
      ,CHFp_CONST_REAL(auVa1AU)
      ,CHFp_CONST_REAL(auNSW1AU) )
{
 CH_TIMELEAF("FORT_SETCONST_TM");
 FORTRAN_NAME( SETCONST_TM ,setconst_tm )(
      CHFt_CONST_REAL(auAlpha)
      ,CHFt_CONST_REAL(auBeta)
      ,CHFt_CONST_REAL(auSmallTM)
      ,CHFt_CONST_REAL(auSigmaD)
      ,CHFt_CONST_REAL(auFD)
      ,CHFt_CONST_REAL(auNH)
      ,CHFt_CONST_REAL(auTion)
      ,CHFt_CONST_REAL(auLcav)
      ,CHFt_CONST_REAL(auUr1AU)
      ,CHFt_CONST_REAL(auVa1AU)
      ,CHFt_CONST_REAL(auNSW1AU) );
}
#endif  // GUARDSETCONST_TM 

#ifndef GUARDSETSCALES_FOR_EPI 
#define GUARDSETSCALES_FOR_EPI 
// Prototype for Fortran procedure SETSCALES_FOR_EPI ...
//
void FORTRAN_NAME( SETSCALES_FOR_EPI ,setscales_for_epi )(
      CHFp_CONST_REAL(auRefLen)
      ,CHFp_CONST_REAL(auRefVel) );

#define FORT_SETSCALES_FOR_EPI FORTRAN_NAME( inlineSETSCALES_FOR_EPI, inlineSETSCALES_FOR_EPI)
#define FORTNT_SETSCALES_FOR_EPI FORTRAN_NAME( SETSCALES_FOR_EPI, setscales_for_epi)

inline void FORTRAN_NAME(inlineSETSCALES_FOR_EPI, inlineSETSCALES_FOR_EPI)(
      CHFp_CONST_REAL(auRefLen)
      ,CHFp_CONST_REAL(auRefVel) )
{
 CH_TIMELEAF("FORT_SETSCALES_FOR_EPI");
 FORTRAN_NAME( SETSCALES_FOR_EPI ,setscales_for_epi )(
      CHFt_CONST_REAL(auRefLen)
      ,CHFt_CONST_REAL(auRefVel) );
}
#endif  // GUARDSETSCALES_FOR_EPI 

#ifndef GUARDCONSTOPRIM_TM 
#define GUARDCONSTOPRIM_TM 
// Prototype for Fortran procedure CONSTOPRIM_TM ...
//
void FORTRAN_NAME( CONSTOPRIM_TM ,constoprim_tm )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) );

#define FORT_CONSTOPRIM_TM FORTRAN_NAME( inlineCONSTOPRIM_TM, inlineCONSTOPRIM_TM)
#define FORTNT_CONSTOPRIM_TM FORTRAN_NAME( CONSTOPRIM_TM, constoprim_tm)

inline void FORTRAN_NAME(inlineCONSTOPRIM_TM, inlineCONSTOPRIM_TM)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONSTOPRIM_TM");
 FORTRAN_NAME( CONSTOPRIM_TM ,constoprim_tm )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONSTOPRIM_TM 

#ifndef GUARDPRIMTOCONS_TM 
#define GUARDPRIMTOCONS_TM 
// Prototype for Fortran procedure PRIMTOCONS_TM ...
//
void FORTRAN_NAME( PRIMTOCONS_TM ,primtocons_tm )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) );

#define FORT_PRIMTOCONS_TM FORTRAN_NAME( inlinePRIMTOCONS_TM, inlinePRIMTOCONS_TM)
#define FORTNT_PRIMTOCONS_TM FORTRAN_NAME( PRIMTOCONS_TM, primtocons_tm)

inline void FORTRAN_NAME(inlinePRIMTOCONS_TM, inlinePRIMTOCONS_TM)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRIMTOCONS_TM");
 FORTRAN_NAME( PRIMTOCONS_TM ,primtocons_tm )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRIMTOCONS_TM 

#ifndef GUARDPOSTPROCESSING_TM 
#define GUARDPOSTPROCESSING_TM 
// Prototype for Fortran procedure POSTPROCESSING_TM ...
//
void FORTRAN_NAME( POSTPROCESSING_TM ,postprocessing_tm )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING_TM FORTRAN_NAME( inlinePOSTPROCESSING_TM, inlinePOSTPROCESSING_TM)
#define FORTNT_POSTPROCESSING_TM FORTRAN_NAME( POSTPROCESSING_TM, postprocessing_tm)

inline void FORTRAN_NAME(inlinePOSTPROCESSING_TM, inlinePOSTPROCESSING_TM)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING_TM");
 FORTRAN_NAME( POSTPROCESSING_TM ,postprocessing_tm )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING_TM 

#ifndef GUARDFLUXE_TM 
#define GUARDFLUXE_TM 
// Prototype for Fortran procedure FLUXE_TM ...
//
void FORTRAN_NAME( FLUXE_TM ,fluxe_tm )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_FLUXE_TM FORTRAN_NAME( inlineFLUXE_TM, inlineFLUXE_TM)
#define FORTNT_FLUXE_TM FORTRAN_NAME( FLUXE_TM, fluxe_tm)

inline void FORTRAN_NAME(inlineFLUXE_TM, inlineFLUXE_TM)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXE_TM");
 FORTRAN_NAME( FLUXE_TM ,fluxe_tm )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXE_TM 

#ifndef GUARDUPWINDSCALARFLUXES 
#define GUARDUPWINDSCALARFLUXES 
// Prototype for Fortran procedure UPWINDSCALARFLUXES ...
//
void FORTRAN_NAME( UPWINDSCALARFLUXES ,upwindscalarfluxes )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) );

#define FORT_UPWINDSCALARFLUXES FORTRAN_NAME( inlineUPWINDSCALARFLUXES, inlineUPWINDSCALARFLUXES)
#define FORTNT_UPWINDSCALARFLUXES FORTRAN_NAME( UPWINDSCALARFLUXES, upwindscalarfluxes)

inline void FORTRAN_NAME(inlineUPWINDSCALARFLUXES, inlineUPWINDSCALARFLUXES)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPWINDSCALARFLUXES");
 FORTRAN_NAME( UPWINDSCALARFLUXES ,upwindscalarfluxes )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPWINDSCALARFLUXES 

#ifndef GUARDSOURCEBREECHETAL2008 
#define GUARDSOURCEBREECHETAL2008 
// Prototype for Fortran procedure SOURCEBREECHETAL2008 ...
//
void FORTRAN_NAME( SOURCEBREECHETAL2008 ,sourcebreechetal2008 )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCEBREECHETAL2008 FORTRAN_NAME( inlineSOURCEBREECHETAL2008, inlineSOURCEBREECHETAL2008)
#define FORTNT_SOURCEBREECHETAL2008 FORTRAN_NAME( SOURCEBREECHETAL2008, sourcebreechetal2008)

inline void FORTRAN_NAME(inlineSOURCEBREECHETAL2008, inlineSOURCEBREECHETAL2008)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEBREECHETAL2008");
 FORTRAN_NAME( SOURCEBREECHETAL2008 ,sourcebreechetal2008 )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEBREECHETAL2008 

#ifndef GUARDSOURCEEPIICENBERG 
#define GUARDSOURCEEPIICENBERG 
// Prototype for Fortran procedure SOURCEEPIICENBERG ...
//
void FORTRAN_NAME( SOURCEEPIICENBERG ,sourceepiicenberg )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCEEPIICENBERG FORTRAN_NAME( inlineSOURCEEPIICENBERG, inlineSOURCEEPIICENBERG)
#define FORTNT_SOURCEEPIICENBERG FORTRAN_NAME( SOURCEEPIICENBERG, sourceepiicenberg)

inline void FORTRAN_NAME(inlineSOURCEEPIICENBERG, inlineSOURCEEPIICENBERG)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEEPIICENBERG");
 FORTRAN_NAME( SOURCEEPIICENBERG ,sourceepiicenberg )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEEPIICENBERG 

#ifndef GUARDSOURCEEPIICENBERG_PI 
#define GUARDSOURCEEPIICENBERG_PI 
// Prototype for Fortran procedure SOURCEEPIICENBERG_PI ...
//
void FORTRAN_NAME( SOURCEEPIICENBERG_PI ,sourceepiicenberg_pi )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_SOURCEEPIICENBERG_PI FORTRAN_NAME( inlineSOURCEEPIICENBERG_PI, inlineSOURCEEPIICENBERG_PI)
#define FORTNT_SOURCEEPIICENBERG_PI FORTRAN_NAME( SOURCEEPIICENBERG_PI, sourceepiicenberg_pi)

inline void FORTRAN_NAME(inlineSOURCEEPIICENBERG_PI, inlineSOURCEEPIICENBERG_PI)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoZ2)
      ,CHFp_CONST_INT(iZ2)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEEPIICENBERG_PI");
 FORTRAN_NAME( SOURCEEPIICENBERG_PI ,sourceepiicenberg_pi )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoH)
      ,CHFt_CONST_INT(iRhoZ2)
      ,CHFt_CONST_INT(iZ2)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEEPIICENBERG_PI 

}

#endif
