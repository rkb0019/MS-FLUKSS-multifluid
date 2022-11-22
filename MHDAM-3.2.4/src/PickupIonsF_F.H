#ifndef _PICKUPIONSF_F_H_
#define _PICKUPIONSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCONST_PI 
#define GUARDSETCONST_PI 
// Prototype for Fortran procedure SETCONST_PI ...
//
void FORTRAN_NAME( SETCONST_PI ,setconst_pi )(
      CHFp_CONST_REAL(auGamma)
      ,CHFp_CONST_REAL(auSmallPI) );

#define FORT_SETCONST_PI FORTRAN_NAME( inlineSETCONST_PI, inlineSETCONST_PI)
#define FORTNT_SETCONST_PI FORTRAN_NAME( SETCONST_PI, setconst_pi)

inline void FORTRAN_NAME(inlineSETCONST_PI, inlineSETCONST_PI)(
      CHFp_CONST_REAL(auGamma)
      ,CHFp_CONST_REAL(auSmallPI) )
{
 CH_TIMELEAF("FORT_SETCONST_PI");
 FORTRAN_NAME( SETCONST_PI ,setconst_pi )(
      CHFt_CONST_REAL(auGamma)
      ,CHFt_CONST_REAL(auSmallPI) );
}
#endif  // GUARDSETCONST_PI 

#ifndef GUARDPOSTPROCESSING_PI 
#define GUARDPOSTPROCESSING_PI 
// Prototype for Fortran procedure POSTPROCESSING_PI ...
//
void FORTRAN_NAME( POSTPROCESSING_PI ,postprocessing_pi )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING_PI FORTRAN_NAME( inlinePOSTPROCESSING_PI, inlinePOSTPROCESSING_PI)
#define FORTNT_POSTPROCESSING_PI FORTRAN_NAME( POSTPROCESSING_PI, postprocessing_pi)

inline void FORTRAN_NAME(inlinePOSTPROCESSING_PI, inlinePOSTPROCESSING_PI)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING_PI");
 FORTRAN_NAME( POSTPROCESSING_PI ,postprocessing_pi )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING_PI 

#ifndef GUARDFLUXE_PI 
#define GUARDFLUXE_PI 
// Prototype for Fortran procedure FLUXE_PI ...
//
void FORTRAN_NAME( FLUXE_PI ,fluxe_pi )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_FLUXE_PI FORTRAN_NAME( inlineFLUXE_PI, inlineFLUXE_PI)
#define FORTNT_FLUXE_PI FORTRAN_NAME( FLUXE_PI, fluxe_pi)

inline void FORTRAN_NAME(inlineFLUXE_PI, inlineFLUXE_PI)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXE_PI");
 FORTRAN_NAME( FLUXE_PI ,fluxe_pi )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXE_PI 

#ifndef GUARDFLUXE_PI_REG 
#define GUARDFLUXE_PI_REG 
// Prototype for Fortran procedure FLUXE_PI_REG ...
//
void FORTRAN_NAME( FLUXE_PI_REG ,fluxe_pi_reg )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_FLUXE_PI_REG FORTRAN_NAME( inlineFLUXE_PI_REG, inlineFLUXE_PI_REG)
#define FORTNT_FLUXE_PI_REG FORTRAN_NAME( FLUXE_PI_REG, fluxe_pi_reg)

inline void FORTRAN_NAME(inlineFLUXE_PI_REG, inlineFLUXE_PI_REG)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXE_PI_REG");
 FORTRAN_NAME( FLUXE_PI_REG ,fluxe_pi_reg )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXE_PI_REG 

#ifndef GUARDUPWINDSCALARFLUXES_PI 
#define GUARDUPWINDSCALARFLUXES_PI 
// Prototype for Fortran procedure UPWINDSCALARFLUXES_PI ...
//
void FORTRAN_NAME( UPWINDSCALARFLUXES_PI ,upwindscalarfluxes_pi )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_UPWINDSCALARFLUXES_PI FORTRAN_NAME( inlineUPWINDSCALARFLUXES_PI, inlineUPWINDSCALARFLUXES_PI)
#define FORTNT_UPWINDSCALARFLUXES_PI FORTRAN_NAME( UPWINDSCALARFLUXES_PI, upwindscalarfluxes_pi)

inline void FORTRAN_NAME(inlineUPWINDSCALARFLUXES_PI, inlineUPWINDSCALARFLUXES_PI)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPWINDSCALARFLUXES_PI");
 FORTRAN_NAME( UPWINDSCALARFLUXES_PI ,upwindscalarfluxes_pi )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPWINDSCALARFLUXES_PI 

#ifndef GUARDUPWINDSCALARFLUXES_PI_REG 
#define GUARDUPWINDSCALARFLUXES_PI_REG 
// Prototype for Fortran procedure UPWINDSCALARFLUXES_PI_REG ...
//
void FORTRAN_NAME( UPWINDSCALARFLUXES_PI_REG ,upwindscalarfluxes_pi_reg )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_UPWINDSCALARFLUXES_PI_REG FORTRAN_NAME( inlineUPWINDSCALARFLUXES_PI_REG, inlineUPWINDSCALARFLUXES_PI_REG)
#define FORTNT_UPWINDSCALARFLUXES_PI_REG FORTRAN_NAME( UPWINDSCALARFLUXES_PI_REG, upwindscalarfluxes_pi_reg)

inline void FORTRAN_NAME(inlineUPWINDSCALARFLUXES_PI_REG, inlineUPWINDSCALARFLUXES_PI_REG)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_UPWINDSCALARFLUXES_PI_REG");
 FORTRAN_NAME( UPWINDSCALARFLUXES_PI_REG ,upwindscalarfluxes_pi_reg )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDUPWINDSCALARFLUXES_PI_REG 

#ifndef GUARDSOURCE_PICKUP_IONS 
#define GUARDSOURCE_PICKUP_IONS 
// Prototype for Fortran procedure SOURCE_PICKUP_IONS ...
//
void FORTRAN_NAME( SOURCE_PICKUP_IONS ,source_pickup_ions )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCE_PICKUP_IONS FORTRAN_NAME( inlineSOURCE_PICKUP_IONS, inlineSOURCE_PICKUP_IONS)
#define FORTNT_SOURCE_PICKUP_IONS FORTRAN_NAME( SOURCE_PICKUP_IONS, source_pickup_ions)

inline void FORTRAN_NAME(inlineSOURCE_PICKUP_IONS, inlineSOURCE_PICKUP_IONS)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCE_PICKUP_IONS");
 FORTRAN_NAME( SOURCE_PICKUP_IONS ,source_pickup_ions )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCE_PICKUP_IONS 

#ifndef GUARDSOURCE_PICKUP_IONS_REG 
#define GUARDSOURCE_PICKUP_IONS_REG 
// Prototype for Fortran procedure SOURCE_PICKUP_IONS_REG ...
//
void FORTRAN_NAME( SOURCE_PICKUP_IONS_REG ,source_pickup_ions_reg )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_FRA1(divU)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCE_PICKUP_IONS_REG FORTRAN_NAME( inlineSOURCE_PICKUP_IONS_REG, inlineSOURCE_PICKUP_IONS_REG)
#define FORTNT_SOURCE_PICKUP_IONS_REG FORTRAN_NAME( SOURCE_PICKUP_IONS_REG, source_pickup_ions_reg)

inline void FORTRAN_NAME(inlineSOURCE_PICKUP_IONS_REG, inlineSOURCE_PICKUP_IONS_REG)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_FRA1(divU)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCE_PICKUP_IONS_REG");
 FORTRAN_NAME( SOURCE_PICKUP_IONS_REG ,source_pickup_ions_reg )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_FRA1(divU)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCE_PICKUP_IONS_REG 

#ifndef GUARDMAXWAVESPEED_PI 
#define GUARDMAXWAVESPEED_PI 
// Prototype for Fortran procedure MAXWAVESPEED_PI ...
//
void FORTRAN_NAME( MAXWAVESPEED_PI ,maxwavespeed_pi )(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_MAXWAVESPEED_PI FORTRAN_NAME( inlineMAXWAVESPEED_PI, inlineMAXWAVESPEED_PI)
#define FORTNT_MAXWAVESPEED_PI FORTRAN_NAME( MAXWAVESPEED_PI, maxwavespeed_pi)

inline void FORTRAN_NAME(inlineMAXWAVESPEED_PI, inlineMAXWAVESPEED_PI)(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXWAVESPEED_PI");
 FORTRAN_NAME( MAXWAVESPEED_PI ,maxwavespeed_pi )(
      CHFt_FRA1(cmax)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXWAVESPEED_PI 

#ifndef GUARDMINDT_SPHERICAL_PI 
#define GUARDMINDT_SPHERICAL_PI 
// Prototype for Fortran procedure MINDT_SPHERICAL_PI ...
//
void FORTRAN_NAME( MINDT_SPHERICAL_PI ,mindt_spherical_pi )(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) );

#define FORT_MINDT_SPHERICAL_PI FORTRAN_NAME( inlineMINDT_SPHERICAL_PI, inlineMINDT_SPHERICAL_PI)
#define FORTNT_MINDT_SPHERICAL_PI FORTRAN_NAME( MINDT_SPHERICAL_PI, mindt_spherical_pi)

inline void FORTRAN_NAME(inlineMINDT_SPHERICAL_PI, inlineMINDT_SPHERICAL_PI)(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_INT(iRhoPI)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINDT_SPHERICAL_PI");
 FORTRAN_NAME( MINDT_SPHERICAL_PI ,mindt_spherical_pi )(
      CHFt_FRA1(dt)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_INT(iRhoPI)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINDT_SPHERICAL_PI 

}

#endif
