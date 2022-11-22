#ifndef _HELIOSPHERICPLASMABCF_F_H_
#define _HELIOSPHERICPLASMABCF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDHELIOINITPLASMASPHERICAL 
#define GUARDHELIOINITPLASMASPHERICAL 
// Prototype for Fortran procedure HELIOINITPLASMASPHERICAL ...
//
void FORTRAN_NAME( HELIOINITPLASMASPHERICAL ,helioinitplasmaspherical )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(subproblem)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOINITPLASMASPHERICAL FORTRAN_NAME( inlineHELIOINITPLASMASPHERICAL, inlineHELIOINITPLASMASPHERICAL)
#define FORTNT_HELIOINITPLASMASPHERICAL FORTRAN_NAME( HELIOINITPLASMASPHERICAL, helioinitplasmaspherical)

inline void FORTRAN_NAME(inlineHELIOINITPLASMASPHERICAL, inlineHELIOINITPLASMASPHERICAL)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(subproblem)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINITPLASMASPHERICAL");
 FORTRAN_NAME( HELIOINITPLASMASPHERICAL ,helioinitplasmaspherical )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(subproblem)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINITPLASMASPHERICAL 

#ifndef GUARDHELIOGSPLASMASPHERICAL 
#define GUARDHELIOGSPLASMASPHERICAL 
// Prototype for Fortran procedure HELIOGSPLASMASPHERICAL ...
//
void FORTRAN_NAME( HELIOGSPLASMASPHERICAL ,heliogsplasmaspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(subproblem)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSPLASMASPHERICAL FORTRAN_NAME( inlineHELIOGSPLASMASPHERICAL, inlineHELIOGSPLASMASPHERICAL)
#define FORTNT_HELIOGSPLASMASPHERICAL FORTRAN_NAME( HELIOGSPLASMASPHERICAL, heliogsplasmaspherical)

inline void FORTRAN_NAME(inlineHELIOGSPLASMASPHERICAL, inlineHELIOGSPLASMASPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(subproblem)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSPLASMASPHERICAL");
 FORTRAN_NAME( HELIOGSPLASMASPHERICAL ,heliogsplasmaspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(subproblem)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSPLASMASPHERICAL 

#ifndef GUARDHELIOINITPLASMASPHDEFAULT 
#define GUARDHELIOINITPLASMASPHDEFAULT 
// Prototype for Fortran procedure helioinitplasmasphdefault ...
//
void FORTRAN_NAME( HELIOINITPLASMASPHDEFAULT ,helioinitplasmasphdefault )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOINITPLASMASPHDEFAULT FORTRAN_NAME( inlineHELIOINITPLASMASPHDEFAULT, inlineHELIOINITPLASMASPHDEFAULT)
#define FORTNT_HELIOINITPLASMASPHDEFAULT FORTRAN_NAME( HELIOINITPLASMASPHDEFAULT, helioinitplasmasphdefault)

inline void FORTRAN_NAME(inlineHELIOINITPLASMASPHDEFAULT, inlineHELIOINITPLASMASPHDEFAULT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINITPLASMASPHDEFAULT");
 FORTRAN_NAME( HELIOINITPLASMASPHDEFAULT ,helioinitplasmasphdefault )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINITPLASMASPHDEFAULT 

#ifndef GUARDHELIOGSPLASMASPHDEFAULT 
#define GUARDHELIOGSPLASMASPHDEFAULT 
// Prototype for Fortran procedure heliogsplasmasphdefault ...
//
void FORTRAN_NAME( HELIOGSPLASMASPHDEFAULT ,heliogsplasmasphdefault )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSPLASMASPHDEFAULT FORTRAN_NAME( inlineHELIOGSPLASMASPHDEFAULT, inlineHELIOGSPLASMASPHDEFAULT)
#define FORTNT_HELIOGSPLASMASPHDEFAULT FORTRAN_NAME( HELIOGSPLASMASPHDEFAULT, heliogsplasmasphdefault)

inline void FORTRAN_NAME(inlineHELIOGSPLASMASPHDEFAULT, inlineHELIOGSPLASMASPHDEFAULT)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSPLASMASPHDEFAULT");
 FORTRAN_NAME( HELIOGSPLASMASPHDEFAULT ,heliogsplasmasphdefault )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSPLASMASPHDEFAULT 

#ifndef GUARDHELIOGSFIRST 
#define GUARDHELIOGSFIRST 
// Prototype for Fortran procedure heliogsfirst ...
//
void FORTRAN_NAME( HELIOGSFIRST ,heliogsfirst )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSFIRST FORTRAN_NAME( inlineHELIOGSFIRST, inlineHELIOGSFIRST)
#define FORTNT_HELIOGSFIRST FORTRAN_NAME( HELIOGSFIRST, heliogsfirst)

inline void FORTRAN_NAME(inlineHELIOGSFIRST, inlineHELIOGSFIRST)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSFIRST");
 FORTRAN_NAME( HELIOGSFIRST ,heliogsfirst )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSFIRST 

#ifndef GUARDHELIOGSSUESS_MIN_FIRST 
#define GUARDHELIOGSSUESS_MIN_FIRST 
// Prototype for Fortran procedure heliogsSuess_min_first ...
//
void FORTRAN_NAME( HELIOGSSUESS_MIN_FIRST ,heliogssuess_min_first )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSSUESS_MIN_FIRST FORTRAN_NAME( inlineHELIOGSSUESS_MIN_FIRST, inlineHELIOGSSUESS_MIN_FIRST)
#define FORTNT_HELIOGSSUESS_MIN_FIRST FORTRAN_NAME( HELIOGSSUESS_MIN_FIRST, heliogssuess_min_first)

inline void FORTRAN_NAME(inlineHELIOGSSUESS_MIN_FIRST, inlineHELIOGSSUESS_MIN_FIRST)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSSUESS_MIN_FIRST");
 FORTRAN_NAME( HELIOGSSUESS_MIN_FIRST ,heliogssuess_min_first )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSSUESS_MIN_FIRST 

#ifndef GUARDHELIOGSSUESS_MIN_THIRD 
#define GUARDHELIOGSSUESS_MIN_THIRD 
// Prototype for Fortran procedure heliogsSuess_min_third ...
//
void FORTRAN_NAME( HELIOGSSUESS_MIN_THIRD ,heliogssuess_min_third )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSSUESS_MIN_THIRD FORTRAN_NAME( inlineHELIOGSSUESS_MIN_THIRD, inlineHELIOGSSUESS_MIN_THIRD)
#define FORTNT_HELIOGSSUESS_MIN_THIRD FORTRAN_NAME( HELIOGSSUESS_MIN_THIRD, heliogssuess_min_third)

inline void FORTRAN_NAME(inlineHELIOGSSUESS_MIN_THIRD, inlineHELIOGSSUESS_MIN_THIRD)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSSUESS_MIN_THIRD");
 FORTRAN_NAME( HELIOGSSUESS_MIN_THIRD ,heliogssuess_min_third )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSSUESS_MIN_THIRD 

#ifndef GUARDHELIOGSTHIRD 
#define GUARDHELIOGSTHIRD 
// Prototype for Fortran procedure heliogsthird ...
//
void FORTRAN_NAME( HELIOGSTHIRD ,heliogsthird )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSTHIRD FORTRAN_NAME( inlineHELIOGSTHIRD, inlineHELIOGSTHIRD)
#define FORTNT_HELIOGSTHIRD FORTRAN_NAME( HELIOGSTHIRD, heliogsthird)

inline void FORTRAN_NAME(inlineHELIOGSTHIRD, inlineHELIOGSTHIRD)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSTHIRD");
 FORTRAN_NAME( HELIOGSTHIRD ,heliogsthird )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSTHIRD 

#ifndef GUARDHELIOGSCIR 
#define GUARDHELIOGSCIR 
// Prototype for Fortran procedure heliogscir ...
//
void FORTRAN_NAME( HELIOGSCIR ,heliogscir )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSCIR FORTRAN_NAME( inlineHELIOGSCIR, inlineHELIOGSCIR)
#define FORTNT_HELIOGSCIR FORTRAN_NAME( HELIOGSCIR, heliogscir)

inline void FORTRAN_NAME(inlineHELIOGSCIR, inlineHELIOGSCIR)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSCIR");
 FORTRAN_NAME( HELIOGSCIR ,heliogscir )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSCIR 

#ifndef GUARDHELIOGSPIZZO 
#define GUARDHELIOGSPIZZO 
// Prototype for Fortran procedure heliogspizzo ...
//
void FORTRAN_NAME( HELIOGSPIZZO ,heliogspizzo )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSPIZZO FORTRAN_NAME( inlineHELIOGSPIZZO, inlineHELIOGSPIZZO)
#define FORTNT_HELIOGSPIZZO FORTRAN_NAME( HELIOGSPIZZO, heliogspizzo)

inline void FORTRAN_NAME(inlineHELIOGSPIZZO, inlineHELIOGSPIZZO)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSPIZZO");
 FORTRAN_NAME( HELIOGSPIZZO ,heliogspizzo )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSPIZZO 

#ifndef GUARDHELIOREINIT_DEFAULT 
#define GUARDHELIOREINIT_DEFAULT 
// Prototype for Fortran procedure HELIOREINIT_DEFAULT ...
//
void FORTRAN_NAME( HELIOREINIT_DEFAULT ,helioreinit_default )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINIT_DEFAULT FORTRAN_NAME( inlineHELIOREINIT_DEFAULT, inlineHELIOREINIT_DEFAULT)
#define FORTNT_HELIOREINIT_DEFAULT FORTRAN_NAME( HELIOREINIT_DEFAULT, helioreinit_default)

inline void FORTRAN_NAME(inlineHELIOREINIT_DEFAULT, inlineHELIOREINIT_DEFAULT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINIT_DEFAULT");
 FORTRAN_NAME( HELIOREINIT_DEFAULT ,helioreinit_default )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_REAL(Rinit)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINIT_DEFAULT 

#ifndef GUARDHELIOREINIT_CYCLE 
#define GUARDHELIOREINIT_CYCLE 
// Prototype for Fortran procedure HELIOREINIT_CYCLE ...
//
void FORTRAN_NAME( HELIOREINIT_CYCLE ,helioreinit_cycle )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINIT_CYCLE FORTRAN_NAME( inlineHELIOREINIT_CYCLE, inlineHELIOREINIT_CYCLE)
#define FORTNT_HELIOREINIT_CYCLE FORTRAN_NAME( HELIOREINIT_CYCLE, helioreinit_cycle)

inline void FORTRAN_NAME(inlineHELIOREINIT_CYCLE, inlineHELIOREINIT_CYCLE)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINIT_CYCLE");
 FORTRAN_NAME( HELIOREINIT_CYCLE ,helioreinit_cycle )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_REAL(Rinit)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINIT_CYCLE 

#ifndef GUARDHELIOREINIT_CYCLE_SPHERICAL 
#define GUARDHELIOREINIT_CYCLE_SPHERICAL 
// Prototype for Fortran procedure HELIOREINIT_CYCLE_SPHERICAL ...
//
void FORTRAN_NAME( HELIOREINIT_CYCLE_SPHERICAL ,helioreinit_cycle_spherical )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOREINIT_CYCLE_SPHERICAL FORTRAN_NAME( inlineHELIOREINIT_CYCLE_SPHERICAL, inlineHELIOREINIT_CYCLE_SPHERICAL)
#define FORTNT_HELIOREINIT_CYCLE_SPHERICAL FORTRAN_NAME( HELIOREINIT_CYCLE_SPHERICAL, helioreinit_cycle_spherical)

inline void FORTRAN_NAME(inlineHELIOREINIT_CYCLE_SPHERICAL, inlineHELIOREINIT_CYCLE_SPHERICAL)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(t)
      ,CHFp_CONST_REAL(Rinit)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOREINIT_CYCLE_SPHERICAL");
 FORTRAN_NAME( HELIOREINIT_CYCLE_SPHERICAL ,helioreinit_cycle_spherical )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(t)
      ,CHFt_CONST_REAL(Rinit)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOREINIT_CYCLE_SPHERICAL 

}

#endif
