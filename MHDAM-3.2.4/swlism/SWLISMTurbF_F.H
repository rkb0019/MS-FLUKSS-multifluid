#ifndef _SWLISMTURBF_F_H_
#define _SWLISMTURBF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETSWLISMTURB 
#define GUARDSETSWLISMTURB 
// Prototype for Fortran procedure SETSWLISMTURB ...
//
void FORTRAN_NAME( SETSWLISMTURB ,setswlismturb )(
      CHFp_CONST_REAL(auLISMDeviation)
      ,CHFp_CONST_REAL(auSWDeviation)
      ,CHFp_CONST_REAL(auR0Turb) );

#define FORT_SETSWLISMTURB FORTRAN_NAME( inlineSETSWLISMTURB, inlineSETSWLISMTURB)
#define FORTNT_SETSWLISMTURB FORTRAN_NAME( SETSWLISMTURB, setswlismturb)

inline void FORTRAN_NAME(inlineSETSWLISMTURB, inlineSETSWLISMTURB)(
      CHFp_CONST_REAL(auLISMDeviation)
      ,CHFp_CONST_REAL(auSWDeviation)
      ,CHFp_CONST_REAL(auR0Turb) )
{
 CH_TIMELEAF("FORT_SETSWLISMTURB");
 FORTRAN_NAME( SETSWLISMTURB ,setswlismturb )(
      CHFt_CONST_REAL(auLISMDeviation)
      ,CHFt_CONST_REAL(auSWDeviation)
      ,CHFt_CONST_REAL(auR0Turb) );
}
#endif  // GUARDSETSWLISMTURB 

#ifndef GUARDSWLISMREINIT_TURB 
#define GUARDSWLISMREINIT_TURB 
// Prototype for Fortran procedure SWLISMREINIT_TURB ...
//
void FORTRAN_NAME( SWLISMREINIT_TURB ,swlismreinit_turb )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SWLISMREINIT_TURB FORTRAN_NAME( inlineSWLISMREINIT_TURB, inlineSWLISMREINIT_TURB)
#define FORTNT_SWLISMREINIT_TURB FORTRAN_NAME( SWLISMREINIT_TURB, swlismreinit_turb)

inline void FORTRAN_NAME(inlineSWLISMREINIT_TURB, inlineSWLISMREINIT_TURB)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWLISMREINIT_TURB");
 FORTRAN_NAME( SWLISMREINIT_TURB ,swlismreinit_turb )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWLISMREINIT_TURB 

#ifndef GUARDSW_TURBR0 
#define GUARDSW_TURBR0 
// Prototype for Fortran procedure SW_TURBR0 ...
//
void FORTRAN_NAME( SW_TURBR0 ,sw_turbr0 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SW_TURBR0 FORTRAN_NAME( inlineSW_TURBR0, inlineSW_TURBR0)
#define FORTNT_SW_TURBR0 FORTRAN_NAME( SW_TURBR0, sw_turbr0)

inline void FORTRAN_NAME(inlineSW_TURBR0, inlineSW_TURBR0)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Phi)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SW_TURBR0");
 FORTRAN_NAME( SW_TURBR0 ,sw_turbr0 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Phi)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSW_TURBR0 

#ifndef GUARDREG3_TURB1 
#define GUARDREG3_TURB1 
// Prototype for Fortran procedure REG3_TURB1 ...
//
void FORTRAN_NAME( REG3_TURB1 ,reg3_turb1 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_BOX(box) );

#define FORT_REG3_TURB1 FORTRAN_NAME( inlineREG3_TURB1, inlineREG3_TURB1)
#define FORTNT_REG3_TURB1 FORTRAN_NAME( REG3_TURB1, reg3_turb1)

inline void FORTRAN_NAME(inlineREG3_TURB1, inlineREG3_TURB1)(
      CHFp_FRA(U)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REG3_TURB1");
 FORTRAN_NAME( REG3_TURB1 ,reg3_turb1 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREG3_TURB1 

#ifndef GUARDREG3_TURB2 
#define GUARDREG3_TURB2 
// Prototype for Fortran procedure REG3_TURB2 ...
//
void FORTRAN_NAME( REG3_TURB2 ,reg3_turb2 )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Phi)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_BOX(box) );

#define FORT_REG3_TURB2 FORTRAN_NAME( inlineREG3_TURB2, inlineREG3_TURB2)
#define FORTNT_REG3_TURB2 FORTRAN_NAME( REG3_TURB2, reg3_turb2)

inline void FORTRAN_NAME(inlineREG3_TURB2, inlineREG3_TURB2)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Phi)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REG3_TURB2");
 FORTRAN_NAME( REG3_TURB2 ,reg3_turb2 )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Phi)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREG3_TURB2 

#ifndef GUARDLISMBC_TURB 
#define GUARDLISMBC_TURB 
// Prototype for Fortran procedure LISMBC_TURB ...
//
void FORTRAN_NAME( LISMBC_TURB ,lismbc_turb )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_LISMBC_TURB FORTRAN_NAME( inlineLISMBC_TURB, inlineLISMBC_TURB)
#define FORTNT_LISMBC_TURB FORTRAN_NAME( LISMBC_TURB, lismbc_turb)

inline void FORTRAN_NAME(inlineLISMBC_TURB, inlineLISMBC_TURB)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LISMBC_TURB");
 FORTRAN_NAME( LISMBC_TURB ,lismbc_turb )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLISMBC_TURB 

}

#endif
