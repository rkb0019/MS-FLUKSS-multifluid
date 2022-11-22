#ifndef _RAMPF_F_H_
#define _RAMPF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETRAMP 
#define GUARDSETRAMP 
// Prototype for Fortran procedure SETRAMP ...
//
void FORTRAN_NAME( SETRAMP ,setramp )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aualpha)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(auxcorner) );

#define FORT_SETRAMP FORTRAN_NAME( inlineSETRAMP, inlineSETRAMP)
#define FORTNT_SETRAMP FORTRAN_NAME( SETRAMP, setramp)

inline void FORTRAN_NAME(inlineSETRAMP, inlineSETRAMP)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aualpha)
      ,CHFp_CONST_REAL(aums)
      ,CHFp_CONST_REAL(auxcorner) )
{
 CH_TIMELEAF("FORT_SETRAMP");
 FORTRAN_NAME( SETRAMP ,setramp )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aualpha)
      ,CHFt_CONST_REAL(aums)
      ,CHFt_CONST_REAL(auxcorner) );
}
#endif  // GUARDSETRAMP 

#ifndef GUARDRAMPINIT 
#define GUARDRAMPINIT 
// Prototype for Fortran procedure RAMPINIT ...
//
void FORTRAN_NAME( RAMPINIT ,rampinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_RAMPINIT FORTRAN_NAME( inlineRAMPINIT, inlineRAMPINIT)
#define FORTNT_RAMPINIT FORTRAN_NAME( RAMPINIT, rampinit)

inline void FORTRAN_NAME(inlineRAMPINIT, inlineRAMPINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RAMPINIT");
 FORTRAN_NAME( RAMPINIT ,rampinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRAMPINIT 

#ifndef GUARDRAMPBC 
#define GUARDRAMPBC 
// Prototype for Fortran procedure RAMPBC ...
//
void FORTRAN_NAME( RAMPBC ,rampbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_RAMPBC FORTRAN_NAME( inlineRAMPBC, inlineRAMPBC)
#define FORTNT_RAMPBC FORTRAN_NAME( RAMPBC, rampbc)

inline void FORTRAN_NAME(inlineRAMPBC, inlineRAMPBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RAMPBC");
 FORTRAN_NAME( RAMPBC ,rampbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRAMPBC 

#ifndef GUARDRAMPSOLIDBC 
#define GUARDRAMPSOLIDBC 
// Prototype for Fortran procedure RAMPSOLIDBC ...
//
void FORTRAN_NAME( RAMPSOLIDBC ,rampsolidbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_RAMPSOLIDBC FORTRAN_NAME( inlineRAMPSOLIDBC, inlineRAMPSOLIDBC)
#define FORTNT_RAMPSOLIDBC FORTRAN_NAME( RAMPSOLIDBC, rampsolidbc)

inline void FORTRAN_NAME(inlineRAMPSOLIDBC, inlineRAMPSOLIDBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Wextrap)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RAMPSOLIDBC");
 FORTRAN_NAME( RAMPSOLIDBC ,rampsolidbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Wextrap)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRAMPSOLIDBC 

#ifndef GUARDRAMPSLOPEBCS 
#define GUARDRAMPSLOPEBCS 
// Prototype for Fortran procedure RAMPSLOPEBCS ...
//
void FORTRAN_NAME( RAMPSLOPEBCS ,rampslopebcs )(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) );

#define FORT_RAMPSLOPEBCS FORTRAN_NAME( inlineRAMPSLOPEBCS, inlineRAMPSLOPEBCS)
#define FORTNT_RAMPSLOPEBCS FORTRAN_NAME( RAMPSLOPEBCS, rampslopebcs)

inline void FORTRAN_NAME(inlineRAMPSLOPEBCS, inlineRAMPSLOPEBCS)(
      CHFp_FRA(dW)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi) )
{
 CH_TIMELEAF("FORT_RAMPSLOPEBCS");
 FORTRAN_NAME( RAMPSLOPEBCS ,rampslopebcs )(
      CHFt_FRA(dW)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi) );
}
#endif  // GUARDRAMPSLOPEBCS 

#ifndef GUARDRAMPGS 
#define GUARDRAMPGS 
// Prototype for Fortran procedure RAMPGS ...
//
void FORTRAN_NAME( RAMPGS ,rampgs )(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_RAMPGS FORTRAN_NAME( inlineRAMPGS, inlineRAMPGS)
#define FORTNT_RAMPGS FORTRAN_NAME( RAMPGS, rampgs)

inline void FORTRAN_NAME(inlineRAMPGS, inlineRAMPGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RAMPGS");
 FORTRAN_NAME( RAMPGS ,rampgs )(
      CHFt_FRA(W)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRAMPGS 

}

#endif
