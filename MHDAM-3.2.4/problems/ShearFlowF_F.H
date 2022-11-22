#ifndef _SHEARFLOWF_F_H_
#define _SHEARFLOWF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETSHEARFLOW 
#define GUARDSETSHEARFLOW 
// Prototype for Fortran procedure SETSHEARFLOW ...
//
void FORTRAN_NAME( SETSHEARFLOW ,setshearflow )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auN1)
      ,CHFp_CONST_REAL(auV1)
      ,CHFp_CONST_REAL(auT1)
      ,CHFp_CONST_REAL(auN2)
      ,CHFp_CONST_REAL(auV2)
      ,CHFp_CONST_REAL(auT2)
      ,CHFp_CONST_REAL(aufs_a) );

#define FORT_SETSHEARFLOW FORTRAN_NAME( inlineSETSHEARFLOW, inlineSETSHEARFLOW)
#define FORTNT_SETSHEARFLOW FORTRAN_NAME( SETSHEARFLOW, setshearflow)

inline void FORTRAN_NAME(inlineSETSHEARFLOW, inlineSETSHEARFLOW)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auN1)
      ,CHFp_CONST_REAL(auV1)
      ,CHFp_CONST_REAL(auT1)
      ,CHFp_CONST_REAL(auN2)
      ,CHFp_CONST_REAL(auV2)
      ,CHFp_CONST_REAL(auT2)
      ,CHFp_CONST_REAL(aufs_a) )
{
 CH_TIMELEAF("FORT_SETSHEARFLOW");
 FORTRAN_NAME( SETSHEARFLOW ,setshearflow )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auN1)
      ,CHFt_CONST_REAL(auV1)
      ,CHFt_CONST_REAL(auT1)
      ,CHFt_CONST_REAL(auN2)
      ,CHFt_CONST_REAL(auV2)
      ,CHFt_CONST_REAL(auT2)
      ,CHFt_CONST_REAL(aufs_a) );
}
#endif  // GUARDSETSHEARFLOW 

#ifndef GUARDSHEARFLOWINIT 
#define GUARDSHEARFLOWINIT 
// Prototype for Fortran procedure SHEARFLOWINIT ...
//
void FORTRAN_NAME( SHEARFLOWINIT ,shearflowinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(jmiddle)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SHEARFLOWINIT FORTRAN_NAME( inlineSHEARFLOWINIT, inlineSHEARFLOWINIT)
#define FORTNT_SHEARFLOWINIT FORTRAN_NAME( SHEARFLOWINIT, shearflowinit)

inline void FORTRAN_NAME(inlineSHEARFLOWINIT, inlineSHEARFLOWINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(jmiddle)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SHEARFLOWINIT");
 FORTRAN_NAME( SHEARFLOWINIT ,shearflowinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(jmiddle)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSHEARFLOWINIT 

#ifndef GUARDSHEARFLOWGS 
#define GUARDSHEARFLOWGS 
// Prototype for Fortran procedure SHEARFLOWGS ...
//
void FORTRAN_NAME( SHEARFLOWGS ,shearflowgs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(jmiddle)
      ,CHFp_BOX(box) );

#define FORT_SHEARFLOWGS FORTRAN_NAME( inlineSHEARFLOWGS, inlineSHEARFLOWGS)
#define FORTNT_SHEARFLOWGS FORTRAN_NAME( SHEARFLOWGS, shearflowgs)

inline void FORTRAN_NAME(inlineSHEARFLOWGS, inlineSHEARFLOWGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(jmiddle)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SHEARFLOWGS");
 FORTRAN_NAME( SHEARFLOWGS ,shearflowgs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(jmiddle)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSHEARFLOWGS 

}

#endif
