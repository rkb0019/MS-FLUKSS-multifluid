#ifndef _BLASTF_F_H_
#define _BLASTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETBLAST 
#define GUARDSETBLAST 
// Prototype for Fortran procedure SETBLAST ...
//
void FORTRAN_NAME( SETBLAST ,setblast )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auBlastp0)
      ,CHFp_CONST_REAL(auBlastp1)
      ,CHFp_CONST_REAL(auBlastXc)
      ,CHFp_CONST_REAL(auBlastYc)
      ,CHFp_CONST_REAL(auBlastZc)
      ,CHFp_CONST_REAL(auBlastBx)
      ,CHFp_CONST_REAL(auBlastBy)
      ,CHFp_CONST_REAL(auBlastr0) );

#define FORT_SETBLAST FORTRAN_NAME( inlineSETBLAST, inlineSETBLAST)
#define FORTNT_SETBLAST FORTRAN_NAME( SETBLAST, setblast)

inline void FORTRAN_NAME(inlineSETBLAST, inlineSETBLAST)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auBlastp0)
      ,CHFp_CONST_REAL(auBlastp1)
      ,CHFp_CONST_REAL(auBlastXc)
      ,CHFp_CONST_REAL(auBlastYc)
      ,CHFp_CONST_REAL(auBlastZc)
      ,CHFp_CONST_REAL(auBlastBx)
      ,CHFp_CONST_REAL(auBlastBy)
      ,CHFp_CONST_REAL(auBlastr0) )
{
 CH_TIMELEAF("FORT_SETBLAST");
 FORTRAN_NAME( SETBLAST ,setblast )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auBlastp0)
      ,CHFt_CONST_REAL(auBlastp1)
      ,CHFt_CONST_REAL(auBlastXc)
      ,CHFt_CONST_REAL(auBlastYc)
      ,CHFt_CONST_REAL(auBlastZc)
      ,CHFt_CONST_REAL(auBlastBx)
      ,CHFt_CONST_REAL(auBlastBy)
      ,CHFt_CONST_REAL(auBlastr0) );
}
#endif  // GUARDSETBLAST 

#ifndef GUARDBLASTINIT 
#define GUARDBLASTINIT 
// Prototype for Fortran procedure BLASTINIT ...
//
void FORTRAN_NAME( BLASTINIT ,blastinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_BLASTINIT FORTRAN_NAME( inlineBLASTINIT, inlineBLASTINIT)
#define FORTNT_BLASTINIT FORTRAN_NAME( BLASTINIT, blastinit)

inline void FORTRAN_NAME(inlineBLASTINIT, inlineBLASTINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_BLASTINIT");
 FORTRAN_NAME( BLASTINIT ,blastinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDBLASTINIT 

}

#endif
