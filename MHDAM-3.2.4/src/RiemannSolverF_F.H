#ifndef _RIEMANNSOLVERF_F_H_
#define _RIEMANNSOLVERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETRSCONST 
#define GUARDSETRSCONST 
// Prototype for Fortran procedure SETRSCONST ...
//
void FORTRAN_NAME( SETRSCONST ,setrsconst )(
      CHFp_CONST_REAL(auDel2)
      ,CHFp_CONST_REAL(auSmallB)
      ,CHFp_CONST_INT(auAver)
      ,CHFp_CONST_INT(auLF) );

#define FORT_SETRSCONST FORTRAN_NAME( inlineSETRSCONST, inlineSETRSCONST)
#define FORTNT_SETRSCONST FORTRAN_NAME( SETRSCONST, setrsconst)

inline void FORTRAN_NAME(inlineSETRSCONST, inlineSETRSCONST)(
      CHFp_CONST_REAL(auDel2)
      ,CHFp_CONST_REAL(auSmallB)
      ,CHFp_CONST_INT(auAver)
      ,CHFp_CONST_INT(auLF) )
{
 CH_TIMELEAF("FORT_SETRSCONST");
 FORTRAN_NAME( SETRSCONST ,setrsconst )(
      CHFt_CONST_REAL(auDel2)
      ,CHFt_CONST_REAL(auSmallB)
      ,CHFt_CONST_INT(auAver)
      ,CHFt_CONST_INT(auLF) );
}
#endif  // GUARDSETRSCONST 

}

#endif
