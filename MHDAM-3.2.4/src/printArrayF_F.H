#ifndef _PRINTARRAYF_F_H_
#define _PRINTARRAYF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPRINTARRAY 
#define GUARDPRINTARRAY 
// Prototype for Fortran procedure PRINTARRAY ...
//
void FORTRAN_NAME( PRINTARRAY ,printarray )(
      CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(index)
      ,CHFp_CONST_INT(flag)
      ,CHFp_BOX(box) );

#define FORT_PRINTARRAY FORTRAN_NAME( inlinePRINTARRAY, inlinePRINTARRAY)
#define FORTNT_PRINTARRAY FORTRAN_NAME( PRINTARRAY, printarray)

inline void FORTRAN_NAME(inlinePRINTARRAY, inlinePRINTARRAY)(
      CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(index)
      ,CHFp_CONST_INT(flag)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRINTARRAY");
 FORTRAN_NAME( PRINTARRAY ,printarray )(
      CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(index)
      ,CHFt_CONST_INT(flag)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRINTARRAY 

}

#endif
