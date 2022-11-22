#ifndef _KELVINF_F_H_
#define _KELVINF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETKELVIN 
#define GUARDSETKELVIN 
// Prototype for Fortran procedure SETKELVIN ...
//
void FORTRAN_NAME( SETKELVIN ,setkelvin )(
      CHFp_CONST_REAL(augamma) );

#define FORT_SETKELVIN FORTRAN_NAME( inlineSETKELVIN, inlineSETKELVIN)
#define FORTNT_SETKELVIN FORTRAN_NAME( SETKELVIN, setkelvin)

inline void FORTRAN_NAME(inlineSETKELVIN, inlineSETKELVIN)(
      CHFp_CONST_REAL(augamma) )
{
 CH_TIMELEAF("FORT_SETKELVIN");
 FORTRAN_NAME( SETKELVIN ,setkelvin )(
      CHFt_CONST_REAL(augamma) );
}
#endif  // GUARDSETKELVIN 

#ifndef GUARDKELVININIT 
#define GUARDKELVININIT 
// Prototype for Fortran procedure KELVININIT ...
//
void FORTRAN_NAME( KELVININIT ,kelvininit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_KELVININIT FORTRAN_NAME( inlineKELVININIT, inlineKELVININIT)
#define FORTNT_KELVININIT FORTRAN_NAME( KELVININIT, kelvininit)

inline void FORTRAN_NAME(inlineKELVININIT, inlineKELVININIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_KELVININIT");
 FORTRAN_NAME( KELVININIT ,kelvininit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDKELVININIT 

}

#endif
