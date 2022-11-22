#ifndef _ORSZAGF_F_H_
#define _ORSZAGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETORSZAG 
#define GUARDSETORSZAG 
// Prototype for Fortran procedure SETORSZAG ...
//
void FORTRAN_NAME( SETORSZAG ,setorszag )(
      CHFp_CONST_REAL(augamma) );

#define FORT_SETORSZAG FORTRAN_NAME( inlineSETORSZAG, inlineSETORSZAG)
#define FORTNT_SETORSZAG FORTRAN_NAME( SETORSZAG, setorszag)

inline void FORTRAN_NAME(inlineSETORSZAG, inlineSETORSZAG)(
      CHFp_CONST_REAL(augamma) )
{
 CH_TIMELEAF("FORT_SETORSZAG");
 FORTRAN_NAME( SETORSZAG ,setorszag )(
      CHFt_CONST_REAL(augamma) );
}
#endif  // GUARDSETORSZAG 

#ifndef GUARDORSZAGINIT 
#define GUARDORSZAGINIT 
// Prototype for Fortran procedure ORSZAGINIT ...
//
void FORTRAN_NAME( ORSZAGINIT ,orszaginit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_ORSZAGINIT FORTRAN_NAME( inlineORSZAGINIT, inlineORSZAGINIT)
#define FORTNT_ORSZAGINIT FORTRAN_NAME( ORSZAGINIT, orszaginit)

inline void FORTRAN_NAME(inlineORSZAGINIT, inlineORSZAGINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ORSZAGINIT");
 FORTRAN_NAME( ORSZAGINIT ,orszaginit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDORSZAGINIT 

}

#endif
