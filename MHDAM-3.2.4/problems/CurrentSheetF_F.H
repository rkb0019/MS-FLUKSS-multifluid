#ifndef _CURRENTSHEETF_F_H_
#define _CURRENTSHEETF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCSHEET 
#define GUARDSETCSHEET 
// Prototype for Fortran procedure SETCSHEET ...
//
void FORTRAN_NAME( SETCSHEET ,setcsheet )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auCSheetU0)
      ,CHFp_CONST_REAL(auCSheetB0)
      ,CHFp_CONST_REAL(auCSheetP0) );

#define FORT_SETCSHEET FORTRAN_NAME( inlineSETCSHEET, inlineSETCSHEET)
#define FORTNT_SETCSHEET FORTRAN_NAME( SETCSHEET, setcsheet)

inline void FORTRAN_NAME(inlineSETCSHEET, inlineSETCSHEET)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auCSheetU0)
      ,CHFp_CONST_REAL(auCSheetB0)
      ,CHFp_CONST_REAL(auCSheetP0) )
{
 CH_TIMELEAF("FORT_SETCSHEET");
 FORTRAN_NAME( SETCSHEET ,setcsheet )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auCSheetU0)
      ,CHFt_CONST_REAL(auCSheetB0)
      ,CHFt_CONST_REAL(auCSheetP0) );
}
#endif  // GUARDSETCSHEET 

#ifndef GUARDCSHEETINIT 
#define GUARDCSHEETINIT 
// Prototype for Fortran procedure CSHEETINIT ...
//
void FORTRAN_NAME( CSHEETINIT ,csheetinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_CSHEETINIT FORTRAN_NAME( inlineCSHEETINIT, inlineCSHEETINIT)
#define FORTNT_CSHEETINIT FORTRAN_NAME( CSHEETINIT, csheetinit)

inline void FORTRAN_NAME(inlineCSHEETINIT, inlineCSHEETINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CSHEETINIT");
 FORTRAN_NAME( CSHEETINIT ,csheetinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCSHEETINIT 

}

#endif
