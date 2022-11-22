#ifndef _OBLIQUESHOCKTUBEF_F_H_
#define _OBLIQUESHOCKTUBEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETOBLIQUESHOCKTUBE 
#define GUARDSETOBLIQUESHOCKTUBE 
// Prototype for Fortran procedure SETOBLIQUESHOCKTUBE ...
//
void FORTRAN_NAME( SETOBLIQUESHOCKTUBE ,setobliqueshocktube )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auPressureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auYs) );

#define FORT_SETOBLIQUESHOCKTUBE FORTRAN_NAME( inlineSETOBLIQUESHOCKTUBE, inlineSETOBLIQUESHOCKTUBE)
#define FORTNT_SETOBLIQUESHOCKTUBE FORTRAN_NAME( SETOBLIQUESHOCKTUBE, setobliqueshocktube)

inline void FORTRAN_NAME(inlineSETOBLIQUESHOCKTUBE, inlineSETOBLIQUESHOCKTUBE)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auPressureL)
      ,CHFp_CONST_REAL(auPressureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auYs) )
{
 CH_TIMELEAF("FORT_SETOBLIQUESHOCKTUBE");
 FORTRAN_NAME( SETOBLIQUESHOCKTUBE ,setobliqueshocktube )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auDensityL)
      ,CHFt_CONST_REAL(auDensityR)
      ,CHFt_CONST_REAL(auPressureL)
      ,CHFt_CONST_REAL(auPressureR)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auVelxR)
      ,CHFt_CONST_REAL(auVelyL)
      ,CHFt_CONST_REAL(auVelyR)
      ,CHFt_CONST_REAL(auVelzL)
      ,CHFt_CONST_REAL(auVelzR)
      ,CHFt_CONST_REAL(auBxL)
      ,CHFt_CONST_REAL(auBxR)
      ,CHFt_CONST_REAL(auByL)
      ,CHFt_CONST_REAL(auByR)
      ,CHFt_CONST_REAL(auBzL)
      ,CHFt_CONST_REAL(auBzR)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_REAL(auYs) );
}
#endif  // GUARDSETOBLIQUESHOCKTUBE 

#ifndef GUARDOBLIQUESHOCKTUBEINIT 
#define GUARDOBLIQUESHOCKTUBEINIT 
// Prototype for Fortran procedure OBLIQUESHOCKTUBEINIT ...
//
void FORTRAN_NAME( OBLIQUESHOCKTUBEINIT ,obliqueshocktubeinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(tanangle)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_OBLIQUESHOCKTUBEINIT FORTRAN_NAME( inlineOBLIQUESHOCKTUBEINIT, inlineOBLIQUESHOCKTUBEINIT)
#define FORTNT_OBLIQUESHOCKTUBEINIT FORTRAN_NAME( OBLIQUESHOCKTUBEINIT, obliqueshocktubeinit)

inline void FORTRAN_NAME(inlineOBLIQUESHOCKTUBEINIT, inlineOBLIQUESHOCKTUBEINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(tanangle)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OBLIQUESHOCKTUBEINIT");
 FORTRAN_NAME( OBLIQUESHOCKTUBEINIT ,obliqueshocktubeinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(tanangle)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOBLIQUESHOCKTUBEINIT 

#ifndef GUARDOBLIQUESHOCKTUBEGS 
#define GUARDOBLIQUESHOCKTUBEGS 
// Prototype for Fortran procedure OBLIQUESHOCKTUBEGS ...
//
void FORTRAN_NAME( OBLIQUESHOCKTUBEGS ,obliqueshocktubegs )(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(tanangle)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_OBLIQUESHOCKTUBEGS FORTRAN_NAME( inlineOBLIQUESHOCKTUBEGS, inlineOBLIQUESHOCKTUBEGS)
#define FORTNT_OBLIQUESHOCKTUBEGS FORTRAN_NAME( OBLIQUESHOCKTUBEGS, obliqueshocktubegs)

inline void FORTRAN_NAME(inlineOBLIQUESHOCKTUBEGS, inlineOBLIQUESHOCKTUBEGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(tanangle)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OBLIQUESHOCKTUBEGS");
 FORTRAN_NAME( OBLIQUESHOCKTUBEGS ,obliqueshocktubegs )(
      CHFt_FRA(W)
      ,CHFt_CONST_REAL(tanangle)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOBLIQUESHOCKTUBEGS 

}

#endif
