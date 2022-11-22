#ifndef _HELIOTILTF_F_H_
#define _HELIOTILTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDHELIOTILTBC 
#define GUARDHELIOTILTBC 
// Prototype for Fortran procedure HELIOTILTBC ...
//
void FORTRAN_NAME( HELIOTILTBC ,heliotiltbc )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_HELIOTILTBC FORTRAN_NAME( inlineHELIOTILTBC, inlineHELIOTILTBC)
#define FORTNT_HELIOTILTBC FORTRAN_NAME( HELIOTILTBC, heliotiltbc)

inline void FORTRAN_NAME(inlineHELIOTILTBC, inlineHELIOTILTBC)(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOTILTBC");
 FORTRAN_NAME( HELIOTILTBC ,heliotiltbc )(
      CHFt_FRA(F)
      ,CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOTILTBC 

#ifndef GUARDHELIOTILTGS 
#define GUARDHELIOTILTGS 
// Prototype for Fortran procedure HELIOTILTGS ...
//
void FORTRAN_NAME( HELIOTILTGS ,heliotiltgs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_BOX(box) );

#define FORT_HELIOTILTGS FORTRAN_NAME( inlineHELIOTILTGS, inlineHELIOTILTGS)
#define FORTNT_HELIOTILTGS FORTRAN_NAME( HELIOTILTGS, heliotiltgs)

inline void FORTRAN_NAME(inlineHELIOTILTGS, inlineHELIOTILTGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOTILTGS");
 FORTRAN_NAME( HELIOTILTGS ,heliotiltgs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOTILTGS 

#ifndef GUARDHELIOTILTGSSPHERICAL 
#define GUARDHELIOTILTGSSPHERICAL 
// Prototype for Fortran procedure HELIOTILTGSSPHERICAL ...
//
void FORTRAN_NAME( HELIOTILTGSSPHERICAL ,heliotiltgsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOTILTGSSPHERICAL FORTRAN_NAME( inlineHELIOTILTGSSPHERICAL, inlineHELIOTILTGSSPHERICAL)
#define FORTNT_HELIOTILTGSSPHERICAL FORTRAN_NAME( HELIOTILTGSSPHERICAL, heliotiltgsspherical)

inline void FORTRAN_NAME(inlineHELIOTILTGSSPHERICAL, inlineHELIOTILTGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOTILTGSSPHERICAL");
 FORTRAN_NAME( HELIOTILTGSSPHERICAL ,heliotiltgsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOTILTGSSPHERICAL 

#ifndef GUARDHELIOTILTINIT 
#define GUARDHELIOTILTINIT 
// Prototype for Fortran procedure HELIOTILTINIT ...
//
void FORTRAN_NAME( HELIOTILTINIT ,heliotiltinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_BOX(box) );

#define FORT_HELIOTILTINIT FORTRAN_NAME( inlineHELIOTILTINIT, inlineHELIOTILTINIT)
#define FORTNT_HELIOTILTINIT FORTRAN_NAME( HELIOTILTINIT, heliotiltinit)

inline void FORTRAN_NAME(inlineHELIOTILTINIT, inlineHELIOTILTINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOTILTINIT");
 FORTRAN_NAME( HELIOTILTINIT ,heliotiltinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOTILTINIT 

#ifndef GUARDHELIOTILTINITSPHERICAL 
#define GUARDHELIOTILTINITSPHERICAL 
// Prototype for Fortran procedure HELIOTILTINITSPHERICAL ...
//
void FORTRAN_NAME( HELIOTILTINITSPHERICAL ,heliotiltinitspherical )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOTILTINITSPHERICAL FORTRAN_NAME( inlineHELIOTILTINITSPHERICAL, inlineHELIOTILTINITSPHERICAL)
#define FORTNT_HELIOTILTINITSPHERICAL FORTRAN_NAME( HELIOTILTINITSPHERICAL, heliotiltinitspherical)

inline void FORTRAN_NAME(inlineHELIOTILTINITSPHERICAL, inlineHELIOTILTINITSPHERICAL)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOTILTINITSPHERICAL");
 FORTRAN_NAME( HELIOTILTINITSPHERICAL ,heliotiltinitspherical )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOTILTINITSPHERICAL 

}

#endif
