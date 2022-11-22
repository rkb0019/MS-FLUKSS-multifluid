#ifndef _DIVCLEANINGF_F_H_
#define _DIVCLEANINGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRECALCULATEB 
#define GUARDRECALCULATEB 
// Prototype for Fortran procedure recalculateB ...
//
void FORTRAN_NAME( RECALCULATEB ,recalculateb )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UOld)
      ,CHFp_CONST_FRA1(EX)
      ,CHFp_CONST_FRA1(EY)
      ,CHFp_CONST_FRA1(EZ)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(BBox) );

#define FORT_RECALCULATEB FORTRAN_NAME( inlineRECALCULATEB, inlineRECALCULATEB)
#define FORTNT_RECALCULATEB FORTRAN_NAME( RECALCULATEB, recalculateb)

inline void FORTRAN_NAME(inlineRECALCULATEB, inlineRECALCULATEB)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UOld)
      ,CHFp_CONST_FRA1(EX)
      ,CHFp_CONST_FRA1(EY)
      ,CHFp_CONST_FRA1(EZ)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(BBox) )
{
 CH_TIMELEAF("FORT_RECALCULATEB");
 FORTRAN_NAME( RECALCULATEB ,recalculateb )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(UOld)
      ,CHFt_CONST_FRA1(EX)
      ,CHFt_CONST_FRA1(EY)
      ,CHFt_CONST_FRA1(EZ)
      ,CHFt_CONST_INT(iBX)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(BBox) );
}
#endif  // GUARDRECALCULATEB 

#ifndef GUARDELECTRICFIELDBS 
#define GUARDELECTRICFIELDBS 
// Prototype for Fortran procedure electricFieldBS ...
//
void FORTRAN_NAME( ELECTRICFIELDBS ,electricfieldbs )(
      CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) );

#define FORT_ELECTRICFIELDBS FORTRAN_NAME( inlineELECTRICFIELDBS, inlineELECTRICFIELDBS)
#define FORTNT_ELECTRICFIELDBS FORTRAN_NAME( ELECTRICFIELDBS, electricfieldbs)

inline void FORTRAN_NAME(inlineELECTRICFIELDBS, inlineELECTRICFIELDBS)(
      CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) )
{
 CH_TIMELEAF("FORT_ELECTRICFIELDBS");
 FORTRAN_NAME( ELECTRICFIELDBS ,electricfieldbs )(
      CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA(G)
      ,CHFt_CONST_FRA(H)
      ,CHFt_FRA1(EX)
      ,CHFt_FRA1(EY)
      ,CHFt_FRA1(EZ)
      ,CHFt_CONST_INT(iBX)
      ,CHFt_BOX(BBox) );
}
#endif  // GUARDELECTRICFIELDBS 

#ifndef GUARDELECTRICFIELDGS0 
#define GUARDELECTRICFIELDGS0 
// Prototype for Fortran procedure electricFieldGS0 ...
//
void FORTRAN_NAME( ELECTRICFIELDGS0 ,electricfieldgs0 )(
      CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iUX)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) );

#define FORT_ELECTRICFIELDGS0 FORTRAN_NAME( inlineELECTRICFIELDGS0, inlineELECTRICFIELDGS0)
#define FORTNT_ELECTRICFIELDGS0 FORTRAN_NAME( ELECTRICFIELDGS0, electricfieldgs0)

inline void FORTRAN_NAME(inlineELECTRICFIELDGS0, inlineELECTRICFIELDGS0)(
      CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iUX)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) )
{
 CH_TIMELEAF("FORT_ELECTRICFIELDGS0");
 FORTRAN_NAME( ELECTRICFIELDGS0 ,electricfieldgs0 )(
      CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA(G)
      ,CHFt_CONST_FRA(H)
      ,CHFt_FRA1(EX)
      ,CHFt_FRA1(EY)
      ,CHFt_FRA1(EZ)
      ,CHFt_CONST_INT(iUX)
      ,CHFt_CONST_INT(iBX)
      ,CHFt_BOX(BBox) );
}
#endif  // GUARDELECTRICFIELDGS0 

#ifndef GUARDELECTRICFIELDGS1 
#define GUARDELECTRICFIELDGS1 
// Prototype for Fortran procedure electricFieldGS1 ...
//
void FORTRAN_NAME( ELECTRICFIELDGS1 ,electricfieldgs1 )(
      CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iUX)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) );

#define FORT_ELECTRICFIELDGS1 FORTRAN_NAME( inlineELECTRICFIELDGS1, inlineELECTRICFIELDGS1)
#define FORTNT_ELECTRICFIELDGS1 FORTRAN_NAME( ELECTRICFIELDGS1, electricfieldgs1)

inline void FORTRAN_NAME(inlineELECTRICFIELDGS1, inlineELECTRICFIELDGS1)(
      CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA(F)
      ,CHFp_CONST_FRA(G)
      ,CHFp_CONST_FRA(H)
      ,CHFp_FRA1(EX)
      ,CHFp_FRA1(EY)
      ,CHFp_FRA1(EZ)
      ,CHFp_CONST_INT(iUX)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) )
{
 CH_TIMELEAF("FORT_ELECTRICFIELDGS1");
 FORTRAN_NAME( ELECTRICFIELDGS1 ,electricfieldgs1 )(
      CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA(F)
      ,CHFt_CONST_FRA(G)
      ,CHFt_CONST_FRA(H)
      ,CHFt_FRA1(EX)
      ,CHFt_FRA1(EY)
      ,CHFt_FRA1(EZ)
      ,CHFt_CONST_INT(iUX)
      ,CHFt_CONST_INT(iBX)
      ,CHFt_BOX(BBox) );
}
#endif  // GUARDELECTRICFIELDGS1 

#ifndef GUARDCHECKDIVB 
#define GUARDCHECKDIVB 
// Prototype for Fortran procedure checkDivB ...
//
void FORTRAN_NAME( CHECKDIVB ,checkdivb )(
      CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) );

#define FORT_CHECKDIVB FORTRAN_NAME( inlineCHECKDIVB, inlineCHECKDIVB)
#define FORTNT_CHECKDIVB FORTRAN_NAME( CHECKDIVB, checkdivb)

inline void FORTRAN_NAME(inlineCHECKDIVB, inlineCHECKDIVB)(
      CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iBX)
      ,CHFp_BOX(BBox) )
{
 CH_TIMELEAF("FORT_CHECKDIVB");
 FORTRAN_NAME( CHECKDIVB ,checkdivb )(
      CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iBX)
      ,CHFt_BOX(BBox) );
}
#endif  // GUARDCHECKDIVB 

}

#endif
