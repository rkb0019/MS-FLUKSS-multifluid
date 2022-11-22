#ifndef _AVERAGEEDGEF_F_H_
#define _AVERAGEEDGEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGEEDGE 
#define GUARDAVERAGEEDGE 
// Prototype for Fortran procedure AVERAGEEDGE ...
//
void FORTRAN_NAME( AVERAGEEDGE ,averageedge )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_BOX(refBox) );

#define FORT_AVERAGEEDGE FORTRAN_NAME( inlineAVERAGEEDGE, inlineAVERAGEEDGE)
#define FORTNT_AVERAGEEDGE FORTRAN_NAME( AVERAGEEDGE, averageedge)

inline void FORTRAN_NAME(inlineAVERAGEEDGE, inlineAVERAGEEDGE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_AVERAGEEDGE");
 FORTRAN_NAME( AVERAGEEDGE ,averageedge )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(nRef)
      ,CHFt_CONST_REAL(dWeight)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDAVERAGEEDGE 

#ifndef GUARDAVERAGEEDGEHARMONIC 
#define GUARDAVERAGEEDGEHARMONIC 
// Prototype for Fortran procedure AVERAGEEDGEHARMONIC ...
//
void FORTRAN_NAME( AVERAGEEDGEHARMONIC ,averageedgeharmonic )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_BOX(refBox) );

#define FORT_AVERAGEEDGEHARMONIC FORTRAN_NAME( inlineAVERAGEEDGEHARMONIC, inlineAVERAGEEDGEHARMONIC)
#define FORTNT_AVERAGEEDGEHARMONIC FORTRAN_NAME( AVERAGEEDGEHARMONIC, averageedgeharmonic)

inline void FORTRAN_NAME(inlineAVERAGEEDGEHARMONIC, inlineAVERAGEEDGEHARMONIC)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_INT(nRef)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_BOX(refBox) )
{
 CH_TIMELEAF("FORT_AVERAGEEDGEHARMONIC");
 FORTRAN_NAME( AVERAGEEDGEHARMONIC ,averageedgeharmonic )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_INT(nRef)
      ,CHFt_CONST_REAL(dWeight)
      ,CHFt_BOX(refBox) );
}
#endif  // GUARDAVERAGEEDGEHARMONIC 

#ifndef GUARDFINETOCOARSE 
#define GUARDFINETOCOARSE 
// Prototype for Fortran procedure FineToCoarse ...
//
void FORTRAN_NAME( FINETOCOARSE ,finetocoarse )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_CONST_INT(nRef) );

#define FORT_FINETOCOARSE FORTRAN_NAME( inlineFINETOCOARSE, inlineFINETOCOARSE)
#define FORTNT_FINETOCOARSE FORTRAN_NAME( FINETOCOARSE, finetocoarse)

inline void FORTRAN_NAME(inlineFINETOCOARSE, inlineFINETOCOARSE)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_BOX(crseBox)
      ,CHFp_CONST_REAL(dWeight)
      ,CHFp_CONST_INT(nRef) )
{
 CH_TIMELEAF("FORT_FINETOCOARSE");
 FORTRAN_NAME( FINETOCOARSE ,finetocoarse )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_BOX(crseBox)
      ,CHFt_CONST_REAL(dWeight)
      ,CHFt_CONST_INT(nRef) );
}
#endif  // GUARDFINETOCOARSE 

}

#endif
