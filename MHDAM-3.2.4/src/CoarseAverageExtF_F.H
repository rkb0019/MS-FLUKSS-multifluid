#ifndef _COARSEAVERAGEEXTF_F_H_
#define _COARSEAVERAGEEXTF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDAVERAGEWITHVOLUMES 
#define GUARDAVERAGEWITHVOLUMES 
// Prototype for Fortran procedure averagewithvolumes ...
//
void FORTRAN_NAME( AVERAGEWITHVOLUMES ,averagewithvolumes )(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(volf)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) );

#define FORT_AVERAGEWITHVOLUMES FORTRAN_NAME( inlineAVERAGEWITHVOLUMES, inlineAVERAGEWITHVOLUMES)
#define FORTNT_AVERAGEWITHVOLUMES FORTRAN_NAME( AVERAGEWITHVOLUMES, averagewithvolumes)

inline void FORTRAN_NAME(inlineAVERAGEWITHVOLUMES, inlineAVERAGEWITHVOLUMES)(
      CHFp_FRA(coarse)
      ,CHFp_CONST_FRA(fine)
      ,CHFp_CONST_FRA1(volf)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_AVERAGEWITHVOLUMES");
 FORTRAN_NAME( AVERAGEWITHVOLUMES ,averagewithvolumes )(
      CHFt_FRA(coarse)
      ,CHFt_CONST_FRA(fine)
      ,CHFt_CONST_FRA1(volf)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDAVERAGEWITHVOLUMES 

}

#endif
