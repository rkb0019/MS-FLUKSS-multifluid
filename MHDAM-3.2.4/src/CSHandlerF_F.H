#ifndef _CSHANDLERF_F_H_
#define _CSHANDLERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMULTIPLYFLUXESBYAREA 
#define GUARDMULTIPLYFLUXESBYAREA 
// Prototype for Fortran procedure MULTIPLYFLUXESBYAREA ...
//
void FORTRAN_NAME( MULTIPLYFLUXESBYAREA ,multiplyfluxesbyarea )(
      CHFp_FRA(F)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_BOX(box) );

#define FORT_MULTIPLYFLUXESBYAREA FORTRAN_NAME( inlineMULTIPLYFLUXESBYAREA, inlineMULTIPLYFLUXESBYAREA)
#define FORTNT_MULTIPLYFLUXESBYAREA FORTRAN_NAME( MULTIPLYFLUXESBYAREA, multiplyfluxesbyarea)

inline void FORTRAN_NAME(inlineMULTIPLYFLUXESBYAREA, inlineMULTIPLYFLUXESBYAREA)(
      CHFp_FRA(F)
      ,CHFp_CONST_INT(iBgn)
      ,CHFp_CONST_INT(iEnd)
      ,CHFp_CONST_FRA1(areas)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MULTIPLYFLUXESBYAREA");
 FORTRAN_NAME( MULTIPLYFLUXESBYAREA ,multiplyfluxesbyarea )(
      CHFt_FRA(F)
      ,CHFt_CONST_INT(iBgn)
      ,CHFt_CONST_INT(iEnd)
      ,CHFt_CONST_FRA1(areas)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMULTIPLYFLUXESBYAREA 

}

#endif
