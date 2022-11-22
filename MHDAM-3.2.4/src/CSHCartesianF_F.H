#ifndef _CSHCARTESIANF_F_H_
#define _CSHCARTESIANF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDTRANSFLUXES_CYL_UPDSTATE 
#define GUARDTRANSFLUXES_CYL_UPDSTATE 
// Prototype for Fortran procedure TRANSFLUXES_CYL_UPDSTATE ...
//
void FORTRAN_NAME( TRANSFLUXES_CYL_UPDSTATE ,transfluxes_cyl_updstate )(
      CHFp_FRA(F)
      ,CHFp_CONST_REAL(dy)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_TRANSFLUXES_CYL_UPDSTATE FORTRAN_NAME( inlineTRANSFLUXES_CYL_UPDSTATE, inlineTRANSFLUXES_CYL_UPDSTATE)
#define FORTNT_TRANSFLUXES_CYL_UPDSTATE FORTRAN_NAME( TRANSFLUXES_CYL_UPDSTATE, transfluxes_cyl_updstate)

inline void FORTRAN_NAME(inlineTRANSFLUXES_CYL_UPDSTATE, inlineTRANSFLUXES_CYL_UPDSTATE)(
      CHFp_FRA(F)
      ,CHFp_CONST_REAL(dy)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_TRANSFLUXES_CYL_UPDSTATE");
 FORTRAN_NAME( TRANSFLUXES_CYL_UPDSTATE ,transfluxes_cyl_updstate )(
      CHFt_FRA(F)
      ,CHFt_CONST_REAL(dy)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTRANSFLUXES_CYL_UPDSTATE 

#ifndef GUARDSCALINGFACTOR_CYL 
#define GUARDSCALINGFACTOR_CYL 
// Prototype for Fortran procedure SCALINGFACTOR_CYL ...
//
void FORTRAN_NAME( SCALINGFACTOR_CYL ,scalingfactor_cyl )(
      CHFp_FRA1(scalef)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SCALINGFACTOR_CYL FORTRAN_NAME( inlineSCALINGFACTOR_CYL, inlineSCALINGFACTOR_CYL)
#define FORTNT_SCALINGFACTOR_CYL FORTRAN_NAME( SCALINGFACTOR_CYL, scalingfactor_cyl)

inline void FORTRAN_NAME(inlineSCALINGFACTOR_CYL, inlineSCALINGFACTOR_CYL)(
      CHFp_FRA1(scalef)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SCALINGFACTOR_CYL");
 FORTRAN_NAME( SCALINGFACTOR_CYL ,scalingfactor_cyl )(
      CHFt_FRA1(scalef)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSCALINGFACTOR_CYL 

}

#endif
