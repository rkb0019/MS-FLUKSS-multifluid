#ifndef _VFROEF_F_H_
#define _VFROEF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDVFROE_STATES 
#define GUARDVFROE_STATES 
// Prototype for Fortran procedure VFROE_STATES ...
//
void FORTRAN_NAME( VFROE_STATES ,vfroe_states )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_VFROE_STATES FORTRAN_NAME( inlineVFROE_STATES, inlineVFROE_STATES)
#define FORTNT_VFROE_STATES FORTRAN_NAME( VFROE_STATES, vfroe_states)

inline void FORTRAN_NAME(inlineVFROE_STATES, inlineVFROE_STATES)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VFROE_STATES");
 FORTRAN_NAME( VFROE_STATES ,vfroe_states )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVFROE_STATES 

#ifndef GUARDVFROE_FLUXES 
#define GUARDVFROE_FLUXES 
// Prototype for Fortran procedure VFROE_FLUXES ...
//
void FORTRAN_NAME( VFROE_FLUXES ,vfroe_fluxes )(
      CHFp_FRA(Fgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_VFROE_FLUXES FORTRAN_NAME( inlineVFROE_FLUXES, inlineVFROE_FLUXES)
#define FORTNT_VFROE_FLUXES FORTRAN_NAME( VFROE_FLUXES, vfroe_fluxes)

inline void FORTRAN_NAME(inlineVFROE_FLUXES, inlineVFROE_FLUXES)(
      CHFp_FRA(Fgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VFROE_FLUXES");
 FORTRAN_NAME( VFROE_FLUXES ,vfroe_fluxes )(
      CHFt_FRA(Fgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVFROE_FLUXES 

}

#endif
