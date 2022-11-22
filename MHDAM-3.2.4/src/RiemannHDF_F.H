#ifndef _RIEMANNHDF_F_H_
#define _RIEMANNHDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDRIEMANNHD_STATES 
#define GUARDRIEMANNHD_STATES 
// Prototype for Fortran procedure RIEMANNHD_STATES ...
//
void FORTRAN_NAME( RIEMANNHD_STATES ,riemannhd_states )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNHD_STATES FORTRAN_NAME( inlineRIEMANNHD_STATES, inlineRIEMANNHD_STATES)
#define FORTNT_RIEMANNHD_STATES FORTRAN_NAME( RIEMANNHD_STATES, riemannhd_states)

inline void FORTRAN_NAME(inlineRIEMANNHD_STATES, inlineRIEMANNHD_STATES)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNHD_STATES");
 FORTRAN_NAME( RIEMANNHD_STATES ,riemannhd_states )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNHD_STATES 

#ifndef GUARDRIEMANNHD_FLUXES 
#define GUARDRIEMANNHD_FLUXES 
// Prototype for Fortran procedure RIEMANNHD_FLUXES ...
//
void FORTRAN_NAME( RIEMANNHD_FLUXES ,riemannhd_fluxes )(
      CHFp_FRA(Fgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNHD_FLUXES FORTRAN_NAME( inlineRIEMANNHD_FLUXES, inlineRIEMANNHD_FLUXES)
#define FORTNT_RIEMANNHD_FLUXES FORTRAN_NAME( RIEMANNHD_FLUXES, riemannhd_fluxes)

inline void FORTRAN_NAME(inlineRIEMANNHD_FLUXES, inlineRIEMANNHD_FLUXES)(
      CHFp_FRA(Fgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNHD_FLUXES");
 FORTRAN_NAME( RIEMANNHD_FLUXES ,riemannhd_fluxes )(
      CHFt_FRA(Fgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNHD_FLUXES 

}

#endif
