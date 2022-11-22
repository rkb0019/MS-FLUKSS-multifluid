#ifndef _HELIOSPHERICF_F_H_
#define _HELIOSPHERICF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETHELIOS 
#define GUARDSETHELIOS 
// Prototype for Fortran procedure SETHELIOS ...
//
void FORTRAN_NAME( SETHELIOS ,sethelios )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aulismN)
      ,CHFp_CONST_REAL(aulismV)
      ,CHFp_CONST_REAL(aulismUX)
      ,CHFp_CONST_REAL(aulismUY)
      ,CHFp_CONST_REAL(aulismUZ)
      ,CHFp_CONST_REAL(aulismT)
      ,CHFp_CONST_REAL(aulismB)
      ,CHFp_CONST_REAL(aulismBX)
      ,CHFp_CONST_REAL(aulismBY)
      ,CHFp_CONST_REAL(aulismBZ)
      ,CHFp_CONST_REAL(ausunXC)
      ,CHFp_CONST_REAL(ausunYC)
      ,CHFp_CONST_REAL(ausunZC)
      ,CHFp_CONST_REAL(ausunR)
      ,CHFp_CONST_REAL(ausunN)
      ,CHFp_CONST_REAL(ausunV)
      ,CHFp_CONST_REAL(ausunT)
      ,CHFp_CONST_REAL(ausunB)
      ,CHFp_CONST_REAL(ausunTILT)
      ,CHFp_CONST_REAL(ausunIntBCRadius)
      ,CHFp_CONST_REAL(auinitR)
      ,CHFp_CONST_REAL(aunetN)
      ,CHFp_CONST_REAL(audirBrN)
      ,CHFp_CONST_REAL(auTMLIM)
      ,CHFp_CONST_REAL(auRegSafeZTop)
      ,CHFp_CONST_REAL(auRegSafeZBot)
      ,CHFp_CONST_INT(ausunBmonopolar) );

#define FORT_SETHELIOS FORTRAN_NAME( inlineSETHELIOS, inlineSETHELIOS)
#define FORTNT_SETHELIOS FORTRAN_NAME( SETHELIOS, sethelios)

inline void FORTRAN_NAME(inlineSETHELIOS, inlineSETHELIOS)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aulismN)
      ,CHFp_CONST_REAL(aulismV)
      ,CHFp_CONST_REAL(aulismUX)
      ,CHFp_CONST_REAL(aulismUY)
      ,CHFp_CONST_REAL(aulismUZ)
      ,CHFp_CONST_REAL(aulismT)
      ,CHFp_CONST_REAL(aulismB)
      ,CHFp_CONST_REAL(aulismBX)
      ,CHFp_CONST_REAL(aulismBY)
      ,CHFp_CONST_REAL(aulismBZ)
      ,CHFp_CONST_REAL(ausunXC)
      ,CHFp_CONST_REAL(ausunYC)
      ,CHFp_CONST_REAL(ausunZC)
      ,CHFp_CONST_REAL(ausunR)
      ,CHFp_CONST_REAL(ausunN)
      ,CHFp_CONST_REAL(ausunV)
      ,CHFp_CONST_REAL(ausunT)
      ,CHFp_CONST_REAL(ausunB)
      ,CHFp_CONST_REAL(ausunTILT)
      ,CHFp_CONST_REAL(ausunIntBCRadius)
      ,CHFp_CONST_REAL(auinitR)
      ,CHFp_CONST_REAL(aunetN)
      ,CHFp_CONST_REAL(audirBrN)
      ,CHFp_CONST_REAL(auTMLIM)
      ,CHFp_CONST_REAL(auRegSafeZTop)
      ,CHFp_CONST_REAL(auRegSafeZBot)
      ,CHFp_CONST_INT(ausunBmonopolar) )
{
 CH_TIMELEAF("FORT_SETHELIOS");
 FORTRAN_NAME( SETHELIOS ,sethelios )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aulismN)
      ,CHFt_CONST_REAL(aulismV)
      ,CHFt_CONST_REAL(aulismUX)
      ,CHFt_CONST_REAL(aulismUY)
      ,CHFt_CONST_REAL(aulismUZ)
      ,CHFt_CONST_REAL(aulismT)
      ,CHFt_CONST_REAL(aulismB)
      ,CHFt_CONST_REAL(aulismBX)
      ,CHFt_CONST_REAL(aulismBY)
      ,CHFt_CONST_REAL(aulismBZ)
      ,CHFt_CONST_REAL(ausunXC)
      ,CHFt_CONST_REAL(ausunYC)
      ,CHFt_CONST_REAL(ausunZC)
      ,CHFt_CONST_REAL(ausunR)
      ,CHFt_CONST_REAL(ausunN)
      ,CHFt_CONST_REAL(ausunV)
      ,CHFt_CONST_REAL(ausunT)
      ,CHFt_CONST_REAL(ausunB)
      ,CHFt_CONST_REAL(ausunTILT)
      ,CHFt_CONST_REAL(ausunIntBCRadius)
      ,CHFt_CONST_REAL(auinitR)
      ,CHFt_CONST_REAL(aunetN)
      ,CHFt_CONST_REAL(audirBrN)
      ,CHFt_CONST_REAL(auTMLIM)
      ,CHFt_CONST_REAL(auRegSafeZTop)
      ,CHFt_CONST_REAL(auRegSafeZBot)
      ,CHFt_CONST_INT(ausunBmonopolar) );
}
#endif  // GUARDSETHELIOS 

#ifndef GUARDSETFASTSLOWSW 
#define GUARDSETFASTSLOWSW 
// Prototype for Fortran procedure setfastslowsw ...
//
void FORTRAN_NAME( SETFASTSLOWSW ,setfastslowsw )(
      CHFp_CONST_REAL(aur)
      ,CHFp_CONST_REAL(auVfast)
      ,CHFp_CONST_REAL(auVslow)
      ,CHFp_CONST_REAL(auTfast)
      ,CHFp_CONST_REAL(auTslow)
      ,CHFp_CONST_REAL(auNfast)
      ,CHFp_CONST_REAL(auNslow)
      ,CHFp_CONST_REAL(aufs_a)
      ,CHFp_CONST_REAL(aufs_b)
      ,CHFp_CONST_REAL(autilt_min)
      ,CHFp_CONST_REAL(autilt_max)
      ,CHFp_CONST_REAL(auslow_min)
      ,CHFp_CONST_REAL(auslow_max) );

#define FORT_SETFASTSLOWSW FORTRAN_NAME( inlineSETFASTSLOWSW, inlineSETFASTSLOWSW)
#define FORTNT_SETFASTSLOWSW FORTRAN_NAME( SETFASTSLOWSW, setfastslowsw)

inline void FORTRAN_NAME(inlineSETFASTSLOWSW, inlineSETFASTSLOWSW)(
      CHFp_CONST_REAL(aur)
      ,CHFp_CONST_REAL(auVfast)
      ,CHFp_CONST_REAL(auVslow)
      ,CHFp_CONST_REAL(auTfast)
      ,CHFp_CONST_REAL(auTslow)
      ,CHFp_CONST_REAL(auNfast)
      ,CHFp_CONST_REAL(auNslow)
      ,CHFp_CONST_REAL(aufs_a)
      ,CHFp_CONST_REAL(aufs_b)
      ,CHFp_CONST_REAL(autilt_min)
      ,CHFp_CONST_REAL(autilt_max)
      ,CHFp_CONST_REAL(auslow_min)
      ,CHFp_CONST_REAL(auslow_max) )
{
 CH_TIMELEAF("FORT_SETFASTSLOWSW");
 FORTRAN_NAME( SETFASTSLOWSW ,setfastslowsw )(
      CHFt_CONST_REAL(aur)
      ,CHFt_CONST_REAL(auVfast)
      ,CHFt_CONST_REAL(auVslow)
      ,CHFt_CONST_REAL(auTfast)
      ,CHFt_CONST_REAL(auTslow)
      ,CHFt_CONST_REAL(auNfast)
      ,CHFt_CONST_REAL(auNslow)
      ,CHFt_CONST_REAL(aufs_a)
      ,CHFt_CONST_REAL(aufs_b)
      ,CHFt_CONST_REAL(autilt_min)
      ,CHFt_CONST_REAL(autilt_max)
      ,CHFt_CONST_REAL(auslow_min)
      ,CHFt_CONST_REAL(auslow_max) );
}
#endif  // GUARDSETFASTSLOWSW 

#ifndef GUARDHELIOBC 
#define GUARDHELIOBC 
// Prototype for Fortran procedure HELIOBC ...
//
void FORTRAN_NAME( HELIOBC ,heliobc )(
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

#define FORT_HELIOBC FORTRAN_NAME( inlineHELIOBC, inlineHELIOBC)
#define FORTNT_HELIOBC FORTRAN_NAME( HELIOBC, heliobc)

inline void FORTRAN_NAME(inlineHELIOBC, inlineHELIOBC)(
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
 CH_TIMELEAF("FORT_HELIOBC");
 FORTRAN_NAME( HELIOBC ,heliobc )(
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
#endif  // GUARDHELIOBC 

#ifndef GUARDNEUTRALS_SUNBC3D 
#define GUARDNEUTRALS_SUNBC3D 
// Prototype for Fortran procedure NEUTRALS_SUNBC3D ...
//
void FORTRAN_NAME( NEUTRALS_SUNBC3D ,neutrals_sunbc3d )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UR0)
      ,CHFp_CONST_INTVECT(sunIJK)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(b) );

#define FORT_NEUTRALS_SUNBC3D FORTRAN_NAME( inlineNEUTRALS_SUNBC3D, inlineNEUTRALS_SUNBC3D)
#define FORTNT_NEUTRALS_SUNBC3D FORTRAN_NAME( NEUTRALS_SUNBC3D, neutrals_sunbc3d)

inline void FORTRAN_NAME(inlineNEUTRALS_SUNBC3D, inlineNEUTRALS_SUNBC3D)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(UR0)
      ,CHFp_CONST_INTVECT(sunIJK)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(b) )
{
 CH_TIMELEAF("FORT_NEUTRALS_SUNBC3D");
 FORTRAN_NAME( NEUTRALS_SUNBC3D ,neutrals_sunbc3d )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(UR0)
      ,CHFt_CONST_INTVECT(sunIJK)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(b) );
}
#endif  // GUARDNEUTRALS_SUNBC3D 

#ifndef GUARDHELIOGS 
#define GUARDHELIOGS 
// Prototype for Fortran procedure HELIOGS ...
//
void FORTRAN_NAME( HELIOGS ,heliogs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) );

#define FORT_HELIOGS FORTRAN_NAME( inlineHELIOGS, inlineHELIOGS)
#define FORTNT_HELIOGS FORTRAN_NAME( HELIOGS, heliogs)

inline void FORTRAN_NAME(inlineHELIOGS, inlineHELIOGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGS");
 FORTRAN_NAME( HELIOGS ,heliogs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGS 

#ifndef GUARDHELIOINIT 
#define GUARDHELIOINIT 
// Prototype for Fortran procedure HELIOINIT ...
//
void FORTRAN_NAME( HELIOINIT ,helioinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) );

#define FORT_HELIOINIT FORTRAN_NAME( inlineHELIOINIT, inlineHELIOINIT)
#define FORTNT_HELIOINIT FORTRAN_NAME( HELIOINIT, helioinit)

inline void FORTRAN_NAME(inlineHELIOINIT, inlineHELIOINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINIT");
 FORTRAN_NAME( HELIOINIT ,helioinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINIT 

#ifndef GUARDBFLUX 
#define GUARDBFLUX 
// Prototype for Fortran procedure BFLUX ...
//
void FORTRAN_NAME( BFLUX ,bflux )(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(Bn_in)
      ,CHFp_CONST_REAL(Bt1_in)
      ,CHFp_CONST_REAL(Bt2_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_CONST_REAL(Bn_out)
      ,CHFp_CONST_REAL(Bt1_out)
      ,CHFp_CONST_REAL(Bt2_out)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBNORM)
      ,CHFp_REAL(FLUXBTAN1)
      ,CHFp_REAL(FLUXBTAN2)
      ,CHFp_REAL(BnValue) );

#define FORT_BFLUX FORTRAN_NAME( inlineBFLUX, inlineBFLUX)
#define FORTNT_BFLUX FORTRAN_NAME( BFLUX, bflux)

inline void FORTRAN_NAME(inlineBFLUX, inlineBFLUX)(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(Bn_in)
      ,CHFp_CONST_REAL(Bt1_in)
      ,CHFp_CONST_REAL(Bt2_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_CONST_REAL(Bn_out)
      ,CHFp_CONST_REAL(Bt1_out)
      ,CHFp_CONST_REAL(Bt2_out)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBNORM)
      ,CHFp_REAL(FLUXBTAN1)
      ,CHFp_REAL(FLUXBTAN2)
      ,CHFp_REAL(BnValue) )
{
 CH_TIMELEAF("FORT_BFLUX");
 FORTRAN_NAME( BFLUX ,bflux )(
      CHFt_CONST_REAL(rho_in)
      ,CHFt_CONST_REAL(vn_in)
      ,CHFt_CONST_REAL(vt1_in)
      ,CHFt_CONST_REAL(vt2_in)
      ,CHFt_CONST_REAL(p_in)
      ,CHFt_CONST_REAL(Bn_in)
      ,CHFt_CONST_REAL(Bt1_in)
      ,CHFt_CONST_REAL(Bt2_in)
      ,CHFt_CONST_REAL(rho_out)
      ,CHFt_CONST_REAL(vn_out)
      ,CHFt_CONST_REAL(vt1_out)
      ,CHFt_CONST_REAL(vt2_out)
      ,CHFt_CONST_REAL(p_out)
      ,CHFt_CONST_REAL(Bn_out)
      ,CHFt_CONST_REAL(Bt1_out)
      ,CHFt_CONST_REAL(Bt2_out)
      ,CHFt_REAL(FLUXRHO)
      ,CHFt_REAL(FLUXRUNORM)
      ,CHFt_REAL(FLUXRUTAN1)
      ,CHFt_REAL(FLUXRUTAN2)
      ,CHFt_REAL(FLUXE)
      ,CHFt_REAL(FLUXBNORM)
      ,CHFt_REAL(FLUXBTAN1)
      ,CHFt_REAL(FLUXBTAN2)
      ,CHFt_REAL(BnValue) );
}
#endif  // GUARDBFLUX 

#ifndef GUARDBFLUXR 
#define GUARDBFLUXR 
// Prototype for Fortran procedure BFLUXR ...
//
void FORTRAN_NAME( BFLUXR ,bfluxr )(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(Bn_in)
      ,CHFp_CONST_REAL(Bt1_in)
      ,CHFp_CONST_REAL(Bt2_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_CONST_REAL(Bn_out)
      ,CHFp_CONST_REAL(Bt1_out)
      ,CHFp_CONST_REAL(Bt2_out)
      ,CHFp_CONST_REAL(DX)
      ,CHFp_CONST_REAL(DY)
      ,CHFp_CONST_REAL(DZ)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBNORM)
      ,CHFp_REAL(FLUXBTAN1)
      ,CHFp_REAL(FLUXBTAN2)
      ,CHFp_REAL(BnValue) );

#define FORT_BFLUXR FORTRAN_NAME( inlineBFLUXR, inlineBFLUXR)
#define FORTNT_BFLUXR FORTRAN_NAME( BFLUXR, bfluxr)

inline void FORTRAN_NAME(inlineBFLUXR, inlineBFLUXR)(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(Bn_in)
      ,CHFp_CONST_REAL(Bt1_in)
      ,CHFp_CONST_REAL(Bt2_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_CONST_REAL(Bn_out)
      ,CHFp_CONST_REAL(Bt1_out)
      ,CHFp_CONST_REAL(Bt2_out)
      ,CHFp_CONST_REAL(DX)
      ,CHFp_CONST_REAL(DY)
      ,CHFp_CONST_REAL(DZ)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_REAL(FLUXBNORM)
      ,CHFp_REAL(FLUXBTAN1)
      ,CHFp_REAL(FLUXBTAN2)
      ,CHFp_REAL(BnValue) )
{
 CH_TIMELEAF("FORT_BFLUXR");
 FORTRAN_NAME( BFLUXR ,bfluxr )(
      CHFt_CONST_REAL(rho_in)
      ,CHFt_CONST_REAL(vn_in)
      ,CHFt_CONST_REAL(vt1_in)
      ,CHFt_CONST_REAL(vt2_in)
      ,CHFt_CONST_REAL(p_in)
      ,CHFt_CONST_REAL(Bn_in)
      ,CHFt_CONST_REAL(Bt1_in)
      ,CHFt_CONST_REAL(Bt2_in)
      ,CHFt_CONST_REAL(rho_out)
      ,CHFt_CONST_REAL(vn_out)
      ,CHFt_CONST_REAL(vt1_out)
      ,CHFt_CONST_REAL(vt2_out)
      ,CHFt_CONST_REAL(p_out)
      ,CHFt_CONST_REAL(Bn_out)
      ,CHFt_CONST_REAL(Bt1_out)
      ,CHFt_CONST_REAL(Bt2_out)
      ,CHFt_CONST_REAL(DX)
      ,CHFt_CONST_REAL(DY)
      ,CHFt_CONST_REAL(DZ)
      ,CHFt_REAL(FLUXRHO)
      ,CHFt_REAL(FLUXRUNORM)
      ,CHFt_REAL(FLUXRUTAN1)
      ,CHFt_REAL(FLUXRUTAN2)
      ,CHFt_REAL(FLUXE)
      ,CHFt_REAL(FLUXBNORM)
      ,CHFt_REAL(FLUXBTAN1)
      ,CHFt_REAL(FLUXBTAN2)
      ,CHFt_REAL(BnValue) );
}
#endif  // GUARDBFLUXR 

#ifndef GUARDHELIOINITSPHERICAL 
#define GUARDHELIOINITSPHERICAL 
// Prototype for Fortran procedure HELIOINITSPHERICAL ...
//
void FORTRAN_NAME( HELIOINITSPHERICAL ,helioinitspherical )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOINITSPHERICAL FORTRAN_NAME( inlineHELIOINITSPHERICAL, inlineHELIOINITSPHERICAL)
#define FORTNT_HELIOINITSPHERICAL FORTRAN_NAME( HELIOINITSPHERICAL, helioinitspherical)

inline void FORTRAN_NAME(inlineHELIOINITSPHERICAL, inlineHELIOINITSPHERICAL)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOINITSPHERICAL");
 FORTRAN_NAME( HELIOINITSPHERICAL ,helioinitspherical )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOINITSPHERICAL 

#ifndef GUARDHELIOGSSPHERICAL 
#define GUARDHELIOGSSPHERICAL 
// Prototype for Fortran procedure HELIOGSSPHERICAL ...
//
void FORTRAN_NAME( HELIOGSSPHERICAL ,heliogsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSSPHERICAL FORTRAN_NAME( inlineHELIOGSSPHERICAL, inlineHELIOGSSPHERICAL)
#define FORTNT_HELIOGSSPHERICAL FORTRAN_NAME( HELIOGSSPHERICAL, heliogsspherical)

inline void FORTRAN_NAME(inlineHELIOGSSPHERICAL, inlineHELIOGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSSPHERICAL");
 FORTRAN_NAME( HELIOGSSPHERICAL ,heliogsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(jsize)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSSPHERICAL 

#ifndef GUARDHELIOGSNEUTRALS 
#define GUARDHELIOGSNEUTRALS 
// Prototype for Fortran procedure HELIOGSNEUTRALS ...
//
void FORTRAN_NAME( HELIOGSNEUTRALS ,heliogsneutrals )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) );

#define FORT_HELIOGSNEUTRALS FORTRAN_NAME( inlineHELIOGSNEUTRALS, inlineHELIOGSNEUTRALS)
#define FORTNT_HELIOGSNEUTRALS FORTRAN_NAME( HELIOGSNEUTRALS, heliogsneutrals)

inline void FORTRAN_NAME(inlineHELIOGSNEUTRALS, inlineHELIOGSNEUTRALS)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(t)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGSNEUTRALS");
 FORTRAN_NAME( HELIOGSNEUTRALS ,heliogsneutrals )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(t)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGSNEUTRALS 

#ifndef GUARDHELIOBCSPHERICAL 
#define GUARDHELIOBCSPHERICAL 
// Prototype for Fortran procedure HELIOBCSPHERICAL ...
//
void FORTRAN_NAME( HELIOBCSPHERICAL ,heliobcspherical )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_HELIOBCSPHERICAL FORTRAN_NAME( inlineHELIOBCSPHERICAL, inlineHELIOBCSPHERICAL)
#define FORTNT_HELIOBCSPHERICAL FORTRAN_NAME( HELIOBCSPHERICAL, heliobcspherical)

inline void FORTRAN_NAME(inlineHELIOBCSPHERICAL, inlineHELIOBCSPHERICAL)(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOBCSPHERICAL");
 FORTRAN_NAME( HELIOBCSPHERICAL ,heliobcspherical )(
      CHFt_FRA(F)
      ,CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOBCSPHERICAL 

#ifndef GUARDNEUTRALS_SUNBCSPHERICAL 
#define GUARDNEUTRALS_SUNBCSPHERICAL 
// Prototype for Fortran procedure NEUTRALS_SUNBCSPHERICAL ...
//
void FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL ,neutrals_sunbcspherical )(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(a_sign) );

#define FORT_NEUTRALS_SUNBCSPHERICAL FORTRAN_NAME( inlineNEUTRALS_SUNBCSPHERICAL, inlineNEUTRALS_SUNBCSPHERICAL)
#define FORTNT_NEUTRALS_SUNBCSPHERICAL FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL, neutrals_sunbcspherical)

inline void FORTRAN_NAME(inlineNEUTRALS_SUNBCSPHERICAL, inlineNEUTRALS_SUNBCSPHERICAL)(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(a_sign) )
{
 CH_TIMELEAF("FORT_NEUTRALS_SUNBCSPHERICAL");
 FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL ,neutrals_sunbcspherical )(
      CHFt_CONST_FRA(a_from)
      ,CHFt_FRA(a_to)
      ,CHFt_CONST_INT(jsize)
      ,CHFt_CONST_INT(a_sign) );
}
#endif  // GUARDNEUTRALS_SUNBCSPHERICAL 

#ifndef GUARDNEUTRALS_SUNBCSPHERICAL2 
#define GUARDNEUTRALS_SUNBCSPHERICAL2 
// Prototype for Fortran procedure NEUTRALS_SUNBCSPHERICAL2 ...
//
void FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL2 ,neutrals_sunbcspherical2 )(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(a_sign) );

#define FORT_NEUTRALS_SUNBCSPHERICAL2 FORTRAN_NAME( inlineNEUTRALS_SUNBCSPHERICAL2, inlineNEUTRALS_SUNBCSPHERICAL2)
#define FORTNT_NEUTRALS_SUNBCSPHERICAL2 FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL2, neutrals_sunbcspherical2)

inline void FORTRAN_NAME(inlineNEUTRALS_SUNBCSPHERICAL2, inlineNEUTRALS_SUNBCSPHERICAL2)(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(jsize)
      ,CHFp_CONST_INT(a_sign) )
{
 CH_TIMELEAF("FORT_NEUTRALS_SUNBCSPHERICAL2");
 FORTRAN_NAME( NEUTRALS_SUNBCSPHERICAL2 ,neutrals_sunbcspherical2 )(
      CHFt_CONST_FRA(a_from)
      ,CHFt_FRA(a_to)
      ,CHFt_CONST_INT(jsize)
      ,CHFt_CONST_INT(a_sign) );
}
#endif  // GUARDNEUTRALS_SUNBCSPHERICAL2 

#ifndef GUARDTRANSFORMTOROTATINGFRAME 
#define GUARDTRANSFORMTOROTATINGFRAME 
// Prototype for Fortran procedure TRANSFORMTOROTATINGFRAME ...
//
void FORTRAN_NAME( TRANSFORMTOROTATINGFRAME ,transformtorotatingframe )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(WI)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_TRANSFORMTOROTATINGFRAME FORTRAN_NAME( inlineTRANSFORMTOROTATINGFRAME, inlineTRANSFORMTOROTATINGFRAME)
#define FORTNT_TRANSFORMTOROTATINGFRAME FORTRAN_NAME( TRANSFORMTOROTATINGFRAME, transformtorotatingframe)

inline void FORTRAN_NAME(inlineTRANSFORMTOROTATINGFRAME, inlineTRANSFORMTOROTATINGFRAME)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(WI)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_TRANSFORMTOROTATINGFRAME");
 FORTRAN_NAME( TRANSFORMTOROTATINGFRAME ,transformtorotatingframe )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(WI)
      ,CHFt_CONST_REAL(time0)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDTRANSFORMTOROTATINGFRAME 

#ifndef GUARDSUNGRAVITYSPH 
#define GUARDSUNGRAVITYSPH 
// Prototype for Fortran procedure SUNGRAVITYSPH ...
//
void FORTRAN_NAME( SUNGRAVITYSPH ,sungravitysph )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SUNGRAVITYSPH FORTRAN_NAME( inlineSUNGRAVITYSPH, inlineSUNGRAVITYSPH)
#define FORTNT_SUNGRAVITYSPH FORTRAN_NAME( SUNGRAVITYSPH, sungravitysph)

inline void FORTRAN_NAME(inlineSUNGRAVITYSPH, inlineSUNGRAVITYSPH)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUNGRAVITYSPH");
 FORTRAN_NAME( SUNGRAVITYSPH ,sungravitysph )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUNGRAVITYSPH 

#ifndef GUARDSUNHEATINGSPH 
#define GUARDSUNHEATINGSPH 
// Prototype for Fortran procedure SUNHEATINGSPH ...
//
void FORTRAN_NAME( SUNHEATINGSPH ,sunheatingsph )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SUNHEATINGSPH FORTRAN_NAME( inlineSUNHEATINGSPH, inlineSUNHEATINGSPH)
#define FORTNT_SUNHEATINGSPH FORTRAN_NAME( SUNHEATINGSPH, sunheatingsph)

inline void FORTRAN_NAME(inlineSUNHEATINGSPH, inlineSUNHEATINGSPH)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUNHEATINGSPH");
 FORTRAN_NAME( SUNHEATINGSPH ,sunheatingsph )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUNHEATINGSPH 

}

#endif
