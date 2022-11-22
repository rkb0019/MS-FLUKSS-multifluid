#ifndef _SWLISMF_F_H_
#define _SWLISMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETSWLISM 
#define GUARDSETSWLISM 
// Prototype for Fortran procedure SETSWLISM ...
//
void FORTRAN_NAME( SETSWLISM ,setswlism )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aulismN)
      ,CHFp_CONST_REAL(aulismV)
      ,CHFp_CONST_REAL(aulismT)
      ,CHFp_CONST_REAL(aulismB)
      ,CHFp_CONST_REAL(ausunXC)
      ,CHFp_CONST_REAL(ausunYC)
      ,CHFp_CONST_REAL(ausunZC)
      ,CHFp_CONST_REAL(ausunR)
      ,CHFp_CONST_REAL(ausunN)
      ,CHFp_CONST_REAL(ausunV)
      ,CHFp_CONST_REAL(ausunT)
      ,CHFp_CONST_REAL(auinitR)
      ,CHFp_CONST_REAL(aunetN)
      ,CHFp_CONST_REAL(auTMLIM)
      ,CHFp_CONST_REAL(auRegSafeZ) );

#define FORT_SETSWLISM FORTRAN_NAME( inlineSETSWLISM, inlineSETSWLISM)
#define FORTNT_SETSWLISM FORTRAN_NAME( SETSWLISM, setswlism)

inline void FORTRAN_NAME(inlineSETSWLISM, inlineSETSWLISM)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aulismN)
      ,CHFp_CONST_REAL(aulismV)
      ,CHFp_CONST_REAL(aulismT)
      ,CHFp_CONST_REAL(aulismB)
      ,CHFp_CONST_REAL(ausunXC)
      ,CHFp_CONST_REAL(ausunYC)
      ,CHFp_CONST_REAL(ausunZC)
      ,CHFp_CONST_REAL(ausunR)
      ,CHFp_CONST_REAL(ausunN)
      ,CHFp_CONST_REAL(ausunV)
      ,CHFp_CONST_REAL(ausunT)
      ,CHFp_CONST_REAL(auinitR)
      ,CHFp_CONST_REAL(aunetN)
      ,CHFp_CONST_REAL(auTMLIM)
      ,CHFp_CONST_REAL(auRegSafeZ) )
{
 CH_TIMELEAF("FORT_SETSWLISM");
 FORTRAN_NAME( SETSWLISM ,setswlism )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aulismN)
      ,CHFt_CONST_REAL(aulismV)
      ,CHFt_CONST_REAL(aulismT)
      ,CHFt_CONST_REAL(aulismB)
      ,CHFt_CONST_REAL(ausunXC)
      ,CHFt_CONST_REAL(ausunYC)
      ,CHFt_CONST_REAL(ausunZC)
      ,CHFt_CONST_REAL(ausunR)
      ,CHFt_CONST_REAL(ausunN)
      ,CHFt_CONST_REAL(ausunV)
      ,CHFt_CONST_REAL(ausunT)
      ,CHFt_CONST_REAL(auinitR)
      ,CHFt_CONST_REAL(aunetN)
      ,CHFt_CONST_REAL(auTMLIM)
      ,CHFt_CONST_REAL(auRegSafeZ) );
}
#endif  // GUARDSETSWLISM 

#ifndef GUARDSWLISMINIT 
#define GUARDSWLISMINIT 
// Prototype for Fortran procedure SWLISMINIT ...
//
void FORTRAN_NAME( SWLISMINIT ,swlisminit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_SWLISMINIT FORTRAN_NAME( inlineSWLISMINIT, inlineSWLISMINIT)
#define FORTNT_SWLISMINIT FORTRAN_NAME( SWLISMINIT, swlisminit)

inline void FORTRAN_NAME(inlineSWLISMINIT, inlineSWLISMINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWLISMINIT");
 FORTRAN_NAME( SWLISMINIT ,swlisminit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWLISMINIT 

#ifndef GUARDREGTRACER_REINIT 
#define GUARDREGTRACER_REINIT 
// Prototype for Fortran procedure REGTRACER_REINIT ...
//
void FORTRAN_NAME( REGTRACER_REINIT ,regtracer_reinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) );

#define FORT_REGTRACER_REINIT FORTRAN_NAME( inlineREGTRACER_REINIT, inlineREGTRACER_REINIT)
#define FORTNT_REGTRACER_REINIT FORTRAN_NAME( REGTRACER_REINIT, regtracer_reinit)

inline void FORTRAN_NAME(inlineREGTRACER_REINIT, inlineREGTRACER_REINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REGTRACER_REINIT");
 FORTRAN_NAME( REGTRACER_REINIT ,regtracer_reinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREGTRACER_REINIT 

#ifndef GUARDSWLISMREINIT 
#define GUARDSWLISMREINIT 
// Prototype for Fortran procedure SWLISMREINIT ...
//
void FORTRAN_NAME( SWLISMREINIT ,swlismreinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_SWLISMREINIT FORTRAN_NAME( inlineSWLISMREINIT, inlineSWLISMREINIT)
#define FORTNT_SWLISMREINIT FORTRAN_NAME( SWLISMREINIT, swlismreinit)

inline void FORTRAN_NAME(inlineSWLISMREINIT, inlineSWLISMREINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWLISMREINIT");
 FORTRAN_NAME( SWLISMREINIT ,swlismreinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWLISMREINIT 

#ifndef GUARDSWLISMDEDNERGS 
#define GUARDSWLISMDEDNERGS 
// Prototype for Fortran procedure SWLISMDEDNERGS ...
//
void FORTRAN_NAME( SWLISMDEDNERGS ,swlismdednergs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_INT(jtop)
      ,CHFp_BOX(box) );

#define FORT_SWLISMDEDNERGS FORTRAN_NAME( inlineSWLISMDEDNERGS, inlineSWLISMDEDNERGS)
#define FORTNT_SWLISMDEDNERGS FORTRAN_NAME( SWLISMDEDNERGS, swlismdednergs)

inline void FORTRAN_NAME(inlineSWLISMDEDNERGS, inlineSWLISMDEDNERGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_INT(jtop)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWLISMDEDNERGS");
 FORTRAN_NAME( SWLISMDEDNERGS ,swlismdednergs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_INT(jtop)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWLISMDEDNERGS 

#ifndef GUARDLISMINIT 
#define GUARDLISMINIT 
// Prototype for Fortran procedure LISMINIT ...
//
void FORTRAN_NAME( LISMINIT ,lisminit )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_LISMINIT FORTRAN_NAME( inlineLISMINIT, inlineLISMINIT)
#define FORTNT_LISMINIT FORTRAN_NAME( LISMINIT, lisminit)

inline void FORTRAN_NAME(inlineLISMINIT, inlineLISMINIT)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LISMINIT");
 FORTRAN_NAME( LISMINIT ,lisminit )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLISMINIT 

#ifndef GUARDDEFINE_REGIONS_2F 
#define GUARDDEFINE_REGIONS_2F 
// Prototype for Fortran procedure DEFINE_REGIONS_2F ...
//
void FORTRAN_NAME( DEFINE_REGIONS_2F ,define_regions_2f )(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_DEFINE_REGIONS_2F FORTRAN_NAME( inlineDEFINE_REGIONS_2F, inlineDEFINE_REGIONS_2F)
#define FORTNT_DEFINE_REGIONS_2F FORTRAN_NAME( DEFINE_REGIONS_2F, define_regions_2f)

inline void FORTRAN_NAME(inlineDEFINE_REGIONS_2F, inlineDEFINE_REGIONS_2F)(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEFINE_REGIONS_2F");
 FORTRAN_NAME( DEFINE_REGIONS_2F ,define_regions_2f )(
      CHFt_CONST_FRA(W)
      ,CHFt_FIA1(REG)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEFINE_REGIONS_2F 

#ifndef GUARDDEFINE_REGIONS_TRACER 
#define GUARDDEFINE_REGIONS_TRACER 
// Prototype for Fortran procedure DEFINE_REGIONS_TRACER ...
//
void FORTRAN_NAME( DEFINE_REGIONS_TRACER ,define_regions_tracer )(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_DEFINE_REGIONS_TRACER FORTRAN_NAME( inlineDEFINE_REGIONS_TRACER, inlineDEFINE_REGIONS_TRACER)
#define FORTNT_DEFINE_REGIONS_TRACER FORTRAN_NAME( DEFINE_REGIONS_TRACER, define_regions_tracer)

inline void FORTRAN_NAME(inlineDEFINE_REGIONS_TRACER, inlineDEFINE_REGIONS_TRACER)(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEFINE_REGIONS_TRACER");
 FORTRAN_NAME( DEFINE_REGIONS_TRACER ,define_regions_tracer )(
      CHFt_CONST_FRA(W)
      ,CHFt_FIA1(REG)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEFINE_REGIONS_TRACER 

#ifndef GUARDDEFINE_REGIONS_TRACER_4F 
#define GUARDDEFINE_REGIONS_TRACER_4F 
// Prototype for Fortran procedure DEFINE_REGIONS_TRACER_4F ...
//
void FORTRAN_NAME( DEFINE_REGIONS_TRACER_4F ,define_regions_tracer_4f )(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_DEFINE_REGIONS_TRACER_4F FORTRAN_NAME( inlineDEFINE_REGIONS_TRACER_4F, inlineDEFINE_REGIONS_TRACER_4F)
#define FORTNT_DEFINE_REGIONS_TRACER_4F FORTRAN_NAME( DEFINE_REGIONS_TRACER_4F, define_regions_tracer_4f)

inline void FORTRAN_NAME(inlineDEFINE_REGIONS_TRACER_4F, inlineDEFINE_REGIONS_TRACER_4F)(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEFINE_REGIONS_TRACER_4F");
 FORTRAN_NAME( DEFINE_REGIONS_TRACER_4F ,define_regions_tracer_4f )(
      CHFt_CONST_FRA(W)
      ,CHFt_FIA1(REG)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEFINE_REGIONS_TRACER_4F 

#ifndef GUARDSWLISMBC 
#define GUARDSWLISMBC 
// Prototype for Fortran procedure SWLISMBC ...
//
void FORTRAN_NAME( SWLISMBC ,swlismbc )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_SWLISMBC FORTRAN_NAME( inlineSWLISMBC, inlineSWLISMBC)
#define FORTNT_SWLISMBC FORTRAN_NAME( SWLISMBC, swlismbc)

inline void FORTRAN_NAME(inlineSWLISMBC, inlineSWLISMBC)(
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
 CH_TIMELEAF("FORT_SWLISMBC");
 FORTRAN_NAME( SWLISMBC ,swlismbc )(
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
#endif  // GUARDSWLISMBC 

#ifndef GUARDDEDNERBC 
#define GUARDDEDNERBC 
// Prototype for Fortran procedure DEDNERBC ...
//
void FORTRAN_NAME( DEDNERBC ,dednerbc )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_DEDNERBC FORTRAN_NAME( inlineDEDNERBC, inlineDEDNERBC)
#define FORTNT_DEDNERBC FORTRAN_NAME( DEDNERBC, dednerbc)

inline void FORTRAN_NAME(inlineDEDNERBC, inlineDEDNERBC)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEDNERBC");
 FORTRAN_NAME( DEDNERBC ,dednerbc )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEDNERBC 

#ifndef GUARDBFLUXMHD 
#define GUARDBFLUXMHD 
// Prototype for Fortran procedure BFLUXMHD ...
//
void FORTRAN_NAME( BFLUXMHD ,bfluxmhd )(
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
      ,CHFp_REAL(BnValue)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_BFLUXMHD FORTRAN_NAME( inlineBFLUXMHD, inlineBFLUXMHD)
#define FORTNT_BFLUXMHD FORTRAN_NAME( BFLUXMHD, bfluxmhd)

inline void FORTRAN_NAME(inlineBFLUXMHD, inlineBFLUXMHD)(
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
      ,CHFp_REAL(BnValue)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_BFLUXMHD");
 FORTRAN_NAME( BFLUXMHD ,bfluxmhd )(
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
      ,CHFt_REAL(BnValue)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDBFLUXMHD 

#ifndef GUARDBFLUXGAS 
#define GUARDBFLUXGAS 
// Prototype for Fortran procedure BFLUXGAS ...
//
void FORTRAN_NAME( BFLUXGAS ,bfluxgas )(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign) );

#define FORT_BFLUXGAS FORTRAN_NAME( inlineBFLUXGAS, inlineBFLUXGAS)
#define FORTNT_BFLUXGAS FORTRAN_NAME( BFLUXGAS, bfluxgas)

inline void FORTRAN_NAME(inlineBFLUXGAS, inlineBFLUXGAS)(
      CHFp_CONST_REAL(rho_in)
      ,CHFp_CONST_REAL(vn_in)
      ,CHFp_CONST_REAL(vt1_in)
      ,CHFp_CONST_REAL(vt2_in)
      ,CHFp_CONST_REAL(p_in)
      ,CHFp_CONST_REAL(rho_out)
      ,CHFp_CONST_REAL(vn_out)
      ,CHFp_CONST_REAL(vt1_out)
      ,CHFp_CONST_REAL(vt2_out)
      ,CHFp_CONST_REAL(p_out)
      ,CHFp_REAL(FLUXRHO)
      ,CHFp_REAL(FLUXRUNORM)
      ,CHFp_REAL(FLUXRUTAN1)
      ,CHFp_REAL(FLUXRUTAN2)
      ,CHFp_REAL(FLUXE)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(lohisign) )
{
 CH_TIMELEAF("FORT_BFLUXGAS");
 FORTRAN_NAME( BFLUXGAS ,bfluxgas )(
      CHFt_CONST_REAL(rho_in)
      ,CHFt_CONST_REAL(vn_in)
      ,CHFt_CONST_REAL(vt1_in)
      ,CHFt_CONST_REAL(vt2_in)
      ,CHFt_CONST_REAL(p_in)
      ,CHFt_CONST_REAL(rho_out)
      ,CHFt_CONST_REAL(vn_out)
      ,CHFt_CONST_REAL(vt1_out)
      ,CHFt_CONST_REAL(vt2_out)
      ,CHFt_CONST_REAL(p_out)
      ,CHFt_REAL(FLUXRHO)
      ,CHFt_REAL(FLUXRUNORM)
      ,CHFt_REAL(FLUXRUTAN1)
      ,CHFt_REAL(FLUXRUTAN2)
      ,CHFt_REAL(FLUXE)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(lohisign) );
}
#endif  // GUARDBFLUXGAS 

#ifndef GUARDNEUTRALS_SUNBC 
#define GUARDNEUTRALS_SUNBC 
// Prototype for Fortran procedure NEUTRALS_SUNBC ...
//
void FORTRAN_NAME( NEUTRALS_SUNBC ,neutrals_sunbc )(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INTVECT(sunIJK)
      ,CHFp_CONST_REAL(dx) );

#define FORT_NEUTRALS_SUNBC FORTRAN_NAME( inlineNEUTRALS_SUNBC, inlineNEUTRALS_SUNBC)
#define FORTNT_NEUTRALS_SUNBC FORTRAN_NAME( NEUTRALS_SUNBC, neutrals_sunbc)

inline void FORTRAN_NAME(inlineNEUTRALS_SUNBC, inlineNEUTRALS_SUNBC)(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INTVECT(sunIJK)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_NEUTRALS_SUNBC");
 FORTRAN_NAME( NEUTRALS_SUNBC ,neutrals_sunbc )(
      CHFt_CONST_FRA(a_from)
      ,CHFt_FRA(a_to)
      ,CHFt_CONST_INTVECT(sunIJK)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDNEUTRALS_SUNBC 

#ifndef GUARDLOCKEDCELLS_REGION 
#define GUARDLOCKEDCELLS_REGION 
// Prototype for Fortran procedure LOCKEDCELLS_REGION ...
//
void FORTRAN_NAME( LOCKEDCELLS_REGION ,lockedcells_region )(
      CHFp_FIA1(flag)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) );

#define FORT_LOCKEDCELLS_REGION FORTRAN_NAME( inlineLOCKEDCELLS_REGION, inlineLOCKEDCELLS_REGION)
#define FORTNT_LOCKEDCELLS_REGION FORTRAN_NAME( LOCKEDCELLS_REGION, lockedcells_region)

inline void FORTRAN_NAME(inlineLOCKEDCELLS_REGION, inlineLOCKEDCELLS_REGION)(
      CHFp_FIA1(flag)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(loBox)
      ,CHFp_CONST_INT(hasLo)
      ,CHFp_BOX(hiBox)
      ,CHFp_CONST_INT(hasHi)
      ,CHFp_BOX(centerBox) )
{
 CH_TIMELEAF("FORT_LOCKEDCELLS_REGION");
 FORTRAN_NAME( LOCKEDCELLS_REGION ,lockedcells_region )(
      CHFt_FIA1(flag)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(loBox)
      ,CHFt_CONST_INT(hasLo)
      ,CHFt_BOX(hiBox)
      ,CHFt_CONST_INT(hasHi)
      ,CHFt_BOX(centerBox) );
}
#endif  // GUARDLOCKEDCELLS_REGION 

#ifndef GUARDLISMINITPOLAR 
#define GUARDLISMINITPOLAR 
// Prototype for Fortran procedure LISMINITPOLAR ...
//
void FORTRAN_NAME( LISMINITPOLAR ,lisminitpolar )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_LISMINITPOLAR FORTRAN_NAME( inlineLISMINITPOLAR, inlineLISMINITPOLAR)
#define FORTNT_LISMINITPOLAR FORTRAN_NAME( LISMINITPOLAR, lisminitpolar)

inline void FORTRAN_NAME(inlineLISMINITPOLAR, inlineLISMINITPOLAR)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iRegTr)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_LISMINITPOLAR");
 FORTRAN_NAME( LISMINITPOLAR ,lisminitpolar )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iRegTr)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDLISMINITPOLAR 

#ifndef GUARDSWINITPOLAR 
#define GUARDSWINITPOLAR 
// Prototype for Fortran procedure SWINITPOLAR ...
//
void FORTRAN_NAME( SWINITPOLAR ,swinitpolar )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SWINITPOLAR FORTRAN_NAME( inlineSWINITPOLAR, inlineSWINITPOLAR)
#define FORTNT_SWINITPOLAR FORTRAN_NAME( SWINITPOLAR, swinitpolar)

inline void FORTRAN_NAME(inlineSWINITPOLAR, inlineSWINITPOLAR)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWINITPOLAR");
 FORTRAN_NAME( SWINITPOLAR ,swinitpolar )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHP)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWINITPOLAR 

#ifndef GUARDPOLARGSAXIS 
#define GUARDPOLARGSAXIS 
// Prototype for Fortran procedure POLARGSAXIS ...
//
void FORTRAN_NAME( POLARGSAXIS ,polargsaxis )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(jbnd)
      ,CHFp_CONST_INT(nGS)
      ,CHFp_BOX(box) );

#define FORT_POLARGSAXIS FORTRAN_NAME( inlinePOLARGSAXIS, inlinePOLARGSAXIS)
#define FORTNT_POLARGSAXIS FORTRAN_NAME( POLARGSAXIS, polargsaxis)

inline void FORTRAN_NAME(inlinePOLARGSAXIS, inlinePOLARGSAXIS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(jbnd)
      ,CHFp_CONST_INT(nGS)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POLARGSAXIS");
 FORTRAN_NAME( POLARGSAXIS ,polargsaxis )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHP)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(jbnd)
      ,CHFt_CONST_INT(nGS)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOLARGSAXIS 

#ifndef GUARDSWLISMINITPOLAR 
#define GUARDSWLISMINITPOLAR 
// Prototype for Fortran procedure SWLISMINITPOLAR ...
//
void FORTRAN_NAME( SWLISMINITPOLAR ,swlisminitpolar )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SWLISMINITPOLAR FORTRAN_NAME( inlineSWLISMINITPOLAR, inlineSWLISMINITPOLAR)
#define FORTNT_SWLISMINITPOLAR FORTRAN_NAME( SWLISMINITPOLAR, swlisminitpolar)

inline void FORTRAN_NAME(inlineSWLISMINITPOLAR, inlineSWLISMINITPOLAR)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(iHP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SWLISMINITPOLAR");
 FORTRAN_NAME( SWLISMINITPOLAR ,swlisminitpolar )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(iHP)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSWLISMINITPOLAR 

#ifndef GUARDSWLISMBCPOLAR 
#define GUARDSWLISMBCPOLAR 
// Prototype for Fortran procedure SWLISMBCPOLAR ...
//
void FORTRAN_NAME( SWLISMBCPOLAR ,swlismbcpolar )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_SWLISMBCPOLAR FORTRAN_NAME( inlineSWLISMBCPOLAR, inlineSWLISMBCPOLAR)
#define FORTNT_SWLISMBCPOLAR FORTRAN_NAME( SWLISMBCPOLAR, swlismbcpolar)

inline void FORTRAN_NAME(inlineSWLISMBCPOLAR, inlineSWLISMBCPOLAR)(
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
 CH_TIMELEAF("FORT_SWLISMBCPOLAR");
 FORTRAN_NAME( SWLISMBCPOLAR ,swlismbcpolar )(
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
#endif  // GUARDSWLISMBCPOLAR 

#ifndef GUARDNEUTRALS_SUNBCPOLAR 
#define GUARDNEUTRALS_SUNBCPOLAR 
// Prototype for Fortran procedure NEUTRALS_SUNBCPOLAR ...
//
void FORTRAN_NAME( NEUTRALS_SUNBCPOLAR ,neutrals_sunbcpolar )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_FRA(UR0)
      ,CHFp_CONST_INT(jmax)
      ,CHFp_BOX(b) );

#define FORT_NEUTRALS_SUNBCPOLAR FORTRAN_NAME( inlineNEUTRALS_SUNBCPOLAR, inlineNEUTRALS_SUNBCPOLAR)
#define FORTNT_NEUTRALS_SUNBCPOLAR FORTRAN_NAME( NEUTRALS_SUNBCPOLAR, neutrals_sunbcpolar)

inline void FORTRAN_NAME(inlineNEUTRALS_SUNBCPOLAR, inlineNEUTRALS_SUNBCPOLAR)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_FRA(UR0)
      ,CHFp_CONST_INT(jmax)
      ,CHFp_BOX(b) )
{
 CH_TIMELEAF("FORT_NEUTRALS_SUNBCPOLAR");
 FORTRAN_NAME( NEUTRALS_SUNBCPOLAR ,neutrals_sunbcpolar )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_FRA(UR0)
      ,CHFt_CONST_INT(jmax)
      ,CHFt_BOX(b) );
}
#endif  // GUARDNEUTRALS_SUNBCPOLAR 

#ifndef GUARDSWLISMBCDEDNER 
#define GUARDSWLISMBCDEDNER 
// Prototype for Fortran procedure SWLISMBCDEDNER ...
//
void FORTRAN_NAME( SWLISMBCDEDNER ,swlismbcdedner )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_SWLISMBCDEDNER FORTRAN_NAME( inlineSWLISMBCDEDNER, inlineSWLISMBCDEDNER)
#define FORTNT_SWLISMBCDEDNER FORTRAN_NAME( SWLISMBCDEDNER, swlismbcdedner)

inline void FORTRAN_NAME(inlineSWLISMBCDEDNER, inlineSWLISMBCDEDNER)(
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
 CH_TIMELEAF("FORT_SWLISMBCDEDNER");
 FORTRAN_NAME( SWLISMBCDEDNER ,swlismbcdedner )(
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
#endif  // GUARDSWLISMBCDEDNER 

#ifndef GUARDBFLUXDEDNER 
#define GUARDBFLUXDEDNER 
// Prototype for Fortran procedure BFLUXDEDNER ...
//
void FORTRAN_NAME( BFLUXDEDNER ,bfluxdedner )(
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

#define FORT_BFLUXDEDNER FORTRAN_NAME( inlineBFLUXDEDNER, inlineBFLUXDEDNER)
#define FORTNT_BFLUXDEDNER FORTRAN_NAME( BFLUXDEDNER, bfluxdedner)

inline void FORTRAN_NAME(inlineBFLUXDEDNER, inlineBFLUXDEDNER)(
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
 CH_TIMELEAF("FORT_BFLUXDEDNER");
 FORTRAN_NAME( BFLUXDEDNER ,bfluxdedner )(
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
#endif  // GUARDBFLUXDEDNER 

}

#endif
