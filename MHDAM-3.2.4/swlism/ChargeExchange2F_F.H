#ifndef _CHARGEEXCHANGE2F_F_H_
#define _CHARGEEXCHANGE2F_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCHARGEEX_PARS 
#define GUARDSETCHARGEEX_PARS 
// Prototype for Fortran procedure SETCHARGEEX_PARS ...
//
void FORTRAN_NAME( SETCHARGEEX_PARS ,setchargeex_pars )(
      CHFp_CONST_REAL(auVel)
      ,CHFp_CONST_REAL(auLen)
      ,CHFp_CONST_REAL(auNDen) );

#define FORT_SETCHARGEEX_PARS FORTRAN_NAME( inlineSETCHARGEEX_PARS, inlineSETCHARGEEX_PARS)
#define FORTNT_SETCHARGEEX_PARS FORTRAN_NAME( SETCHARGEEX_PARS, setchargeex_pars)

inline void FORTRAN_NAME(inlineSETCHARGEEX_PARS, inlineSETCHARGEEX_PARS)(
      CHFp_CONST_REAL(auVel)
      ,CHFp_CONST_REAL(auLen)
      ,CHFp_CONST_REAL(auNDen) )
{
 CH_TIMELEAF("FORT_SETCHARGEEX_PARS");
 FORTRAN_NAME( SETCHARGEEX_PARS ,setchargeex_pars )(
      CHFt_CONST_REAL(auVel)
      ,CHFt_CONST_REAL(auLen)
      ,CHFt_CONST_REAL(auNDen) );
}
#endif  // GUARDSETCHARGEEX_PARS 

#ifndef GUARDPHOTOIONIZE_SPH 
#define GUARDPHOTOIONIZE_SPH 
// Prototype for Fortran procedure PHOTOIONIZE_SPH ...
//
void FORTRAN_NAME( PHOTOIONIZE_SPH ,photoionize_sph )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_PHOTOIONIZE_SPH FORTRAN_NAME( inlinePHOTOIONIZE_SPH, inlinePHOTOIONIZE_SPH)
#define FORTNT_PHOTOIONIZE_SPH FORTRAN_NAME( PHOTOIONIZE_SPH, photoionize_sph)

inline void FORTRAN_NAME(inlinePHOTOIONIZE_SPH, inlinePHOTOIONIZE_SPH)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_PHOTOIONIZE_SPH");
 FORTRAN_NAME( PHOTOIONIZE_SPH ,photoionize_sph )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDPHOTOIONIZE_SPH 

#ifndef GUARDCHARGE_EXCHANGE_EXPSPH_CONSTH 
#define GUARDCHARGE_EXCHANGE_EXPSPH_CONSTH 
// Prototype for Fortran procedure CHARGE_EXCHANGE_EXPSPH_CONSTH ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_EXPSPH_CONSTH ,charge_exchange_expsph_consth )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_R1D(Wnet,iWnethi0)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_CHARGE_EXCHANGE_EXPSPH_CONSTH FORTRAN_NAME( inlineCHARGE_EXCHANGE_EXPSPH_CONSTH, inlineCHARGE_EXCHANGE_EXPSPH_CONSTH)
#define FORTNT_CHARGE_EXCHANGE_EXPSPH_CONSTH FORTRAN_NAME( CHARGE_EXCHANGE_EXPSPH_CONSTH, charge_exchange_expsph_consth)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_EXPSPH_CONSTH, inlineCHARGE_EXCHANGE_EXPSPH_CONSTH)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_R1D(Wnet,iWnethi0)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_EXPSPH_CONSTH");
 FORTRAN_NAME( CHARGE_EXCHANGE_EXPSPH_CONSTH ,charge_exchange_expsph_consth )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_R1D(Wnet,iWnethi0)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDCHARGE_EXCHANGE_EXPSPH_CONSTH 

#ifndef GUARDCHARGE_EXCHANGE_CONSTH 
#define GUARDCHARGE_EXCHANGE_CONSTH 
// Prototype for Fortran procedure CHARGE_EXCHANGE_CONSTH ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_CONSTH ,charge_exchange_consth )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_R1D(Wnet,iWnethi0)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_CHARGE_EXCHANGE_CONSTH FORTRAN_NAME( inlineCHARGE_EXCHANGE_CONSTH, inlineCHARGE_EXCHANGE_CONSTH)
#define FORTNT_CHARGE_EXCHANGE_CONSTH FORTRAN_NAME( CHARGE_EXCHANGE_CONSTH, charge_exchange_consth)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_CONSTH, inlineCHARGE_EXCHANGE_CONSTH)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_R1D(Wnet,iWnethi0)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_CONSTH");
 FORTRAN_NAME( CHARGE_EXCHANGE_CONSTH ,charge_exchange_consth )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_R1D(Wnet,iWnethi0)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDCHARGE_EXCHANGE_CONSTH 

#ifndef GUARDCHARGE_EXCHANGE_2F 
#define GUARDCHARGE_EXCHANGE_2F 
// Prototype for Fortran procedure CHARGE_EXCHANGE_2F ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_2F ,charge_exchange_2f )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_2F FORTRAN_NAME( inlineCHARGE_EXCHANGE_2F, inlineCHARGE_EXCHANGE_2F)
#define FORTNT_CHARGE_EXCHANGE_2F FORTRAN_NAME( CHARGE_EXCHANGE_2F, charge_exchange_2f)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_2F, inlineCHARGE_EXCHANGE_2F)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_2F");
 FORTRAN_NAME( CHARGE_EXCHANGE_2F ,charge_exchange_2f )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_2F 

#ifndef GUARDCHARGE_EXCHANGE_3F 
#define GUARDCHARGE_EXCHANGE_3F 
// Prototype for Fortran procedure CHARGE_EXCHANGE_3F ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_3F ,charge_exchange_3f )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_3F FORTRAN_NAME( inlineCHARGE_EXCHANGE_3F, inlineCHARGE_EXCHANGE_3F)
#define FORTNT_CHARGE_EXCHANGE_3F FORTRAN_NAME( CHARGE_EXCHANGE_3F, charge_exchange_3f)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_3F, inlineCHARGE_EXCHANGE_3F)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_3F");
 FORTRAN_NAME( CHARGE_EXCHANGE_3F ,charge_exchange_3f )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_3F 

#ifndef GUARDCHARGE_EXCHANGE_4F 
#define GUARDCHARGE_EXCHANGE_4F 
// Prototype for Fortran procedure CHARGE_EXCHANGE_4F ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_4F ,charge_exchange_4f )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_4F FORTRAN_NAME( inlineCHARGE_EXCHANGE_4F, inlineCHARGE_EXCHANGE_4F)
#define FORTNT_CHARGE_EXCHANGE_4F FORTRAN_NAME( CHARGE_EXCHANGE_4F, charge_exchange_4f)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_4F, inlineCHARGE_EXCHANGE_4F)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_4F");
 FORTRAN_NAME( CHARGE_EXCHANGE_4F ,charge_exchange_4f )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_4F 

#ifndef GUARDCHARGE_EXCHANGE_MF 
#define GUARDCHARGE_EXCHANGE_MF 
// Prototype for Fortran procedure CHARGE_EXCHANGE_MF ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_MF ,charge_exchange_mf )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_MF FORTRAN_NAME( inlineCHARGE_EXCHANGE_MF, inlineCHARGE_EXCHANGE_MF)
#define FORTNT_CHARGE_EXCHANGE_MF FORTRAN_NAME( CHARGE_EXCHANGE_MF, charge_exchange_mf)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_MF, inlineCHARGE_EXCHANGE_MF)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_MF");
 FORTRAN_NAME( CHARGE_EXCHANGE_MF ,charge_exchange_mf )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_MF 

#ifndef GUARDCHARGE_EXCHANGE_2FPI 
#define GUARDCHARGE_EXCHANGE_2FPI 
// Prototype for Fortran procedure CHARGE_EXCHANGE_2FPI ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_2FPI ,charge_exchange_2fpi )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_2FPI FORTRAN_NAME( inlineCHARGE_EXCHANGE_2FPI, inlineCHARGE_EXCHANGE_2FPI)
#define FORTNT_CHARGE_EXCHANGE_2FPI FORTRAN_NAME( CHARGE_EXCHANGE_2FPI, charge_exchange_2fpi)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_2FPI, inlineCHARGE_EXCHANGE_2FPI)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_2FPI");
 FORTRAN_NAME( CHARGE_EXCHANGE_2FPI ,charge_exchange_2fpi )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoH)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_2FPI 

}

#endif
