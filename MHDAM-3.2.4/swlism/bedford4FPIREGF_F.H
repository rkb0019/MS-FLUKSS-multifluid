#ifndef _BEDFORD4FPIREGF_F_H_
#define _BEDFORD4FPIREGF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCHARGE_EXCHANGE_4FPI 
#define GUARDCHARGE_EXCHANGE_4FPI 
// Prototype for Fortran procedure CHARGE_EXCHANGE_4FPI ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_4FPI ,charge_exchange_4fpi )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_4FPI FORTRAN_NAME( inlineCHARGE_EXCHANGE_4FPI, inlineCHARGE_EXCHANGE_4FPI)
#define FORTNT_CHARGE_EXCHANGE_4FPI FORTRAN_NAME( CHARGE_EXCHANGE_4FPI, charge_exchange_4fpi)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_4FPI, inlineCHARGE_EXCHANGE_4FPI)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_4FPI");
 FORTRAN_NAME( CHARGE_EXCHANGE_4FPI ,charge_exchange_4fpi )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_4FPI 

}

#endif
