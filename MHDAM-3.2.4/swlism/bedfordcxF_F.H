#ifndef _BEDFORDCXF_F_H_
#define _BEDFORDCXF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCHARGE_EXCHANGE_2FPI_REG 
#define GUARDCHARGE_EXCHANGE_2FPI_REG 
// Prototype for Fortran procedure CHARGE_EXCHANGE_2FPI_REG ...
//
void FORTRAN_NAME( CHARGE_EXCHANGE_2FPI_REG ,charge_exchange_2fpi_reg )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) );

#define FORT_CHARGE_EXCHANGE_2FPI_REG FORTRAN_NAME( inlineCHARGE_EXCHANGE_2FPI_REG, inlineCHARGE_EXCHANGE_2FPI_REG)
#define FORTNT_CHARGE_EXCHANGE_2FPI_REG FORTRAN_NAME( CHARGE_EXCHANGE_2FPI_REG, charge_exchange_2fpi_reg)

inline void FORTRAN_NAME(inlineCHARGE_EXCHANGE_2FPI_REG, inlineCHARGE_EXCHANGE_2FPI_REG)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FIA1(REG)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRhoH)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARGE_EXCHANGE_2FPI_REG");
 FORTRAN_NAME( CHARGE_EXCHANGE_2FPI_REG ,charge_exchange_2fpi_reg )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FIA1(REG)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRhoH)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARGE_EXCHANGE_2FPI_REG 

}

#endif
