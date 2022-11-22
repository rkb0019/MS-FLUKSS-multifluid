#ifndef _ROTORF_F_H_
#define _ROTORF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETROTOR 
#define GUARDSETROTOR 
// Prototype for Fortran procedure SETROTOR ...
//
void FORTRAN_NAME( SETROTOR ,setrotor )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aurotorU)
      ,CHFp_CONST_REAL(aurotorP)
      ,CHFp_CONST_REAL(aurotorR)
      ,CHFp_CONST_REAL(aurotorB)
      ,CHFp_CONST_REAL(aurotorXC)
      ,CHFp_CONST_REAL(aurotorYC)
      ,CHFp_CONST_REAL(aurotorR0)
      ,CHFp_CONST_REAL(aurotorR1) );

#define FORT_SETROTOR FORTRAN_NAME( inlineSETROTOR, inlineSETROTOR)
#define FORTNT_SETROTOR FORTRAN_NAME( SETROTOR, setrotor)

inline void FORTRAN_NAME(inlineSETROTOR, inlineSETROTOR)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aurotorU)
      ,CHFp_CONST_REAL(aurotorP)
      ,CHFp_CONST_REAL(aurotorR)
      ,CHFp_CONST_REAL(aurotorB)
      ,CHFp_CONST_REAL(aurotorXC)
      ,CHFp_CONST_REAL(aurotorYC)
      ,CHFp_CONST_REAL(aurotorR0)
      ,CHFp_CONST_REAL(aurotorR1) )
{
 CH_TIMELEAF("FORT_SETROTOR");
 FORTRAN_NAME( SETROTOR ,setrotor )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aurotorU)
      ,CHFt_CONST_REAL(aurotorP)
      ,CHFt_CONST_REAL(aurotorR)
      ,CHFt_CONST_REAL(aurotorB)
      ,CHFt_CONST_REAL(aurotorXC)
      ,CHFt_CONST_REAL(aurotorYC)
      ,CHFt_CONST_REAL(aurotorR0)
      ,CHFt_CONST_REAL(aurotorR1) );
}
#endif  // GUARDSETROTOR 

#ifndef GUARDROTORINIT 
#define GUARDROTORINIT 
// Prototype for Fortran procedure ROTORINIT ...
//
void FORTRAN_NAME( ROTORINIT ,rotorinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_ROTORINIT FORTRAN_NAME( inlineROTORINIT, inlineROTORINIT)
#define FORTNT_ROTORINIT FORTRAN_NAME( ROTORINIT, rotorinit)

inline void FORTRAN_NAME(inlineROTORINIT, inlineROTORINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ROTORINIT");
 FORTRAN_NAME( ROTORINIT ,rotorinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDROTORINIT 

}

#endif
