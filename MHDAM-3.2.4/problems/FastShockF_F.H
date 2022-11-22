#ifndef _FASTSHOCKF_F_H_
#define _FASTSHOCKF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETFASTSHOCK 
#define GUARDSETFASTSHOCK 
// Prototype for Fortran procedure SETFASTSHOCK ...
//
void FORTRAN_NAME( SETFASTSHOCK ,setfastshock )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auVelModuleL)
      ,CHFp_CONST_REAL(auVelAngleL)
      ,CHFp_CONST_REAL(auMagModuleL)
      ,CHFp_CONST_REAL(auMagAngleL)
      ,CHFp_CONST_REAL(auPressRatioL)
      ,CHFp_CONST_REAL(auFastMachL)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auAlfvenMag)
      ,CHFp_CONST_REAL(auAlfvenK) );

#define FORT_SETFASTSHOCK FORTRAN_NAME( inlineSETFASTSHOCK, inlineSETFASTSHOCK)
#define FORTNT_SETFASTSHOCK FORTRAN_NAME( SETFASTSHOCK, setfastshock)

inline void FORTRAN_NAME(inlineSETFASTSHOCK, inlineSETFASTSHOCK)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auVelModuleL)
      ,CHFp_CONST_REAL(auVelAngleL)
      ,CHFp_CONST_REAL(auMagModuleL)
      ,CHFp_CONST_REAL(auMagAngleL)
      ,CHFp_CONST_REAL(auPressRatioL)
      ,CHFp_CONST_REAL(auFastMachL)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_REAL(auAlfvenMag)
      ,CHFp_CONST_REAL(auAlfvenK) )
{
 CH_TIMELEAF("FORT_SETFASTSHOCK");
 FORTRAN_NAME( SETFASTSHOCK ,setfastshock )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auVelModuleL)
      ,CHFt_CONST_REAL(auVelAngleL)
      ,CHFt_CONST_REAL(auMagModuleL)
      ,CHFt_CONST_REAL(auMagAngleL)
      ,CHFt_CONST_REAL(auPressRatioL)
      ,CHFt_CONST_REAL(auFastMachL)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_REAL(auAlfvenMag)
      ,CHFt_CONST_REAL(auAlfvenK) );
}
#endif  // GUARDSETFASTSHOCK 

#ifndef GUARDFASTSHOCKINIT 
#define GUARDFASTSHOCKINIT 
// Prototype for Fortran procedure FASTSHOCKINIT ...
//
void FORTRAN_NAME( FASTSHOCKINIT ,fastshockinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_FASTSHOCKINIT FORTRAN_NAME( inlineFASTSHOCKINIT, inlineFASTSHOCKINIT)
#define FORTNT_FASTSHOCKINIT FORTRAN_NAME( FASTSHOCKINIT, fastshockinit)

inline void FORTRAN_NAME(inlineFASTSHOCKINIT, inlineFASTSHOCKINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FASTSHOCKINIT");
 FORTRAN_NAME( FASTSHOCKINIT ,fastshockinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFASTSHOCKINIT 

#ifndef GUARDFASTSHOCKBC 
#define GUARDFASTSHOCKBC 
// Prototype for Fortran procedure FASTSHOCKBC ...
//
void FORTRAN_NAME( FASTSHOCKBC ,fastshockbc )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_FASTSHOCKBC FORTRAN_NAME( inlineFASTSHOCKBC, inlineFASTSHOCKBC)
#define FORTNT_FASTSHOCKBC FORTRAN_NAME( FASTSHOCKBC, fastshockbc)

inline void FORTRAN_NAME(inlineFASTSHOCKBC, inlineFASTSHOCKBC)(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FASTSHOCKBC");
 FORTRAN_NAME( FASTSHOCKBC ,fastshockbc )(
      CHFt_FRA(F)
      ,CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFASTSHOCKBC 

#ifndef GUARDFASTSHOCKGS 
#define GUARDFASTSHOCKGS 
// Prototype for Fortran procedure FASTSHOCKGS ...
//
void FORTRAN_NAME( FASTSHOCKGS ,fastshockgs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_FASTSHOCKGS FORTRAN_NAME( inlineFASTSHOCKGS, inlineFASTSHOCKGS)
#define FORTNT_FASTSHOCKGS FORTRAN_NAME( FASTSHOCKGS, fastshockgs)

inline void FORTRAN_NAME(inlineFASTSHOCKGS, inlineFASTSHOCKGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FASTSHOCKGS");
 FORTRAN_NAME( FASTSHOCKGS ,fastshockgs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFASTSHOCKGS 

#ifndef GUARDSTEADYFASTSHOCK 
#define GUARDSTEADYFASTSHOCK 
// Prototype for Fortran procedure SteadyFastShock ...
//
void FORTRAN_NAME( STEADYFASTSHOCK ,steadyfastshock )(
      CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(phi_UL)
      ,CHFp_CONST_REAL(BL)
      ,CHFp_CONST_REAL(phi_BL)
      ,CHFp_CONST_REAL(QL)
      ,CHFp_CONST_REAL(M_0xL)
      ,CHFp_CONST_REAL(gamma)
      ,CHFp_REAL(RhoL)
      ,CHFp_REAL(u_xL)
      ,CHFp_REAL(u_yL)
      ,CHFp_REAL(B_xL)
      ,CHFp_REAL(B_yL)
      ,CHFp_REAL(pL)
      ,CHFp_REAL(RhoR)
      ,CHFp_REAL(u_xR)
      ,CHFp_REAL(u_yR)
      ,CHFp_REAL(B_xR)
      ,CHFp_REAL(B_yR)
      ,CHFp_REAL(pR) );

#define FORT_STEADYFASTSHOCK FORTRAN_NAME( inlineSTEADYFASTSHOCK, inlineSTEADYFASTSHOCK)
#define FORTNT_STEADYFASTSHOCK FORTRAN_NAME( STEADYFASTSHOCK, steadyfastshock)

inline void FORTRAN_NAME(inlineSTEADYFASTSHOCK, inlineSTEADYFASTSHOCK)(
      CHFp_CONST_REAL(UL)
      ,CHFp_CONST_REAL(phi_UL)
      ,CHFp_CONST_REAL(BL)
      ,CHFp_CONST_REAL(phi_BL)
      ,CHFp_CONST_REAL(QL)
      ,CHFp_CONST_REAL(M_0xL)
      ,CHFp_CONST_REAL(gamma)
      ,CHFp_REAL(RhoL)
      ,CHFp_REAL(u_xL)
      ,CHFp_REAL(u_yL)
      ,CHFp_REAL(B_xL)
      ,CHFp_REAL(B_yL)
      ,CHFp_REAL(pL)
      ,CHFp_REAL(RhoR)
      ,CHFp_REAL(u_xR)
      ,CHFp_REAL(u_yR)
      ,CHFp_REAL(B_xR)
      ,CHFp_REAL(B_yR)
      ,CHFp_REAL(pR) )
{
 CH_TIMELEAF("FORT_STEADYFASTSHOCK");
 FORTRAN_NAME( STEADYFASTSHOCK ,steadyfastshock )(
      CHFt_CONST_REAL(UL)
      ,CHFt_CONST_REAL(phi_UL)
      ,CHFt_CONST_REAL(BL)
      ,CHFt_CONST_REAL(phi_BL)
      ,CHFt_CONST_REAL(QL)
      ,CHFt_CONST_REAL(M_0xL)
      ,CHFt_CONST_REAL(gamma)
      ,CHFt_REAL(RhoL)
      ,CHFt_REAL(u_xL)
      ,CHFt_REAL(u_yL)
      ,CHFt_REAL(B_xL)
      ,CHFt_REAL(B_yL)
      ,CHFt_REAL(pL)
      ,CHFt_REAL(RhoR)
      ,CHFt_REAL(u_xR)
      ,CHFt_REAL(u_yR)
      ,CHFt_REAL(B_xR)
      ,CHFt_REAL(B_yR)
      ,CHFt_REAL(pR) );
}
#endif  // GUARDSTEADYFASTSHOCK 

#ifndef GUARDPOCUBIC 
#define GUARDPOCUBIC 
// Prototype for Fortran procedure poCubic ...
//
void FORTRAN_NAME( POCUBIC ,pocubic )(
      CHFp_CONST_REAL(a0)
      ,CHFp_CONST_REAL(a1)
      ,CHFp_CONST_REAL(a2)
      ,CHFp_CONST_REAL(a3)
      ,CHFp_REAL(x1r)
      ,CHFp_REAL(x1i)
      ,CHFp_REAL(x2r)
      ,CHFp_REAL(x2i)
      ,CHFp_REAL(x3r)
      ,CHFp_REAL(x3i)
      ,CHFp_INT(iRetCode) );

#define FORT_POCUBIC FORTRAN_NAME( inlinePOCUBIC, inlinePOCUBIC)
#define FORTNT_POCUBIC FORTRAN_NAME( POCUBIC, pocubic)

inline void FORTRAN_NAME(inlinePOCUBIC, inlinePOCUBIC)(
      CHFp_CONST_REAL(a0)
      ,CHFp_CONST_REAL(a1)
      ,CHFp_CONST_REAL(a2)
      ,CHFp_CONST_REAL(a3)
      ,CHFp_REAL(x1r)
      ,CHFp_REAL(x1i)
      ,CHFp_REAL(x2r)
      ,CHFp_REAL(x2i)
      ,CHFp_REAL(x3r)
      ,CHFp_REAL(x3i)
      ,CHFp_INT(iRetCode) )
{
 CH_TIMELEAF("FORT_POCUBIC");
 FORTRAN_NAME( POCUBIC ,pocubic )(
      CHFt_CONST_REAL(a0)
      ,CHFt_CONST_REAL(a1)
      ,CHFt_CONST_REAL(a2)
      ,CHFt_CONST_REAL(a3)
      ,CHFt_REAL(x1r)
      ,CHFt_REAL(x1i)
      ,CHFt_REAL(x2r)
      ,CHFt_REAL(x2i)
      ,CHFt_REAL(x3r)
      ,CHFt_REAL(x3i)
      ,CHFt_INT(iRetCode) );
}
#endif  // GUARDPOCUBIC 

#ifndef GUARDPOQUADRATIC 
#define GUARDPOQUADRATIC 
// Prototype for Fortran procedure poQuadratic ...
//
void FORTRAN_NAME( POQUADRATIC ,poquadratic )(
      CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_REAL(c)
      ,CHFp_REAL(ReSmallRoot)
      ,CHFp_REAL(ImSmallRoot)
      ,CHFp_REAL(ReLargeRoot)
      ,CHFp_REAL(ImLargeRoot)
      ,CHFp_INT(iRetCode) );

#define FORT_POQUADRATIC FORTRAN_NAME( inlinePOQUADRATIC, inlinePOQUADRATIC)
#define FORTNT_POQUADRATIC FORTRAN_NAME( POQUADRATIC, poquadratic)

inline void FORTRAN_NAME(inlinePOQUADRATIC, inlinePOQUADRATIC)(
      CHFp_CONST_REAL(a)
      ,CHFp_CONST_REAL(b)
      ,CHFp_CONST_REAL(c)
      ,CHFp_REAL(ReSmallRoot)
      ,CHFp_REAL(ImSmallRoot)
      ,CHFp_REAL(ReLargeRoot)
      ,CHFp_REAL(ImLargeRoot)
      ,CHFp_INT(iRetCode) )
{
 CH_TIMELEAF("FORT_POQUADRATIC");
 FORTRAN_NAME( POQUADRATIC ,poquadratic )(
      CHFt_CONST_REAL(a)
      ,CHFt_CONST_REAL(b)
      ,CHFt_CONST_REAL(c)
      ,CHFt_REAL(ReSmallRoot)
      ,CHFt_REAL(ImSmallRoot)
      ,CHFt_REAL(ReLargeRoot)
      ,CHFt_REAL(ImLargeRoot)
      ,CHFt_INT(iRetCode) );
}
#endif  // GUARDPOQUADRATIC 

}

#endif
