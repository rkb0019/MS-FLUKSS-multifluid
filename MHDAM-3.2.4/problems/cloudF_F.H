#ifndef _CLOUDF_F_H_
#define _CLOUDF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETCLOUD 
#define GUARDSETCLOUD 
// Prototype for Fortran procedure SETCLOUD ...
//
void FORTRAN_NAME( SETCLOUD ,setcloud )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aucloudRL)
      ,CHFp_CONST_REAL(aucloudUL)
      ,CHFp_CONST_REAL(aucloudVL)
      ,CHFp_CONST_REAL(aucloudWL)
      ,CHFp_CONST_REAL(aucloudPL)
      ,CHFp_CONST_REAL(aucloudBXL)
      ,CHFp_CONST_REAL(aucloudBYL)
      ,CHFp_CONST_REAL(aucloudBZL)
      ,CHFp_CONST_REAL(aucloudRR)
      ,CHFp_CONST_REAL(aucloudUR)
      ,CHFp_CONST_REAL(aucloudPR)
      ,CHFp_CONST_REAL(aucloudBYR)
      ,CHFp_CONST_REAL(aucloudBZR)
      ,CHFp_CONST_REAL(aucloudXS)
      ,CHFp_CONST_REAL(aucloudXC)
      ,CHFp_CONST_REAL(aucloudYC)
      ,CHFp_CONST_REAL(aucloudZC)
      ,CHFp_CONST_REAL(aucloudR0)
      ,CHFp_CONST_REAL(aucloudRho) );

#define FORT_SETCLOUD FORTRAN_NAME( inlineSETCLOUD, inlineSETCLOUD)
#define FORTNT_SETCLOUD FORTRAN_NAME( SETCLOUD, setcloud)

inline void FORTRAN_NAME(inlineSETCLOUD, inlineSETCLOUD)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(aucloudRL)
      ,CHFp_CONST_REAL(aucloudUL)
      ,CHFp_CONST_REAL(aucloudVL)
      ,CHFp_CONST_REAL(aucloudWL)
      ,CHFp_CONST_REAL(aucloudPL)
      ,CHFp_CONST_REAL(aucloudBXL)
      ,CHFp_CONST_REAL(aucloudBYL)
      ,CHFp_CONST_REAL(aucloudBZL)
      ,CHFp_CONST_REAL(aucloudRR)
      ,CHFp_CONST_REAL(aucloudUR)
      ,CHFp_CONST_REAL(aucloudPR)
      ,CHFp_CONST_REAL(aucloudBYR)
      ,CHFp_CONST_REAL(aucloudBZR)
      ,CHFp_CONST_REAL(aucloudXS)
      ,CHFp_CONST_REAL(aucloudXC)
      ,CHFp_CONST_REAL(aucloudYC)
      ,CHFp_CONST_REAL(aucloudZC)
      ,CHFp_CONST_REAL(aucloudR0)
      ,CHFp_CONST_REAL(aucloudRho) )
{
 CH_TIMELEAF("FORT_SETCLOUD");
 FORTRAN_NAME( SETCLOUD ,setcloud )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(aucloudRL)
      ,CHFt_CONST_REAL(aucloudUL)
      ,CHFt_CONST_REAL(aucloudVL)
      ,CHFt_CONST_REAL(aucloudWL)
      ,CHFt_CONST_REAL(aucloudPL)
      ,CHFt_CONST_REAL(aucloudBXL)
      ,CHFt_CONST_REAL(aucloudBYL)
      ,CHFt_CONST_REAL(aucloudBZL)
      ,CHFt_CONST_REAL(aucloudRR)
      ,CHFt_CONST_REAL(aucloudUR)
      ,CHFt_CONST_REAL(aucloudPR)
      ,CHFt_CONST_REAL(aucloudBYR)
      ,CHFt_CONST_REAL(aucloudBZR)
      ,CHFt_CONST_REAL(aucloudXS)
      ,CHFt_CONST_REAL(aucloudXC)
      ,CHFt_CONST_REAL(aucloudYC)
      ,CHFt_CONST_REAL(aucloudZC)
      ,CHFt_CONST_REAL(aucloudR0)
      ,CHFt_CONST_REAL(aucloudRho) );
}
#endif  // GUARDSETCLOUD 

#ifndef GUARDCLOUDINIT 
#define GUARDCLOUDINIT 
// Prototype for Fortran procedure CLOUDINIT ...
//
void FORTRAN_NAME( CLOUDINIT ,cloudinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_CLOUDINIT FORTRAN_NAME( inlineCLOUDINIT, inlineCLOUDINIT)
#define FORTNT_CLOUDINIT FORTRAN_NAME( CLOUDINIT, cloudinit)

inline void FORTRAN_NAME(inlineCLOUDINIT, inlineCLOUDINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CLOUDINIT");
 FORTRAN_NAME( CLOUDINIT ,cloudinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCLOUDINIT 

#ifndef GUARDCLOUDBC 
#define GUARDCLOUDBC 
// Prototype for Fortran procedure CLOUDBC ...
//
void FORTRAN_NAME( CLOUDBC ,cloudbc )(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_CLOUDBC FORTRAN_NAME( inlineCLOUDBC, inlineCLOUDBC)
#define FORTNT_CLOUDBC FORTRAN_NAME( CLOUDBC, cloudbc)

inline void FORTRAN_NAME(inlineCLOUDBC, inlineCLOUDBC)(
      CHFp_FRA(F)
      ,CHFp_FRA1(Bn)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CLOUDBC");
 FORTRAN_NAME( CLOUDBC ,cloudbc )(
      CHFt_FRA(F)
      ,CHFt_FRA1(Bn)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCLOUDBC 

#ifndef GUARDCLOUDGS 
#define GUARDCLOUDGS 
// Prototype for Fortran procedure CLOUDGS ...
//
void FORTRAN_NAME( CLOUDGS ,cloudgs )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_BOX(box) );

#define FORT_CLOUDGS FORTRAN_NAME( inlineCLOUDGS, inlineCLOUDGS)
#define FORTNT_CLOUDGS FORTRAN_NAME( CLOUDGS, cloudgs)

inline void FORTRAN_NAME(inlineCLOUDGS, inlineCLOUDGS)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iCP)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CLOUDGS");
 FORTRAN_NAME( CLOUDGS ,cloudgs )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iCP)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCLOUDGS 

}

#endif
