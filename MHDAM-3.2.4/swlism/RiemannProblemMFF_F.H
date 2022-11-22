#ifndef _RIEMANNPROBLEMMFF_F_H_
#define _RIEMANNPROBLEMMFF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETRIEMANNMF 
#define GUARDSETRIEMANNMF 
// Prototype for Fortran procedure SETRIEMANNMF ...
//
void FORTRAN_NAME( SETRIEMANNMF ,setriemannmf )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auNumdenL)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_INT(auNInitDistr) );

#define FORT_SETRIEMANNMF FORTRAN_NAME( inlineSETRIEMANNMF, inlineSETRIEMANNMF)
#define FORTNT_SETRIEMANNMF FORTRAN_NAME( SETRIEMANNMF, setriemannmf)

inline void FORTRAN_NAME(inlineSETRIEMANNMF, inlineSETRIEMANNMF)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auNumdenL)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_INT(auNInitDistr) )
{
 CH_TIMELEAF("FORT_SETRIEMANNMF");
 FORTRAN_NAME( SETRIEMANNMF ,setriemannmf )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auNumdenL)
      ,CHFt_CONST_REAL(auTemperatureL)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auBxL)
      ,CHFt_CONST_REAL(auNetnumL)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_INT(auNInitDistr) );
}
#endif  // GUARDSETRIEMANNMF 

#ifndef GUARDSETRIEMANNLR 
#define GUARDSETRIEMANNLR 
// Prototype for Fortran procedure SETRIEMANNLR ...
//
void FORTRAN_NAME( SETRIEMANNLR ,setriemannlr )(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auTemperatureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auNetTempL)
      ,CHFp_CONST_REAL(aunetVelx)
      ,CHFp_CONST_REAL(aunetVely)
      ,CHFp_CONST_REAL(aunetVelz)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_INT(auNInitDistr)
      ,CHFp_CONST_REAL(auTMLIM) );

#define FORT_SETRIEMANNLR FORTRAN_NAME( inlineSETRIEMANNLR, inlineSETRIEMANNLR)
#define FORTNT_SETRIEMANNLR FORTRAN_NAME( SETRIEMANNLR, setriemannlr)

inline void FORTRAN_NAME(inlineSETRIEMANNLR, inlineSETRIEMANNLR)(
      CHFp_CONST_REAL(augamma)
      ,CHFp_CONST_REAL(auDensityL)
      ,CHFp_CONST_REAL(auDensityR)
      ,CHFp_CONST_REAL(auTemperatureL)
      ,CHFp_CONST_REAL(auTemperatureR)
      ,CHFp_CONST_REAL(auVelxL)
      ,CHFp_CONST_REAL(auVelxR)
      ,CHFp_CONST_REAL(auVelyL)
      ,CHFp_CONST_REAL(auVelyR)
      ,CHFp_CONST_REAL(auVelzL)
      ,CHFp_CONST_REAL(auVelzR)
      ,CHFp_CONST_REAL(auBxL)
      ,CHFp_CONST_REAL(auBxR)
      ,CHFp_CONST_REAL(auByL)
      ,CHFp_CONST_REAL(auByR)
      ,CHFp_CONST_REAL(auBzL)
      ,CHFp_CONST_REAL(auBzR)
      ,CHFp_CONST_REAL(auNetnumL)
      ,CHFp_CONST_REAL(auNetTempL)
      ,CHFp_CONST_REAL(aunetVelx)
      ,CHFp_CONST_REAL(aunetVely)
      ,CHFp_CONST_REAL(aunetVelz)
      ,CHFp_CONST_REAL(auXs)
      ,CHFp_CONST_INT(auNInitDistr)
      ,CHFp_CONST_REAL(auTMLIM) )
{
 CH_TIMELEAF("FORT_SETRIEMANNLR");
 FORTRAN_NAME( SETRIEMANNLR ,setriemannlr )(
      CHFt_CONST_REAL(augamma)
      ,CHFt_CONST_REAL(auDensityL)
      ,CHFt_CONST_REAL(auDensityR)
      ,CHFt_CONST_REAL(auTemperatureL)
      ,CHFt_CONST_REAL(auTemperatureR)
      ,CHFt_CONST_REAL(auVelxL)
      ,CHFt_CONST_REAL(auVelxR)
      ,CHFt_CONST_REAL(auVelyL)
      ,CHFt_CONST_REAL(auVelyR)
      ,CHFt_CONST_REAL(auVelzL)
      ,CHFt_CONST_REAL(auVelzR)
      ,CHFt_CONST_REAL(auBxL)
      ,CHFt_CONST_REAL(auBxR)
      ,CHFt_CONST_REAL(auByL)
      ,CHFt_CONST_REAL(auByR)
      ,CHFt_CONST_REAL(auBzL)
      ,CHFt_CONST_REAL(auBzR)
      ,CHFt_CONST_REAL(auNetnumL)
      ,CHFt_CONST_REAL(auNetTempL)
      ,CHFt_CONST_REAL(aunetVelx)
      ,CHFt_CONST_REAL(aunetVely)
      ,CHFt_CONST_REAL(aunetVelz)
      ,CHFt_CONST_REAL(auXs)
      ,CHFt_CONST_INT(auNInitDistr)
      ,CHFt_CONST_REAL(auTMLIM) );
}
#endif  // GUARDSETRIEMANNLR 

#ifndef GUARDRIEMANNINITMF 
#define GUARDRIEMANNINITMF 
// Prototype for Fortran procedure RIEMANNINITMF ...
//
void FORTRAN_NAME( RIEMANNINITMF ,riemanninitmf )(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_RIEMANNINITMF FORTRAN_NAME( inlineRIEMANNINITMF, inlineRIEMANNINITMF)
#define FORTNT_RIEMANNINITMF FORTRAN_NAME( RIEMANNINITMF, riemanninitmf)

inline void FORTRAN_NAME(inlineRIEMANNINITMF, inlineRIEMANNINITMF)(
      CHFp_FRA(U)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_RIEMANNINITMF");
 FORTRAN_NAME( RIEMANNINITMF ,riemanninitmf )(
      CHFt_FRA(U)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDRIEMANNINITMF 

#ifndef GUARDFLUXBCMF 
#define GUARDFLUXBCMF 
// Prototype for Fortran procedure FLUXBCMF ...
//
void FORTRAN_NAME( FLUXBCMF ,fluxbcmf )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_FLUXBCMF FORTRAN_NAME( inlineFLUXBCMF, inlineFLUXBCMF)
#define FORTNT_FLUXBCMF FORTRAN_NAME( FLUXBCMF, fluxbcmf)

inline void FORTRAN_NAME(inlineFLUXBCMF, inlineFLUXBCMF)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXBCMF");
 FORTRAN_NAME( FLUXBCMF ,fluxbcmf )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXBCMF 

#ifndef GUARDDEFINE_REGIONS_3F 
#define GUARDDEFINE_REGIONS_3F 
// Prototype for Fortran procedure DEFINE_REGIONS_3F ...
//
void FORTRAN_NAME( DEFINE_REGIONS_3F ,define_regions_3f )(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_BOX(box) );

#define FORT_DEFINE_REGIONS_3F FORTRAN_NAME( inlineDEFINE_REGIONS_3F, inlineDEFINE_REGIONS_3F)
#define FORTNT_DEFINE_REGIONS_3F FORTRAN_NAME( DEFINE_REGIONS_3F, define_regions_3f)

inline void FORTRAN_NAME(inlineDEFINE_REGIONS_3F, inlineDEFINE_REGIONS_3F)(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEFINE_REGIONS_3F");
 FORTRAN_NAME( DEFINE_REGIONS_3F ,define_regions_3f )(
      CHFt_CONST_FRA(W)
      ,CHFt_FIA1(REG)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEFINE_REGIONS_3F 

}

#endif
