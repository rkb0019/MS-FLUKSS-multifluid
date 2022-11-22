#ifndef _DEDNERF_F_H_
#define _DEDNERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETDEDNERCONST 
#define GUARDSETDEDNERCONST 
// Prototype for Fortran procedure SETDEDNERCONST ...
//
void FORTRAN_NAME( SETDEDNERCONST ,setdednerconst )(
      CHFp_CONST_REAL(auCh)
      ,CHFp_CONST_REAL(auCp) );

#define FORT_SETDEDNERCONST FORTRAN_NAME( inlineSETDEDNERCONST, inlineSETDEDNERCONST)
#define FORTNT_SETDEDNERCONST FORTRAN_NAME( SETDEDNERCONST, setdednerconst)

inline void FORTRAN_NAME(inlineSETDEDNERCONST, inlineSETDEDNERCONST)(
      CHFp_CONST_REAL(auCh)
      ,CHFp_CONST_REAL(auCp) )
{
 CH_TIMELEAF("FORT_SETDEDNERCONST");
 FORTRAN_NAME( SETDEDNERCONST ,setdednerconst )(
      CHFt_CONST_REAL(auCh)
      ,CHFt_CONST_REAL(auCp) );
}
#endif  // GUARDSETDEDNERCONST 

#ifndef GUARDDEDNERSOURCETERMS 
#define GUARDDEDNERSOURCETERMS 
// Prototype for Fortran procedure DEDNERSOURCETERMS ...
//
void FORTRAN_NAME( DEDNERSOURCETERMS ,dednersourceterms )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_DEDNERSOURCETERMS FORTRAN_NAME( inlineDEDNERSOURCETERMS, inlineDEDNERSOURCETERMS)
#define FORTNT_DEDNERSOURCETERMS FORTRAN_NAME( DEDNERSOURCETERMS, dednersourceterms)

inline void FORTRAN_NAME(inlineDEDNERSOURCETERMS, inlineDEDNERSOURCETERMS)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_FRA1(divBFra)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEDNERSOURCETERMS");
 FORTRAN_NAME( DEDNERSOURCETERMS ,dednersourceterms )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_FRA1(divBFra)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEDNERSOURCETERMS 

#ifndef GUARDPSIPARABOLIC 
#define GUARDPSIPARABOLIC 
// Prototype for Fortran procedure PSIPARABOLIC ...
//
void FORTRAN_NAME( PSIPARABOLIC ,psiparabolic )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_PSIPARABOLIC FORTRAN_NAME( inlinePSIPARABOLIC, inlinePSIPARABOLIC)
#define FORTNT_PSIPARABOLIC FORTRAN_NAME( PSIPARABOLIC, psiparabolic)

inline void FORTRAN_NAME(inlinePSIPARABOLIC, inlinePSIPARABOLIC)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PSIPARABOLIC");
 FORTRAN_NAME( PSIPARABOLIC ,psiparabolic )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPSIPARABOLIC 

#ifndef GUARDDEDNER2X2F 
#define GUARDDEDNER2X2F 
// Prototype for Fortran procedure DEDNER2X2F ...
//
void FORTRAN_NAME( DEDNER2X2F ,dedner2x2f )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_DEDNER2X2F FORTRAN_NAME( inlineDEDNER2X2F, inlineDEDNER2X2F)
#define FORTNT_DEDNER2X2F FORTRAN_NAME( DEDNER2X2F, dedner2x2f)

inline void FORTRAN_NAME(inlineDEDNER2X2F, inlineDEDNER2X2F)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEDNER2X2F");
 FORTRAN_NAME( DEDNER2X2F ,dedner2x2f )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEDNER2X2F 

#ifndef GUARDDEDNER9X9F 
#define GUARDDEDNER9X9F 
// Prototype for Fortran procedure DEDNER9X9F ...
//
void FORTRAN_NAME( DEDNER9X9F ,dedner9x9f )(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) );

#define FORT_DEDNER9X9F FORTRAN_NAME( inlineDEDNER9X9F, inlineDEDNER9X9F)
#define FORTNT_DEDNER9X9F FORTRAN_NAME( DEDNER9X9F, dedner9x9f)

inline void FORTRAN_NAME(inlineDEDNER9X9F, inlineDEDNER9X9F)(
      CHFp_FRA(Wgdnv)
      ,CHFp_CONST_FRA(WLeft)
      ,CHFp_CONST_FRA(WRight)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRho)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEDNER9X9F");
 FORTRAN_NAME( DEDNER9X9F ,dedner9x9f )(
      CHFt_FRA(Wgdnv)
      ,CHFt_CONST_FRA(WLeft)
      ,CHFt_CONST_FRA(WRight)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRho)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEDNER9X9F 

#ifndef GUARDSOURCEAXISYMMETRICDEDNER 
#define GUARDSOURCEAXISYMMETRICDEDNER 
// Prototype for Fortran procedure SOURCEAXISYMMETRICDEDNER ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRICDEDNER ,sourceaxisymmetricdedner )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRICDEDNER FORTRAN_NAME( inlineSOURCEAXISYMMETRICDEDNER, inlineSOURCEAXISYMMETRICDEDNER)
#define FORTNT_SOURCEAXISYMMETRICDEDNER FORTRAN_NAME( SOURCEAXISYMMETRICDEDNER, sourceaxisymmetricdedner)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRICDEDNER, inlineSOURCEAXISYMMETRICDEDNER)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRICDEDNER");
 FORTRAN_NAME( SOURCEAXISYMMETRICDEDNER ,sourceaxisymmetricdedner )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRICDEDNER 

}

#endif
