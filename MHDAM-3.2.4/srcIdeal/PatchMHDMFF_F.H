#ifndef _PATCHMHDMFF_F_H_
#define _PATCHMHDMFF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCONSTOPRIM_MF 
#define GUARDCONSTOPRIM_MF 
// Prototype for Fortran procedure CONSTOPRIM_MF ...
//
void FORTRAN_NAME( CONSTOPRIM_MF ,constoprim_mf )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_CONSTOPRIM_MF FORTRAN_NAME( inlineCONSTOPRIM_MF, inlineCONSTOPRIM_MF)
#define FORTNT_CONSTOPRIM_MF FORTRAN_NAME( CONSTOPRIM_MF, constoprim_mf)

inline void FORTRAN_NAME(inlineCONSTOPRIM_MF, inlineCONSTOPRIM_MF)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONSTOPRIM_MF");
 FORTRAN_NAME( CONSTOPRIM_MF ,constoprim_mf )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONSTOPRIM_MF 

#ifndef GUARDPRIMTOCONS_MF 
#define GUARDPRIMTOCONS_MF 
// Prototype for Fortran procedure PRIMTOCONS_MF ...
//
void FORTRAN_NAME( PRIMTOCONS_MF ,primtocons_mf )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) );

#define FORT_PRIMTOCONS_MF FORTRAN_NAME( inlinePRIMTOCONS_MF, inlinePRIMTOCONS_MF)
#define FORTNT_PRIMTOCONS_MF FORTRAN_NAME( PRIMTOCONS_MF, primtocons_mf)

inline void FORTRAN_NAME(inlinePRIMTOCONS_MF, inlinePRIMTOCONS_MF)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRIMTOCONS_MF");
 FORTRAN_NAME( PRIMTOCONS_MF ,primtocons_mf )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRIMTOCONS_MF 

#ifndef GUARDSOURCEAXISYMMETRIC_MF 
#define GUARDSOURCEAXISYMMETRIC_MF 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_MF ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_MF ,sourceaxisymmetric_mf )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_MF FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_MF, inlineSOURCEAXISYMMETRIC_MF)
#define FORTNT_SOURCEAXISYMMETRIC_MF FORTRAN_NAME( SOURCEAXISYMMETRIC_MF, sourceaxisymmetric_mf)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_MF, inlineSOURCEAXISYMMETRIC_MF)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_MF");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_MF ,sourceaxisymmetric_mf )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_MF 

#ifndef GUARDSOURCEAXISYMMETRIC_CYL_MF 
#define GUARDSOURCEAXISYMMETRIC_CYL_MF 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_CYL_MF ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL_MF ,sourceaxisymmetric_cyl_mf )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_CYL_MF FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_CYL_MF, inlineSOURCEAXISYMMETRIC_CYL_MF)
#define FORTNT_SOURCEAXISYMMETRIC_CYL_MF FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL_MF, sourceaxisymmetric_cyl_mf)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_CYL_MF, inlineSOURCEAXISYMMETRIC_CYL_MF)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_CYL_MF");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_CYL_MF ,sourceaxisymmetric_cyl_mf )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_CYL_MF 

#ifndef GUARDPOSTPROCESSING_MF 
#define GUARDPOSTPROCESSING_MF 
// Prototype for Fortran procedure POSTPROCESSING_MF ...
//
void FORTRAN_NAME( POSTPROCESSING_MF ,postprocessing_mf )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING_MF FORTRAN_NAME( inlinePOSTPROCESSING_MF, inlinePOSTPROCESSING_MF)
#define FORTNT_POSTPROCESSING_MF FORTRAN_NAME( POSTPROCESSING_MF, postprocessing_mf)

inline void FORTRAN_NAME(inlinePOSTPROCESSING_MF, inlinePOSTPROCESSING_MF)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING_MF");
 FORTRAN_NAME( POSTPROCESSING_MF ,postprocessing_mf )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Uold)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING_MF 

}

#endif
