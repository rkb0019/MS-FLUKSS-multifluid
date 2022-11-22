#ifndef _PATCHEULERF_F_H_
#define _PATCHEULERF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDMAXWAVESPEED_E 
#define GUARDMAXWAVESPEED_E 
// Prototype for Fortran procedure MAXWAVESPEED_E ...
//
void FORTRAN_NAME( MAXWAVESPEED_E ,maxwavespeed_e )(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) );

#define FORT_MAXWAVESPEED_E FORTRAN_NAME( inlineMAXWAVESPEED_E, inlineMAXWAVESPEED_E)
#define FORTNT_MAXWAVESPEED_E FORTRAN_NAME( MAXWAVESPEED_E, maxwavespeed_e)

inline void FORTRAN_NAME(inlineMAXWAVESPEED_E, inlineMAXWAVESPEED_E)(
      CHFp_FRA1(cmax)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MAXWAVESPEED_E");
 FORTRAN_NAME( MAXWAVESPEED_E ,maxwavespeed_e )(
      CHFt_FRA1(cmax)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMAXWAVESPEED_E 

#ifndef GUARDMINDT_SPHERICAL_E 
#define GUARDMINDT_SPHERICAL_E 
// Prototype for Fortran procedure MINDT_SPHERICAL_E ...
//
void FORTRAN_NAME( MINDT_SPHERICAL_E ,mindt_spherical_e )(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_MINDT_SPHERICAL_E FORTRAN_NAME( inlineMINDT_SPHERICAL_E, inlineMINDT_SPHERICAL_E)
#define FORTNT_MINDT_SPHERICAL_E FORTRAN_NAME( MINDT_SPHERICAL_E, mindt_spherical_e)

inline void FORTRAN_NAME(inlineMINDT_SPHERICAL_E, inlineMINDT_SPHERICAL_E)(
      CHFp_FRA1(dt)
      ,CHFp_CONST_FRA(cons)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINDT_SPHERICAL_E");
 FORTRAN_NAME( MINDT_SPHERICAL_E ,mindt_spherical_e )(
      CHFt_FRA1(dt)
      ,CHFt_CONST_FRA(cons)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINDT_SPHERICAL_E 

#ifndef GUARDCONSTOPRIM_E 
#define GUARDCONSTOPRIM_E 
// Prototype for Fortran procedure CONSTOPRIM_E ...
//
void FORTRAN_NAME( CONSTOPRIM_E ,constoprim_e )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) );

#define FORT_CONSTOPRIM_E FORTRAN_NAME( inlineCONSTOPRIM_E, inlineCONSTOPRIM_E)
#define FORTNT_CONSTOPRIM_E FORTRAN_NAME( CONSTOPRIM_E, constoprim_e)

inline void FORTRAN_NAME(inlineCONSTOPRIM_E, inlineCONSTOPRIM_E)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONSTOPRIM_E");
 FORTRAN_NAME( CONSTOPRIM_E ,constoprim_e )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRHOP)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONSTOPRIM_E 

#ifndef GUARDPRIMTOCONS_E 
#define GUARDPRIMTOCONS_E 
// Prototype for Fortran procedure PRIMTOCONS_E ...
//
void FORTRAN_NAME( PRIMTOCONS_E ,primtocons_e )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) );

#define FORT_PRIMTOCONS_E FORTRAN_NAME( inlinePRIMTOCONS_E, inlinePRIMTOCONS_E)
#define FORTNT_PRIMTOCONS_E FORTRAN_NAME( PRIMTOCONS_E, primtocons_e)

inline void FORTRAN_NAME(inlinePRIMTOCONS_E, inlinePRIMTOCONS_E)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PRIMTOCONS_E");
 FORTRAN_NAME( PRIMTOCONS_E ,primtocons_e )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRHOP)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPRIMTOCONS_E 

#ifndef GUARDSOURCEAXISYMMETRIC_E 
#define GUARDSOURCEAXISYMMETRIC_E 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_E ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_E ,sourceaxisymmetric_e )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_E FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_E, inlineSOURCEAXISYMMETRIC_E)
#define FORTNT_SOURCEAXISYMMETRIC_E FORTRAN_NAME( SOURCEAXISYMMETRIC_E, sourceaxisymmetric_e)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_E, inlineSOURCEAXISYMMETRIC_E)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_E");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_E ,sourceaxisymmetric_e )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRHOP)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_E 

#ifndef GUARDSOURCEAXISYMMETRIC_POLAR_E 
#define GUARDSOURCEAXISYMMETRIC_POLAR_E 
// Prototype for Fortran procedure SOURCEAXISYMMETRIC_POLAR_E ...
//
void FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR_E ,sourceaxisymmetric_polar_e )(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SOURCEAXISYMMETRIC_POLAR_E FORTRAN_NAME( inlineSOURCEAXISYMMETRIC_POLAR_E, inlineSOURCEAXISYMMETRIC_POLAR_E)
#define FORTNT_SOURCEAXISYMMETRIC_POLAR_E FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR_E, sourceaxisymmetric_polar_e)

inline void FORTRAN_NAME(inlineSOURCEAXISYMMETRIC_POLAR_E, inlineSOURCEAXISYMMETRIC_POLAR_E)(
      CHFp_FRA(S)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SOURCEAXISYMMETRIC_POLAR_E");
 FORTRAN_NAME( SOURCEAXISYMMETRIC_POLAR_E ,sourceaxisymmetric_polar_e )(
      CHFt_FRA(S)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRHOP)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSOURCEAXISYMMETRIC_POLAR_E 

#ifndef GUARDADDSOURCE_E 
#define GUARDADDSOURCE_E 
// Prototype for Fortran procedure ADDSOURCE_E ...
//
void FORTRAN_NAME( ADDSOURCE_E ,addsource_e )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) );

#define FORT_ADDSOURCE_E FORTRAN_NAME( inlineADDSOURCE_E, inlineADDSOURCE_E)
#define FORTNT_ADDSOURCE_E FORTRAN_NAME( ADDSOURCE_E, addsource_e)

inline void FORTRAN_NAME(inlineADDSOURCE_E, inlineADDSOURCE_E)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADDSOURCE_E");
 FORTRAN_NAME( ADDSOURCE_E ,addsource_e )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(S)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADDSOURCE_E 

#ifndef GUARDPOSTPROCESSING_E 
#define GUARDPOSTPROCESSING_E 
// Prototype for Fortran procedure POSTPROCESSING_E ...
//
void FORTRAN_NAME( POSTPROCESSING_E ,postprocessing_e )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) );

#define FORT_POSTPROCESSING_E FORTRAN_NAME( inlinePOSTPROCESSING_E, inlinePOSTPROCESSING_E)
#define FORTNT_POSTPROCESSING_E FORTRAN_NAME( POSTPROCESSING_E, postprocessing_e)

inline void FORTRAN_NAME(inlinePOSTPROCESSING_E, inlinePOSTPROCESSING_E)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(Uold)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_POSTPROCESSING_E");
 FORTRAN_NAME( POSTPROCESSING_E ,postprocessing_e )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(Uold)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPOSTPROCESSING_E 

#ifndef GUARDFLUXESHANCOCK_E 
#define GUARDFLUXESHANCOCK_E 
// Prototype for Fortran procedure FLUXESHANCOCK_E ...
//
void FORTRAN_NAME( FLUXESHANCOCK_E ,fluxeshancock_e )(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) );

#define FORT_FLUXESHANCOCK_E FORTRAN_NAME( inlineFLUXESHANCOCK_E, inlineFLUXESHANCOCK_E)
#define FORTNT_FLUXESHANCOCK_E FORTRAN_NAME( FLUXESHANCOCK_E, fluxeshancock_e)

inline void FORTRAN_NAME(inlineFLUXESHANCOCK_E, inlineFLUXESHANCOCK_E)(
      CHFp_FRA(F)
      ,CHFp_CONST_FRA(Prim)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iRHO)
      ,CHFp_CONST_INT(iRHOP)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FLUXESHANCOCK_E");
 FORTRAN_NAME( FLUXESHANCOCK_E ,fluxeshancock_e )(
      CHFt_FRA(F)
      ,CHFt_CONST_FRA(Prim)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iRHO)
      ,CHFt_CONST_INT(iRHOP)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFLUXESHANCOCK_E 

#ifndef GUARDCHARANALYSIS_E 
#define GUARDCHARANALYSIS_E 
// Prototype for Fortran procedure CHARANALYSIS_E ...
//
void FORTRAN_NAME( CHARANALYSIS_E ,charanalysis_e )(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_CHARANALYSIS_E FORTRAN_NAME( inlineCHARANALYSIS_E, inlineCHARANALYSIS_E)
#define FORTNT_CHARANALYSIS_E FORTRAN_NAME( CHARANALYSIS_E, charanalysis_e)

inline void FORTRAN_NAME(inlineCHARANALYSIS_E, inlineCHARANALYSIS_E)(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARANALYSIS_E");
 FORTRAN_NAME( CHARANALYSIS_E ,charanalysis_e )(
      CHFt_FRA(dWLeft)
      ,CHFt_FRA(dWRight)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARANALYSIS_E 

#ifndef GUARDCHARSYNTHESIS_E 
#define GUARDCHARSYNTHESIS_E 
// Prototype for Fortran procedure CHARSYNTHESIS_E ...
//
void FORTRAN_NAME( CHARSYNTHESIS_E ,charsynthesis_e )(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) );

#define FORT_CHARSYNTHESIS_E FORTRAN_NAME( inlineCHARSYNTHESIS_E, inlineCHARSYNTHESIS_E)
#define FORTNT_CHARSYNTHESIS_E FORTRAN_NAME( CHARSYNTHESIS_E, charsynthesis_e)

inline void FORTRAN_NAME(inlineCHARSYNTHESIS_E, inlineCHARSYNTHESIS_E)(
      CHFp_FRA(dWLeft)
      ,CHFp_FRA(dWRight)
      ,CHFp_CONST_FRA(W)
      ,CHFp_CONST_INT(dir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHARSYNTHESIS_E");
 FORTRAN_NAME( CHARSYNTHESIS_E ,charsynthesis_e )(
      CHFt_FRA(dWLeft)
      ,CHFt_FRA(dWRight)
      ,CHFt_CONST_FRA(W)
      ,CHFt_CONST_INT(dir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHARSYNTHESIS_E 

}

#endif
