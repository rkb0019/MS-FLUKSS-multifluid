#ifndef _KINETICSOURCES2DF_F_H_
#define _KINETICSOURCES2DF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDPREPAREDATAFORKINETICSC2D 
#define GUARDPREPAREDATAFORKINETICSC2D 
// Prototype for Fortran procedure PREPAREDATAFORKINETICSC2D ...
//
void FORTRAN_NAME( PREPAREDATAFORKINETICSC2D ,preparedataforkineticsc2d )(
      CHFp_CONST_FRA(U)
      ,CHFp_FRA(Ph)
      ,CHFp_BOX(box) );

#define FORT_PREPAREDATAFORKINETICSC2D FORTRAN_NAME( inlinePREPAREDATAFORKINETICSC2D, inlinePREPAREDATAFORKINETICSC2D)
#define FORTNT_PREPAREDATAFORKINETICSC2D FORTRAN_NAME( PREPAREDATAFORKINETICSC2D, preparedataforkineticsc2d)

inline void FORTRAN_NAME(inlinePREPAREDATAFORKINETICSC2D, inlinePREPAREDATAFORKINETICSC2D)(
      CHFp_CONST_FRA(U)
      ,CHFp_FRA(Ph)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_PREPAREDATAFORKINETICSC2D");
 FORTRAN_NAME( PREPAREDATAFORKINETICSC2D ,preparedataforkineticsc2d )(
      CHFt_CONST_FRA(U)
      ,CHFt_FRA(Ph)
      ,CHFt_BOX(box) );
}
#endif  // GUARDPREPAREDATAFORKINETICSC2D 

#ifndef GUARDMODIFYSOURCEAFTERKINETICSC2D 
#define GUARDMODIFYSOURCEAFTERKINETICSC2D 
// Prototype for Fortran procedure MODIFYSOURCEAFTERKINETICSC2D ...
//
void FORTRAN_NAME( MODIFYSOURCEAFTERKINETICSC2D ,modifysourceafterkineticsc2d )(
      CHFp_FRA(S)
      ,CHFp_BOX(box) );

#define FORT_MODIFYSOURCEAFTERKINETICSC2D FORTRAN_NAME( inlineMODIFYSOURCEAFTERKINETICSC2D, inlineMODIFYSOURCEAFTERKINETICSC2D)
#define FORTNT_MODIFYSOURCEAFTERKINETICSC2D FORTRAN_NAME( MODIFYSOURCEAFTERKINETICSC2D, modifysourceafterkineticsc2d)

inline void FORTRAN_NAME(inlineMODIFYSOURCEAFTERKINETICSC2D, inlineMODIFYSOURCEAFTERKINETICSC2D)(
      CHFp_FRA(S)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MODIFYSOURCEAFTERKINETICSC2D");
 FORTRAN_NAME( MODIFYSOURCEAFTERKINETICSC2D ,modifysourceafterkineticsc2d )(
      CHFt_FRA(S)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMODIFYSOURCEAFTERKINETICSC2D 

#ifndef GUARDADDKINETICSOURCE2D 
#define GUARDADDKINETICSOURCE2D 
// Prototype for Fortran procedure ADDKINETICSOURCE2D ...
//
void FORTRAN_NAME( ADDKINETICSOURCE2D ,addkineticsource2d )(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) );

#define FORT_ADDKINETICSOURCE2D FORTRAN_NAME( inlineADDKINETICSOURCE2D, inlineADDKINETICSOURCE2D)
#define FORTNT_ADDKINETICSOURCE2D FORTRAN_NAME( ADDKINETICSOURCE2D, addkineticsource2d)

inline void FORTRAN_NAME(inlineADDKINETICSOURCE2D, inlineADDKINETICSOURCE2D)(
      CHFp_FRA(U)
      ,CHFp_CONST_FRA(S)
      ,CHFp_CONST_REAL(dt)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ADDKINETICSOURCE2D");
 FORTRAN_NAME( ADDKINETICSOURCE2D ,addkineticsource2d )(
      CHFt_FRA(U)
      ,CHFt_CONST_FRA(S)
      ,CHFt_CONST_REAL(dt)
      ,CHFt_BOX(box) );
}
#endif  // GUARDADDKINETICSOURCE2D 

#ifndef GUARDDEFINE_REGIONS_KINETIC 
#define GUARDDEFINE_REGIONS_KINETIC 
// Prototype for Fortran procedure DEFINE_REGIONS_KINETIC ...
//
void FORTRAN_NAME( DEFINE_REGIONS_KINETIC ,define_regions_kinetic )(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) );

#define FORT_DEFINE_REGIONS_KINETIC FORTRAN_NAME( inlineDEFINE_REGIONS_KINETIC, inlineDEFINE_REGIONS_KINETIC)
#define FORTNT_DEFINE_REGIONS_KINETIC FORTRAN_NAME( DEFINE_REGIONS_KINETIC, define_regions_kinetic)

inline void FORTRAN_NAME(inlineDEFINE_REGIONS_KINETIC, inlineDEFINE_REGIONS_KINETIC)(
      CHFp_CONST_FRA(W)
      ,CHFp_FIA1(REG)
      ,CHFp_CONST_REAL(dx)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_DEFINE_REGIONS_KINETIC");
 FORTRAN_NAME( DEFINE_REGIONS_KINETIC ,define_regions_kinetic )(
      CHFt_CONST_FRA(W)
      ,CHFt_FIA1(REG)
      ,CHFt_CONST_REAL(dx)
      ,CHFt_BOX(box) );
}
#endif  // GUARDDEFINE_REGIONS_KINETIC 

}

#endif
