#ifndef _SHOCKSF_F_H_
#define _SHOCKSF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSHOCKBOUNDARY 
#define GUARDSHOCKBOUNDARY 
// Prototype for Fortran procedure SHOCKBOUNDARY ...
//
void FORTRAN_NAME( SHOCKBOUNDARY ,shockboundary )(
      CHFp_CONST_FIA1(REG)
      ,CHFp_FIA1(boundary)
      ,CHFp_BOX(box) );

#define FORT_SHOCKBOUNDARY FORTRAN_NAME( inlineSHOCKBOUNDARY, inlineSHOCKBOUNDARY)
#define FORTNT_SHOCKBOUNDARY FORTRAN_NAME( SHOCKBOUNDARY, shockboundary)

inline void FORTRAN_NAME(inlineSHOCKBOUNDARY, inlineSHOCKBOUNDARY)(
      CHFp_CONST_FIA1(REG)
      ,CHFp_FIA1(boundary)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SHOCKBOUNDARY");
 FORTRAN_NAME( SHOCKBOUNDARY ,shockboundary )(
      CHFt_CONST_FIA1(REG)
      ,CHFt_FIA1(boundary)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSHOCKBOUNDARY 

#ifndef GUARDSHOCKNORMAL 
#define GUARDSHOCKNORMAL 
// Prototype for Fortran procedure SHOCKNORMAL ...
//
void FORTRAN_NAME( SHOCKNORMAL ,shocknormal )(
      CHFp_CONST_FRA(W)
      ,CHFp_BOX(box)
      ,CHFp_FRA(normal)
      ,CHFp_CONST_FIA1(boundary) );

#define FORT_SHOCKNORMAL FORTRAN_NAME( inlineSHOCKNORMAL, inlineSHOCKNORMAL)
#define FORTNT_SHOCKNORMAL FORTRAN_NAME( SHOCKNORMAL, shocknormal)

inline void FORTRAN_NAME(inlineSHOCKNORMAL, inlineSHOCKNORMAL)(
      CHFp_CONST_FRA(W)
      ,CHFp_BOX(box)
      ,CHFp_FRA(normal)
      ,CHFp_CONST_FIA1(boundary) )
{
 CH_TIMELEAF("FORT_SHOCKNORMAL");
 FORTRAN_NAME( SHOCKNORMAL ,shocknormal )(
      CHFt_CONST_FRA(W)
      ,CHFt_BOX(box)
      ,CHFt_FRA(normal)
      ,CHFt_CONST_FIA1(boundary) );
}
#endif  // GUARDSHOCKNORMAL 

#ifndef GUARDSTRONGSCATTER 
#define GUARDSTRONGSCATTER 
// Prototype for Fortran procedure STRONGSCATTER ...
//
void FORTRAN_NAME( STRONGSCATTER ,strongscatter )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_STRONGSCATTER FORTRAN_NAME( inlineSTRONGSCATTER, inlineSTRONGSCATTER)
#define FORTNT_STRONGSCATTER FORTRAN_NAME( STRONGSCATTER, strongscatter)

inline void FORTRAN_NAME(inlineSTRONGSCATTER, inlineSTRONGSCATTER)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_STRONGSCATTER");
 FORTRAN_NAME( STRONGSCATTER ,strongscatter )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDSTRONGSCATTER 

#ifndef GUARDZANKSCATTER 
#define GUARDZANKSCATTER 
// Prototype for Fortran procedure ZANKSCATTER ...
//
void FORTRAN_NAME( ZANKSCATTER ,zankscatter )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_ZANKSCATTER FORTRAN_NAME( inlineZANKSCATTER, inlineZANKSCATTER)
#define FORTNT_ZANKSCATTER FORTRAN_NAME( ZANKSCATTER, zankscatter)

inline void FORTRAN_NAME(inlineZANKSCATTER, inlineZANKSCATTER)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_ZANKSCATTER");
 FORTRAN_NAME( ZANKSCATTER ,zankscatter )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDZANKSCATTER 

#ifndef GUARDZANKSCATTER2 
#define GUARDZANKSCATTER2 
// Prototype for Fortran procedure ZANKSCATTER2 ...
//
void FORTRAN_NAME( ZANKSCATTER2 ,zankscatter2 )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_ZANKSCATTER2 FORTRAN_NAME( inlineZANKSCATTER2, inlineZANKSCATTER2)
#define FORTNT_ZANKSCATTER2 FORTRAN_NAME( ZANKSCATTER2, zankscatter2)

inline void FORTRAN_NAME(inlineZANKSCATTER2, inlineZANKSCATTER2)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_ZANKSCATTER2");
 FORTRAN_NAME( ZANKSCATTER2 ,zankscatter2 )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDZANKSCATTER2 

#ifndef GUARDFLORSCATTER 
#define GUARDFLORSCATTER 
// Prototype for Fortran procedure FLORSCATTER ...
//
void FORTRAN_NAME( FLORSCATTER ,florscatter )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_FLORSCATTER FORTRAN_NAME( inlineFLORSCATTER, inlineFLORSCATTER)
#define FORTNT_FLORSCATTER FORTRAN_NAME( FLORSCATTER, florscatter)

inline void FORTRAN_NAME(inlineFLORSCATTER, inlineFLORSCATTER)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_FLORSCATTER");
 FORTRAN_NAME( FLORSCATTER ,florscatter )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDFLORSCATTER 

#ifndef GUARDVADIMBC 
#define GUARDVADIMBC 
// Prototype for Fortran procedure VADIMBC ...
//
void FORTRAN_NAME( VADIMBC ,vadimbc )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_VADIMBC FORTRAN_NAME( inlineVADIMBC, inlineVADIMBC)
#define FORTNT_VADIMBC FORTRAN_NAME( VADIMBC, vadimbc)

inline void FORTRAN_NAME(inlineVADIMBC, inlineVADIMBC)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_VADIMBC");
 FORTRAN_NAME( VADIMBC ,vadimbc )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDVADIMBC 

#ifndef GUARDVADIMBC_TABLE 
#define GUARDVADIMBC_TABLE 
// Prototype for Fortran procedure VADIMBC_TABLE ...
//
void FORTRAN_NAME( VADIMBC_TABLE ,vadimbc_table )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_VADIMBC_TABLE FORTRAN_NAME( inlineVADIMBC_TABLE, inlineVADIMBC_TABLE)
#define FORTNT_VADIMBC_TABLE FORTRAN_NAME( VADIMBC_TABLE, vadimbc_table)

inline void FORTRAN_NAME(inlineVADIMBC_TABLE, inlineVADIMBC_TABLE)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoPIW)
      ,CHFp_CONST_INT(iRhoPIU)
      ,CHFp_CONST_FIA1(boundary)
      ,CHFp_CONST_FRA(normal)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_VADIMBC_TABLE");
 FORTRAN_NAME( VADIMBC_TABLE ,vadimbc_table )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoPIW)
      ,CHFt_CONST_INT(iRhoPIU)
      ,CHFt_CONST_FIA1(boundary)
      ,CHFt_CONST_FRA(normal)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDVADIMBC_TABLE 

#ifndef GUARDLOAD_TABLE_BC 
#define GUARDLOAD_TABLE_BC 
// Prototype for Fortran procedure LOAD_TABLE_BC ...
//
void FORTRAN_NAME( LOAD_TABLE_BC ,load_table_bc )( );

#define FORT_LOAD_TABLE_BC FORTRAN_NAME( inlineLOAD_TABLE_BC, inlineLOAD_TABLE_BC)
#define FORTNT_LOAD_TABLE_BC FORTRAN_NAME( LOAD_TABLE_BC, load_table_bc)

inline void FORTRAN_NAME(inlineLOAD_TABLE_BC, inlineLOAD_TABLE_BC)( )
{
 CH_TIMELEAF("FORT_LOAD_TABLE_BC");
 FORTRAN_NAME( LOAD_TABLE_BC ,load_table_bc )( );
}
#endif  // GUARDLOAD_TABLE_BC 

}

#endif
