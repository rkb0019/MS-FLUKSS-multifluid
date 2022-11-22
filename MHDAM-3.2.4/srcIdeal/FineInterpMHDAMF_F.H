#ifndef _FINEINTERPMHDAMF_F_H_
#define _FINEINTERPMHDAMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTERPCENTRALDERIV_C 
#define GUARDINTERPCENTRALDERIV_C 
// Prototype for Fortran procedure interpcentralderiv_c ...
//
void FORTRAN_NAME( INTERPCENTRALDERIV_C ,interpcentralderiv_c )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_INTERPCENTRALDERIV_C FORTRAN_NAME( inlineINTERPCENTRALDERIV_C, inlineINTERPCENTRALDERIV_C)
#define FORTNT_INTERPCENTRALDERIV_C FORTRAN_NAME( INTERPCENTRALDERIV_C, interpcentralderiv_c)

inline void FORTRAN_NAME(inlineINTERPCENTRALDERIV_C, inlineINTERPCENTRALDERIV_C)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_INTERPCENTRALDERIV_C");
 FORTRAN_NAME( INTERPCENTRALDERIV_C ,interpcentralderiv_c )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDINTERPCENTRALDERIV_C 

#ifndef GUARDINTERPCENTRALDERIV_V 
#define GUARDINTERPCENTRALDERIV_V 
// Prototype for Fortran procedure interpcentralderiv_v ...
//
void FORTRAN_NAME( INTERPCENTRALDERIV_V ,interpcentralderiv_v )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) );

#define FORT_INTERPCENTRALDERIV_V FORTRAN_NAME( inlineINTERPCENTRALDERIV_V, inlineINTERPCENTRALDERIV_V)
#define FORTNT_INTERPCENTRALDERIV_V FORTRAN_NAME( INTERPCENTRALDERIV_V, interpcentralderiv_v)

inline void FORTRAN_NAME(inlineINTERPCENTRALDERIV_V, inlineINTERPCENTRALDERIV_V)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) )
{
 CH_TIMELEAF("FORT_INTERPCENTRALDERIV_V");
 FORTRAN_NAME( INTERPCENTRALDERIV_V ,interpcentralderiv_v )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(dxarr) );
}
#endif  // GUARDINTERPCENTRALDERIV_V 

#ifndef GUARDINTERPHISIDEDERIV_C 
#define GUARDINTERPHISIDEDERIV_C 
// Prototype for Fortran procedure interphisidederiv_c ...
//
void FORTRAN_NAME( INTERPHISIDEDERIV_C ,interphisidederiv_c )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_INTERPHISIDEDERIV_C FORTRAN_NAME( inlineINTERPHISIDEDERIV_C, inlineINTERPHISIDEDERIV_C)
#define FORTNT_INTERPHISIDEDERIV_C FORTRAN_NAME( INTERPHISIDEDERIV_C, interphisidederiv_c)

inline void FORTRAN_NAME(inlineINTERPHISIDEDERIV_C, inlineINTERPHISIDEDERIV_C)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_INTERPHISIDEDERIV_C");
 FORTRAN_NAME( INTERPHISIDEDERIV_C ,interphisidederiv_c )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDINTERPHISIDEDERIV_C 

#ifndef GUARDINTERPHISIDEDERIV_V 
#define GUARDINTERPHISIDEDERIV_V 
// Prototype for Fortran procedure interphisidederiv_v ...
//
void FORTRAN_NAME( INTERPHISIDEDERIV_V ,interphisidederiv_v )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) );

#define FORT_INTERPHISIDEDERIV_V FORTRAN_NAME( inlineINTERPHISIDEDERIV_V, inlineINTERPHISIDEDERIV_V)
#define FORTNT_INTERPHISIDEDERIV_V FORTRAN_NAME( INTERPHISIDEDERIV_V, interphisidederiv_v)

inline void FORTRAN_NAME(inlineINTERPHISIDEDERIV_V, inlineINTERPHISIDEDERIV_V)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) )
{
 CH_TIMELEAF("FORT_INTERPHISIDEDERIV_V");
 FORTRAN_NAME( INTERPHISIDEDERIV_V ,interphisidederiv_v )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(dxarr) );
}
#endif  // GUARDINTERPHISIDEDERIV_V 

#ifndef GUARDINTERPLOSIDEDERIV_C 
#define GUARDINTERPLOSIDEDERIV_C 
// Prototype for Fortran procedure interplosidederiv_c ...
//
void FORTRAN_NAME( INTERPLOSIDEDERIV_C ,interplosidederiv_c )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) );

#define FORT_INTERPLOSIDEDERIV_C FORTRAN_NAME( inlineINTERPLOSIDEDERIV_C, inlineINTERPLOSIDEDERIV_C)
#define FORTNT_INTERPLOSIDEDERIV_C FORTRAN_NAME( INTERPLOSIDEDERIV_C, interplosidederiv_c)

inline void FORTRAN_NAME(inlineINTERPLOSIDEDERIV_C, inlineINTERPLOSIDEDERIV_C)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_REAL(dx) )
{
 CH_TIMELEAF("FORT_INTERPLOSIDEDERIV_C");
 FORTRAN_NAME( INTERPLOSIDEDERIV_C ,interplosidederiv_c )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_REAL(dx) );
}
#endif  // GUARDINTERPLOSIDEDERIV_C 

#ifndef GUARDINTERPLOSIDEDERIV_V 
#define GUARDINTERPLOSIDEDERIV_V 
// Prototype for Fortran procedure interplosidederiv_v ...
//
void FORTRAN_NAME( INTERPLOSIDEDERIV_V ,interplosidederiv_v )(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) );

#define FORT_INTERPLOSIDEDERIV_V FORTRAN_NAME( inlineINTERPLOSIDEDERIV_V, inlineINTERPLOSIDEDERIV_V)
#define FORTNT_INTERPLOSIDEDERIV_V FORTRAN_NAME( INTERPLOSIDEDERIV_V, interplosidederiv_v)

inline void FORTRAN_NAME(inlineINTERPLOSIDEDERIV_V, inlineINTERPLOSIDEDERIV_V)(
      CHFp_FRA(slope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_FRA1(dxarr) )
{
 CH_TIMELEAF("FORT_INTERPLOSIDEDERIV_V");
 FORTRAN_NAME( INTERPLOSIDEDERIV_V ,interplosidederiv_v )(
      CHFt_FRA(slope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_FRA1(dxarr) );
}
#endif  // GUARDINTERPLOSIDEDERIV_V 

#ifndef GUARDINTERPLIMIT_V 
#define GUARDINTERPLIMIT_V 
// Prototype for Fortran procedure interplimit_v ...
//
void FORTRAN_NAME( INTERPLIMIT_V ,interplimit_v )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_CONST_FRA1(igrid)
      ,CHFp_CONST_FRA1(jgrid)
      ,CHFp_CONST_FRA1(kgrid)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) );

#define FORT_INTERPLIMIT_V FORTRAN_NAME( inlineINTERPLIMIT_V, inlineINTERPLIMIT_V)
#define FORTNT_INTERPLIMIT_V FORTRAN_NAME( INTERPLIMIT_V, interplimit_v)

inline void FORTRAN_NAME(inlineINTERPLIMIT_V, inlineINTERPLIMIT_V)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_CONST_FRA1(igrid)
      ,CHFp_CONST_FRA1(jgrid)
      ,CHFp_CONST_FRA1(kgrid)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) )
{
 CH_TIMELEAF("FORT_INTERPLIMIT_V");
 FORTRAN_NAME( INTERPLIMIT_V ,interplimit_v )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_CONST_FRA1(igrid)
      ,CHFt_CONST_FRA1(jgrid)
      ,CHFt_CONST_FRA1(kgrid)
      ,CHFt_BOX(bcoarse)
      ,CHFt_BOX(bn)
      ,CHFt_BOX(physdomain) );
}
#endif  // GUARDINTERPLIMIT_V 

#ifndef GUARDINTERPLINEAR_V 
#define GUARDINTERPLINEAR_V 
// Prototype for Fortran procedure interplinear_v ...
//
void FORTRAN_NAME( INTERPLINEAR_V ,interplinear_v )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_CONST_FRA1(crseCC)
      ,CHFp_CONST_FRA1(fineCC)
      ,CHFp_CONST_FRA1(volc)
      ,CHFp_CONST_FRA1(volf)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) );

#define FORT_INTERPLINEAR_V FORTRAN_NAME( inlineINTERPLINEAR_V, inlineINTERPLINEAR_V)
#define FORTNT_INTERPLINEAR_V FORTRAN_NAME( INTERPLINEAR_V, interplinear_v)

inline void FORTRAN_NAME(inlineINTERPLINEAR_V, inlineINTERPLINEAR_V)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_CONST_FRA1(crseCC)
      ,CHFp_CONST_FRA1(fineCC)
      ,CHFp_CONST_FRA1(volc)
      ,CHFp_CONST_FRA1(volf)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_BOX(bref) )
{
 CH_TIMELEAF("FORT_INTERPLINEAR_V");
 FORTRAN_NAME( INTERPLINEAR_V ,interplinear_v )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_CONST_FRA1(crseCC)
      ,CHFt_CONST_FRA1(fineCC)
      ,CHFt_CONST_FRA1(volc)
      ,CHFt_CONST_FRA1(volf)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_BOX(bref) );
}
#endif  // GUARDINTERPLINEAR_V 

#ifndef GUARDINTERP_DIVB_0 
#define GUARDINTERP_DIVB_0 
// Prototype for Fortran procedure interp_divB_0 ...
//
void FORTRAN_NAME( INTERP_DIVB_0 ,interp_divb_0 )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INT(iBX) );

#define FORT_INTERP_DIVB_0 FORTRAN_NAME( inlineINTERP_DIVB_0, inlineINTERP_DIVB_0)
#define FORTNT_INTERP_DIVB_0 FORTRAN_NAME( INTERP_DIVB_0, interp_divb_0)

inline void FORTRAN_NAME(inlineINTERP_DIVB_0, inlineINTERP_DIVB_0)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_BOX(bcoarse)
      ,CHFp_CONST_INT(iBX) )
{
 CH_TIMELEAF("FORT_INTERP_DIVB_0");
 FORTRAN_NAME( INTERP_DIVB_0 ,interp_divb_0 )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_BOX(bcoarse)
      ,CHFt_CONST_INT(iBX) );
}
#endif  // GUARDINTERP_DIVB_0 

#ifndef GUARDCORRECT_DIVB 
#define GUARDCORRECT_DIVB 
// Prototype for Fortran procedure correct_divB ...
//
void FORTRAN_NAME( CORRECT_DIVB ,correct_divb )(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_INT(iBX) );

#define FORT_CORRECT_DIVB FORTRAN_NAME( inlineCORRECT_DIVB, inlineCORRECT_DIVB)
#define FORTNT_CORRECT_DIVB FORTRAN_NAME( CORRECT_DIVB, correct_divb)

inline void FORTRAN_NAME(inlineCORRECT_DIVB, inlineCORRECT_DIVB)(
      CHFp_FRA(fine)
      ,CHFp_CONST_FRA(slope)
      ,CHFp_BOX(b)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(ref_ratio)
      ,CHFp_CONST_INT(iBX) )
{
 CH_TIMELEAF("FORT_CORRECT_DIVB");
 FORTRAN_NAME( CORRECT_DIVB ,correct_divb )(
      CHFt_FRA(fine)
      ,CHFt_CONST_FRA(slope)
      ,CHFt_BOX(b)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(ref_ratio)
      ,CHFt_CONST_INT(iBX) );
}
#endif  // GUARDCORRECT_DIVB 

}

#endif
