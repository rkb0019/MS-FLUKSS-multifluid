#ifndef _RECONSTRUCTIONF_F_H_
#define _RECONSTRUCTIONF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDCONSTANTFACEVALUES 
#define GUARDCONSTANTFACEVALUES 
// Prototype for Fortran procedure CONSTANTFACEVALUES ...
//
void FORTRAN_NAME( CONSTANTFACEVALUES ,constantfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CONSTANTFACEVALUES FORTRAN_NAME( inlineCONSTANTFACEVALUES, inlineCONSTANTFACEVALUES)
#define FORTNT_CONSTANTFACEVALUES FORTRAN_NAME( CONSTANTFACEVALUES, constantfacevalues)

inline void FORTRAN_NAME(inlineCONSTANTFACEVALUES, inlineCONSTANTFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CONSTANTFACEVALUES");
 FORTRAN_NAME( CONSTANTFACEVALUES ,constantfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCONSTANTFACEVALUES 

#ifndef GUARDCENTEREDFACEVALUES 
#define GUARDCENTEREDFACEVALUES 
// Prototype for Fortran procedure CENTEREDFACEVALUES ...
//
void FORTRAN_NAME( CENTEREDFACEVALUES ,centeredfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CENTEREDFACEVALUES FORTRAN_NAME( inlineCENTEREDFACEVALUES, inlineCENTEREDFACEVALUES)
#define FORTNT_CENTEREDFACEVALUES FORTRAN_NAME( CENTEREDFACEVALUES, centeredfacevalues)

inline void FORTRAN_NAME(inlineCENTEREDFACEVALUES, inlineCENTEREDFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CENTEREDFACEVALUES");
 FORTRAN_NAME( CENTEREDFACEVALUES ,centeredfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCENTEREDFACEVALUES 

#ifndef GUARDINITIALSLOPES 
#define GUARDINITIALSLOPES 
// Prototype for Fortran procedure INITIALSLOPES ...
//
void FORTRAN_NAME( INITIALSLOPES ,initialslopes )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(Wleft)
      ,CHFp_FRA(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) );

#define FORT_INITIALSLOPES FORTRAN_NAME( inlineINITIALSLOPES, inlineINITIALSLOPES)
#define FORTNT_INITIALSLOPES FORTRAN_NAME( INITIALSLOPES, initialslopes)

inline void FORTRAN_NAME(inlineINITIALSLOPES, inlineINITIALSLOPES)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(Wleft)
      ,CHFp_FRA(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_INITIALSLOPES");
 FORTRAN_NAME( INITIALSLOPES ,initialslopes )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(Wleft)
      ,CHFt_FRA(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_CONST_INT(iBGN)
      ,CHFt_CONST_INT(iEND)
      ,CHFt_BOX(box) );
}
#endif  // GUARDINITIALSLOPES 

#ifndef GUARDFROMSLOPESTOVALUES 
#define GUARDFROMSLOPESTOVALUES 
// Prototype for Fortran procedure FROMSLOPESTOVALUES ...
//
void FORTRAN_NAME( FROMSLOPESTOVALUES ,fromslopestovalues )(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(Wleft)
      ,CHFp_FRA(Wright)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) );

#define FORT_FROMSLOPESTOVALUES FORTRAN_NAME( inlineFROMSLOPESTOVALUES, inlineFROMSLOPESTOVALUES)
#define FORTNT_FROMSLOPESTOVALUES FORTRAN_NAME( FROMSLOPESTOVALUES, fromslopestovalues)

inline void FORTRAN_NAME(inlineFROMSLOPESTOVALUES, inlineFROMSLOPESTOVALUES)(
      CHFp_CONST_FRA(W)
      ,CHFp_FRA(Wleft)
      ,CHFp_FRA(Wright)
      ,CHFp_CONST_INT(iBGN)
      ,CHFp_CONST_INT(iEND)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_FROMSLOPESTOVALUES");
 FORTRAN_NAME( FROMSLOPESTOVALUES ,fromslopestovalues )(
      CHFt_CONST_FRA(W)
      ,CHFt_FRA(Wleft)
      ,CHFt_FRA(Wright)
      ,CHFt_CONST_INT(iBGN)
      ,CHFt_CONST_INT(iEND)
      ,CHFt_BOX(box) );
}
#endif  // GUARDFROMSLOPESTOVALUES 

#ifndef GUARDMINMODFACEVALUES 
#define GUARDMINMODFACEVALUES 
// Prototype for Fortran procedure MINMODFACEVALUES ...
//
void FORTRAN_NAME( MINMODFACEVALUES ,minmodfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MINMODFACEVALUES FORTRAN_NAME( inlineMINMODFACEVALUES, inlineMINMODFACEVALUES)
#define FORTNT_MINMODFACEVALUES FORTRAN_NAME( MINMODFACEVALUES, minmodfacevalues)

inline void FORTRAN_NAME(inlineMINMODFACEVALUES, inlineMINMODFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINMODFACEVALUES");
 FORTRAN_NAME( MINMODFACEVALUES ,minmodfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINMODFACEVALUES 

#ifndef GUARDMINMODFACESLOPES 
#define GUARDMINMODFACESLOPES 
// Prototype for Fortran procedure MINMODFACESLOPES ...
//
void FORTRAN_NAME( MINMODFACESLOPES ,minmodfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MINMODFACESLOPES FORTRAN_NAME( inlineMINMODFACESLOPES, inlineMINMODFACESLOPES)
#define FORTNT_MINMODFACESLOPES FORTRAN_NAME( MINMODFACESLOPES, minmodfaceslopes)

inline void FORTRAN_NAME(inlineMINMODFACESLOPES, inlineMINMODFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINMODFACESLOPES");
 FORTRAN_NAME( MINMODFACESLOPES ,minmodfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINMODFACESLOPES 

#ifndef GUARDMINMODFACEVALUES_V 
#define GUARDMINMODFACEVALUES_V 
// Prototype for Fortran procedure MINMODFACEVALUES_V ...
//
void FORTRAN_NAME( MINMODFACEVALUES_V ,minmodfacevalues_v )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MINMODFACEVALUES_V FORTRAN_NAME( inlineMINMODFACEVALUES_V, inlineMINMODFACEVALUES_V)
#define FORTNT_MINMODFACEVALUES_V FORTRAN_NAME( MINMODFACEVALUES_V, minmodfacevalues_v)

inline void FORTRAN_NAME(inlineMINMODFACEVALUES_V, inlineMINMODFACEVALUES_V)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINMODFACEVALUES_V");
 FORTRAN_NAME( MINMODFACEVALUES_V ,minmodfacevalues_v )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINMODFACEVALUES_V 

#ifndef GUARDMINMODFACESLOPES_V 
#define GUARDMINMODFACESLOPES_V 
// Prototype for Fortran procedure MINMODFACESLOPES_V ...
//
void FORTRAN_NAME( MINMODFACESLOPES_V ,minmodfaceslopes_v )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MINMODFACESLOPES_V FORTRAN_NAME( inlineMINMODFACESLOPES_V, inlineMINMODFACESLOPES_V)
#define FORTNT_MINMODFACESLOPES_V FORTRAN_NAME( MINMODFACESLOPES_V, minmodfaceslopes_v)

inline void FORTRAN_NAME(inlineMINMODFACESLOPES_V, inlineMINMODFACESLOPES_V)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MINMODFACESLOPES_V");
 FORTRAN_NAME( MINMODFACESLOPES_V ,minmodfaceslopes_v )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMINMODFACESLOPES_V 

#ifndef GUARDSUPERBEEFACEVALUES 
#define GUARDSUPERBEEFACEVALUES 
// Prototype for Fortran procedure SUPERBEEFACEVALUES ...
//
void FORTRAN_NAME( SUPERBEEFACEVALUES ,superbeefacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_SUPERBEEFACEVALUES FORTRAN_NAME( inlineSUPERBEEFACEVALUES, inlineSUPERBEEFACEVALUES)
#define FORTNT_SUPERBEEFACEVALUES FORTRAN_NAME( SUPERBEEFACEVALUES, superbeefacevalues)

inline void FORTRAN_NAME(inlineSUPERBEEFACEVALUES, inlineSUPERBEEFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUPERBEEFACEVALUES");
 FORTRAN_NAME( SUPERBEEFACEVALUES ,superbeefacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUPERBEEFACEVALUES 

#ifndef GUARDSUPERBEEFACEVALUES_V 
#define GUARDSUPERBEEFACEVALUES_V 
// Prototype for Fortran procedure SUPERBEEFACEVALUES_V ...
//
void FORTRAN_NAME( SUPERBEEFACEVALUES_V ,superbeefacevalues_v )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_SUPERBEEFACEVALUES_V FORTRAN_NAME( inlineSUPERBEEFACEVALUES_V, inlineSUPERBEEFACEVALUES_V)
#define FORTNT_SUPERBEEFACEVALUES_V FORTRAN_NAME( SUPERBEEFACEVALUES_V, superbeefacevalues_v)

inline void FORTRAN_NAME(inlineSUPERBEEFACEVALUES_V, inlineSUPERBEEFACEVALUES_V)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUPERBEEFACEVALUES_V");
 FORTRAN_NAME( SUPERBEEFACEVALUES_V ,superbeefacevalues_v )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUPERBEEFACEVALUES_V 

#ifndef GUARDSUPERBEEFACESLOPES 
#define GUARDSUPERBEEFACESLOPES 
// Prototype for Fortran procedure SUPERBEEFACESLOPES ...
//
void FORTRAN_NAME( SUPERBEEFACESLOPES ,superbeefaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_SUPERBEEFACESLOPES FORTRAN_NAME( inlineSUPERBEEFACESLOPES, inlineSUPERBEEFACESLOPES)
#define FORTNT_SUPERBEEFACESLOPES FORTRAN_NAME( SUPERBEEFACESLOPES, superbeefaceslopes)

inline void FORTRAN_NAME(inlineSUPERBEEFACESLOPES, inlineSUPERBEEFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUPERBEEFACESLOPES");
 FORTRAN_NAME( SUPERBEEFACESLOPES ,superbeefaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUPERBEEFACESLOPES 

#ifndef GUARDSUPERBEEFACESLOPES_V 
#define GUARDSUPERBEEFACESLOPES_V 
// Prototype for Fortran procedure SUPERBEEFACESLOPES_V ...
//
void FORTRAN_NAME( SUPERBEEFACESLOPES_V ,superbeefaceslopes_v )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_SUPERBEEFACESLOPES_V FORTRAN_NAME( inlineSUPERBEEFACESLOPES_V, inlineSUPERBEEFACESLOPES_V)
#define FORTNT_SUPERBEEFACESLOPES_V FORTRAN_NAME( SUPERBEEFACESLOPES_V, superbeefaceslopes_v)

inline void FORTRAN_NAME(inlineSUPERBEEFACESLOPES_V, inlineSUPERBEEFACESLOPES_V)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUPERBEEFACESLOPES_V");
 FORTRAN_NAME( SUPERBEEFACESLOPES_V ,superbeefaceslopes_v )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUPERBEEFACESLOPES_V 

#ifndef GUARDHYPERBEEFACEVALUES 
#define GUARDHYPERBEEFACEVALUES 
// Prototype for Fortran procedure HYPERBEEFACEVALUES ...
//
void FORTRAN_NAME( HYPERBEEFACEVALUES ,hyperbeefacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_HYPERBEEFACEVALUES FORTRAN_NAME( inlineHYPERBEEFACEVALUES, inlineHYPERBEEFACEVALUES)
#define FORTNT_HYPERBEEFACEVALUES FORTRAN_NAME( HYPERBEEFACEVALUES, hyperbeefacevalues)

inline void FORTRAN_NAME(inlineHYPERBEEFACEVALUES, inlineHYPERBEEFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HYPERBEEFACEVALUES");
 FORTRAN_NAME( HYPERBEEFACEVALUES ,hyperbeefacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHYPERBEEFACEVALUES 

#ifndef GUARDHYPERBEEFACESLOPES 
#define GUARDHYPERBEEFACESLOPES 
// Prototype for Fortran procedure HYPERBEEFACESLOPES ...
//
void FORTRAN_NAME( HYPERBEEFACESLOPES ,hyperbeefaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_HYPERBEEFACESLOPES FORTRAN_NAME( inlineHYPERBEEFACESLOPES, inlineHYPERBEEFACESLOPES)
#define FORTNT_HYPERBEEFACESLOPES FORTRAN_NAME( HYPERBEEFACESLOPES, hyperbeefaceslopes)

inline void FORTRAN_NAME(inlineHYPERBEEFACESLOPES, inlineHYPERBEEFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HYPERBEEFACESLOPES");
 FORTRAN_NAME( HYPERBEEFACESLOPES ,hyperbeefaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHYPERBEEFACESLOPES 

#ifndef GUARDHARMONICFACEVALUES 
#define GUARDHARMONICFACEVALUES 
// Prototype for Fortran procedure HARMONICFACEVALUES ...
//
void FORTRAN_NAME( HARMONICFACEVALUES ,harmonicfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_HARMONICFACEVALUES FORTRAN_NAME( inlineHARMONICFACEVALUES, inlineHARMONICFACEVALUES)
#define FORTNT_HARMONICFACEVALUES FORTRAN_NAME( HARMONICFACEVALUES, harmonicfacevalues)

inline void FORTRAN_NAME(inlineHARMONICFACEVALUES, inlineHARMONICFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HARMONICFACEVALUES");
 FORTRAN_NAME( HARMONICFACEVALUES ,harmonicfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHARMONICFACEVALUES 

#ifndef GUARDHARMONICFACESLOPES 
#define GUARDHARMONICFACESLOPES 
// Prototype for Fortran procedure HARMONICFACESLOPES ...
//
void FORTRAN_NAME( HARMONICFACESLOPES ,harmonicfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_HARMONICFACESLOPES FORTRAN_NAME( inlineHARMONICFACESLOPES, inlineHARMONICFACESLOPES)
#define FORTNT_HARMONICFACESLOPES FORTRAN_NAME( HARMONICFACESLOPES, harmonicfaceslopes)

inline void FORTRAN_NAME(inlineHARMONICFACESLOPES, inlineHARMONICFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HARMONICFACESLOPES");
 FORTRAN_NAME( HARMONICFACESLOPES ,harmonicfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHARMONICFACESLOPES 

#ifndef GUARDVANALBADAFACEVALUES 
#define GUARDVANALBADAFACEVALUES 
// Prototype for Fortran procedure VANALBADAFACEVALUES ...
//
void FORTRAN_NAME( VANALBADAFACEVALUES ,vanalbadafacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VANALBADAFACEVALUES FORTRAN_NAME( inlineVANALBADAFACEVALUES, inlineVANALBADAFACEVALUES)
#define FORTNT_VANALBADAFACEVALUES FORTRAN_NAME( VANALBADAFACEVALUES, vanalbadafacevalues)

inline void FORTRAN_NAME(inlineVANALBADAFACEVALUES, inlineVANALBADAFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANALBADAFACEVALUES");
 FORTRAN_NAME( VANALBADAFACEVALUES ,vanalbadafacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANALBADAFACEVALUES 

#ifndef GUARDVANALBADAFACESLOPES 
#define GUARDVANALBADAFACESLOPES 
// Prototype for Fortran procedure VANALBADAFACESLOPES ...
//
void FORTRAN_NAME( VANALBADAFACESLOPES ,vanalbadafaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VANALBADAFACESLOPES FORTRAN_NAME( inlineVANALBADAFACESLOPES, inlineVANALBADAFACESLOPES)
#define FORTNT_VANALBADAFACESLOPES FORTRAN_NAME( VANALBADAFACESLOPES, vanalbadafaceslopes)

inline void FORTRAN_NAME(inlineVANALBADAFACESLOPES, inlineVANALBADAFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANALBADAFACESLOPES");
 FORTRAN_NAME( VANALBADAFACESLOPES ,vanalbadafaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANALBADAFACESLOPES 

#ifndef GUARDMCFACEVALUES 
#define GUARDMCFACEVALUES 
// Prototype for Fortran procedure MCFACEVALUES ...
//
void FORTRAN_NAME( MCFACEVALUES ,mcfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MCFACEVALUES FORTRAN_NAME( inlineMCFACEVALUES, inlineMCFACEVALUES)
#define FORTNT_MCFACEVALUES FORTRAN_NAME( MCFACEVALUES, mcfacevalues)

inline void FORTRAN_NAME(inlineMCFACEVALUES, inlineMCFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MCFACEVALUES");
 FORTRAN_NAME( MCFACEVALUES ,mcfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMCFACEVALUES 

#ifndef GUARDMCFACESLOPES 
#define GUARDMCFACESLOPES 
// Prototype for Fortran procedure MCFACESLOPES ...
//
void FORTRAN_NAME( MCFACESLOPES ,mcfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MCFACESLOPES FORTRAN_NAME( inlineMCFACESLOPES, inlineMCFACESLOPES)
#define FORTNT_MCFACESLOPES FORTRAN_NAME( MCFACESLOPES, mcfaceslopes)

inline void FORTRAN_NAME(inlineMCFACESLOPES, inlineMCFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MCFACESLOPES");
 FORTRAN_NAME( MCFACESLOPES ,mcfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMCFACESLOPES 

#ifndef GUARDMCFACEVALUES_V 
#define GUARDMCFACEVALUES_V 
// Prototype for Fortran procedure MCFACEVALUES_V ...
//
void FORTRAN_NAME( MCFACEVALUES_V ,mcfacevalues_v )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MCFACEVALUES_V FORTRAN_NAME( inlineMCFACEVALUES_V, inlineMCFACEVALUES_V)
#define FORTNT_MCFACEVALUES_V FORTRAN_NAME( MCFACEVALUES_V, mcfacevalues_v)

inline void FORTRAN_NAME(inlineMCFACEVALUES_V, inlineMCFACEVALUES_V)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MCFACEVALUES_V");
 FORTRAN_NAME( MCFACEVALUES_V ,mcfacevalues_v )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMCFACEVALUES_V 

#ifndef GUARDMCFACESLOPES_V 
#define GUARDMCFACESLOPES_V 
// Prototype for Fortran procedure MCFACESLOPES_V ...
//
void FORTRAN_NAME( MCFACESLOPES_V ,mcfaceslopes_v )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_MCFACESLOPES_V FORTRAN_NAME( inlineMCFACESLOPES_V, inlineMCFACESLOPES_V)
#define FORTNT_MCFACESLOPES_V FORTRAN_NAME( MCFACESLOPES_V, mcfaceslopes_v)

inline void FORTRAN_NAME(inlineMCFACESLOPES_V, inlineMCFACESLOPES_V)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_FRA1(dxarr)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_MCFACESLOPES_V");
 FORTRAN_NAME( MCFACESLOPES_V ,mcfaceslopes_v )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_FRA1(dxarr)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDMCFACESLOPES_V 

#ifndef GUARDVANLEERFACEVALUES 
#define GUARDVANLEERFACEVALUES 
// Prototype for Fortran procedure VANLEERFACEVALUES ...
//
void FORTRAN_NAME( VANLEERFACEVALUES ,vanleerfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VANLEERFACEVALUES FORTRAN_NAME( inlineVANLEERFACEVALUES, inlineVANLEERFACEVALUES)
#define FORTNT_VANLEERFACEVALUES FORTRAN_NAME( VANLEERFACEVALUES, vanleerfacevalues)

inline void FORTRAN_NAME(inlineVANLEERFACEVALUES, inlineVANLEERFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANLEERFACEVALUES");
 FORTRAN_NAME( VANLEERFACEVALUES ,vanleerfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANLEERFACEVALUES 

#ifndef GUARDVANLEERFACESLOPES 
#define GUARDVANLEERFACESLOPES 
// Prototype for Fortran procedure VANLEERFACESLOPES ...
//
void FORTRAN_NAME( VANLEERFACESLOPES ,vanleerfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VANLEERFACESLOPES FORTRAN_NAME( inlineVANLEERFACESLOPES, inlineVANLEERFACESLOPES)
#define FORTNT_VANLEERFACESLOPES FORTRAN_NAME( VANLEERFACESLOPES, vanleerfaceslopes)

inline void FORTRAN_NAME(inlineVANLEERFACESLOPES, inlineVANLEERFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VANLEERFACESLOPES");
 FORTRAN_NAME( VANLEERFACESLOPES ,vanleerfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVANLEERFACESLOPES 

#ifndef GUARDKORENFACEVALUES 
#define GUARDKORENFACEVALUES 
// Prototype for Fortran procedure KORENFACEVALUES ...
//
void FORTRAN_NAME( KORENFACEVALUES ,korenfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_KORENFACEVALUES FORTRAN_NAME( inlineKORENFACEVALUES, inlineKORENFACEVALUES)
#define FORTNT_KORENFACEVALUES FORTRAN_NAME( KORENFACEVALUES, korenfacevalues)

inline void FORTRAN_NAME(inlineKORENFACEVALUES, inlineKORENFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_KORENFACEVALUES");
 FORTRAN_NAME( KORENFACEVALUES ,korenfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDKORENFACEVALUES 

#ifndef GUARDKORENFACESLOPES 
#define GUARDKORENFACESLOPES 
// Prototype for Fortran procedure KORENFACESLOPES ...
//
void FORTRAN_NAME( KORENFACESLOPES ,korenfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_KORENFACESLOPES FORTRAN_NAME( inlineKORENFACESLOPES, inlineKORENFACESLOPES)
#define FORTNT_KORENFACESLOPES FORTRAN_NAME( KORENFACESLOPES, korenfaceslopes)

inline void FORTRAN_NAME(inlineKORENFACESLOPES, inlineKORENFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_KORENFACESLOPES");
 FORTRAN_NAME( KORENFACESLOPES ,korenfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDKORENFACESLOPES 

#ifndef GUARDVENKATAKRISHNANFACEVALUES 
#define GUARDVENKATAKRISHNANFACEVALUES 
// Prototype for Fortran procedure VENKATAKRISHNANFACEVALUES ...
//
void FORTRAN_NAME( VENKATAKRISHNANFACEVALUES ,venkatakrishnanfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VENKATAKRISHNANFACEVALUES FORTRAN_NAME( inlineVENKATAKRISHNANFACEVALUES, inlineVENKATAKRISHNANFACEVALUES)
#define FORTNT_VENKATAKRISHNANFACEVALUES FORTRAN_NAME( VENKATAKRISHNANFACEVALUES, venkatakrishnanfacevalues)

inline void FORTRAN_NAME(inlineVENKATAKRISHNANFACEVALUES, inlineVENKATAKRISHNANFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VENKATAKRISHNANFACEVALUES");
 FORTRAN_NAME( VENKATAKRISHNANFACEVALUES ,venkatakrishnanfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVENKATAKRISHNANFACEVALUES 

#ifndef GUARDVENKATAKRISHNANFACESLOPES 
#define GUARDVENKATAKRISHNANFACESLOPES 
// Prototype for Fortran procedure VENKATAKRISHNANFACESLOPES ...
//
void FORTRAN_NAME( VENKATAKRISHNANFACESLOPES ,venkatakrishnanfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_VENKATAKRISHNANFACESLOPES FORTRAN_NAME( inlineVENKATAKRISHNANFACESLOPES, inlineVENKATAKRISHNANFACESLOPES)
#define FORTNT_VENKATAKRISHNANFACESLOPES FORTRAN_NAME( VENKATAKRISHNANFACESLOPES, venkatakrishnanfaceslopes)

inline void FORTRAN_NAME(inlineVENKATAKRISHNANFACESLOPES, inlineVENKATAKRISHNANFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_VENKATAKRISHNANFACESLOPES");
 FORTRAN_NAME( VENKATAKRISHNANFACESLOPES ,venkatakrishnanfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDVENKATAKRISHNANFACESLOPES 

#ifndef GUARDCHECKPOSITIVITY 
#define GUARDCHECKPOSITIVITY 
// Prototype for Fortran procedure CHECKPOSITIVITY ...
//
void FORTRAN_NAME( CHECKPOSITIVITY ,checkpositivity )(
      CHFp_FIA1(neg)
      ,CHFp_CONST_FRA1(WL)
      ,CHFp_CONST_FRA1(WR)
      ,CHFp_BOX(box) );

#define FORT_CHECKPOSITIVITY FORTRAN_NAME( inlineCHECKPOSITIVITY, inlineCHECKPOSITIVITY)
#define FORTNT_CHECKPOSITIVITY FORTRAN_NAME( CHECKPOSITIVITY, checkpositivity)

inline void FORTRAN_NAME(inlineCHECKPOSITIVITY, inlineCHECKPOSITIVITY)(
      CHFp_FIA1(neg)
      ,CHFp_CONST_FRA1(WL)
      ,CHFp_CONST_FRA1(WR)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHECKPOSITIVITY");
 FORTRAN_NAME( CHECKPOSITIVITY ,checkpositivity )(
      CHFt_FIA1(neg)
      ,CHFt_CONST_FRA1(WL)
      ,CHFt_CONST_FRA1(WR)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHECKPOSITIVITY 

#ifndef GUARDREDUCEORDER 
#define GUARDREDUCEORDER 
// Prototype for Fortran procedure REDUCEORDER ...
//
void FORTRAN_NAME( REDUCEORDER ,reduceorder )(
      CHFp_FRA1(WL)
      ,CHFp_FRA1(WR)
      ,CHFp_CONST_FRA1(W)
      ,CHFp_CONST_FIA1(neg)
      ,CHFp_BOX(box) );

#define FORT_REDUCEORDER FORTRAN_NAME( inlineREDUCEORDER, inlineREDUCEORDER)
#define FORTNT_REDUCEORDER FORTRAN_NAME( REDUCEORDER, reduceorder)

inline void FORTRAN_NAME(inlineREDUCEORDER, inlineREDUCEORDER)(
      CHFp_FRA1(WL)
      ,CHFp_FRA1(WR)
      ,CHFp_CONST_FRA1(W)
      ,CHFp_CONST_FIA1(neg)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_REDUCEORDER");
 FORTRAN_NAME( REDUCEORDER ,reduceorder )(
      CHFt_FRA1(WL)
      ,CHFt_FRA1(WR)
      ,CHFt_CONST_FRA1(W)
      ,CHFt_CONST_FIA1(neg)
      ,CHFt_BOX(box) );
}
#endif  // GUARDREDUCEORDER 

#ifndef GUARDCDFACEVALUES 
#define GUARDCDFACEVALUES 
// Prototype for Fortran procedure CDFACEVALUES ...
//
void FORTRAN_NAME( CDFACEVALUES ,cdfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CDFACEVALUES FORTRAN_NAME( inlineCDFACEVALUES, inlineCDFACEVALUES)
#define FORTNT_CDFACEVALUES FORTRAN_NAME( CDFACEVALUES, cdfacevalues)

inline void FORTRAN_NAME(inlineCDFACEVALUES, inlineCDFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CDFACEVALUES");
 FORTRAN_NAME( CDFACEVALUES ,cdfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCDFACEVALUES 

#ifndef GUARDCDFACESLOPES 
#define GUARDCDFACESLOPES 
// Prototype for Fortran procedure CDFACESLOPES ...
//
void FORTRAN_NAME( CDFACESLOPES ,cdfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CDFACESLOPES FORTRAN_NAME( inlineCDFACESLOPES, inlineCDFACESLOPES)
#define FORTNT_CDFACESLOPES FORTRAN_NAME( CDFACESLOPES, cdfaceslopes)

inline void FORTRAN_NAME(inlineCDFACESLOPES, inlineCDFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CDFACESLOPES");
 FORTRAN_NAME( CDFACESLOPES ,cdfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCDFACESLOPES 

#ifndef GUARDWENO3FACEVALUES 
#define GUARDWENO3FACEVALUES 
// Prototype for Fortran procedure WENO3FACEVALUES ...
//
void FORTRAN_NAME( WENO3FACEVALUES ,weno3facevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_WENO3FACEVALUES FORTRAN_NAME( inlineWENO3FACEVALUES, inlineWENO3FACEVALUES)
#define FORTNT_WENO3FACEVALUES FORTRAN_NAME( WENO3FACEVALUES, weno3facevalues)

inline void FORTRAN_NAME(inlineWENO3FACEVALUES, inlineWENO3FACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WENO3FACEVALUES");
 FORTRAN_NAME( WENO3FACEVALUES ,weno3facevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWENO3FACEVALUES 

#ifndef GUARDWENO3FACESLOPES 
#define GUARDWENO3FACESLOPES 
// Prototype for Fortran procedure WENO3FACESLOPES ...
//
void FORTRAN_NAME( WENO3FACESLOPES ,weno3faceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_WENO3FACESLOPES FORTRAN_NAME( inlineWENO3FACESLOPES, inlineWENO3FACESLOPES)
#define FORTNT_WENO3FACESLOPES FORTRAN_NAME( WENO3FACESLOPES, weno3faceslopes)

inline void FORTRAN_NAME(inlineWENO3FACESLOPES, inlineWENO3FACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WENO3FACESLOPES");
 FORTRAN_NAME( WENO3FACESLOPES ,weno3faceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWENO3FACESLOPES 

#ifndef GUARDWENO3YCFACEVALUES 
#define GUARDWENO3YCFACEVALUES 
// Prototype for Fortran procedure WENO3YCFACEVALUES ...
//
void FORTRAN_NAME( WENO3YCFACEVALUES ,weno3ycfacevalues )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_WENO3YCFACEVALUES FORTRAN_NAME( inlineWENO3YCFACEVALUES, inlineWENO3YCFACEVALUES)
#define FORTNT_WENO3YCFACEVALUES FORTRAN_NAME( WENO3YCFACEVALUES, weno3ycfacevalues)

inline void FORTRAN_NAME(inlineWENO3YCFACEVALUES, inlineWENO3YCFACEVALUES)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WENO3YCFACEVALUES");
 FORTRAN_NAME( WENO3YCFACEVALUES ,weno3ycfacevalues )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWENO3YCFACEVALUES 

#ifndef GUARDWENO3YCFACESLOPES 
#define GUARDWENO3YCFACESLOPES 
// Prototype for Fortran procedure WENO3YCFACESLOPES ...
//
void FORTRAN_NAME( WENO3YCFACESLOPES ,weno3ycfaceslopes )(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_WENO3YCFACESLOPES FORTRAN_NAME( inlineWENO3YCFACESLOPES, inlineWENO3YCFACESLOPES)
#define FORTNT_WENO3YCFACESLOPES FORTRAN_NAME( WENO3YCFACESLOPES, weno3ycfaceslopes)

inline void FORTRAN_NAME(inlineWENO3YCFACESLOPES, inlineWENO3YCFACESLOPES)(
      CHFp_FRA1(Sleft)
      ,CHFp_FRA1(Sright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WENO3YCFACESLOPES");
 FORTRAN_NAME( WENO3YCFACESLOPES ,weno3ycfaceslopes )(
      CHFt_FRA1(Sleft)
      ,CHFt_FRA1(Sright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWENO3YCFACESLOPES 

#ifndef GUARDCHECKTVD 
#define GUARDCHECKTVD 
// Prototype for Fortran procedure CHECKTVD ...
//
void FORTRAN_NAME( CHECKTVD ,checktvd )(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) );

#define FORT_CHECKTVD FORTRAN_NAME( inlineCHECKTVD, inlineCHECKTVD)
#define FORTNT_CHECKTVD FORTRAN_NAME( CHECKTVD, checktvd)

inline void FORTRAN_NAME(inlineCHECKTVD, inlineCHECKTVD)(
      CHFp_CONST_FRA1(W)
      ,CHFp_FRA1(Wleft)
      ,CHFp_FRA1(Wright)
      ,CHFp_CONST_INT(idir)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CHECKTVD");
 FORTRAN_NAME( CHECKTVD ,checktvd )(
      CHFt_CONST_FRA1(W)
      ,CHFt_FRA1(Wleft)
      ,CHFt_FRA1(Wright)
      ,CHFt_CONST_INT(idir)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCHECKTVD 

}

#endif
