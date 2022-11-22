#ifndef _HELIOREALBCPROBLEMF_F_H_
#define _HELIOREALBCPROBLEMF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDODSTRCILPREPARE 
#define GUARDODSTRCILPREPARE 
// Prototype for Fortran procedure ODSTRCILPREPARE ...
//
void FORTRAN_NAME( ODSTRCILPREPARE ,odstrcilprepare )(
      CHFp_FRA(OData)
      ,CHFp_BOX(box) );

#define FORT_ODSTRCILPREPARE FORTRAN_NAME( inlineODSTRCILPREPARE, inlineODSTRCILPREPARE)
#define FORTNT_ODSTRCILPREPARE FORTRAN_NAME( ODSTRCILPREPARE, odstrcilprepare)

inline void FORTRAN_NAME(inlineODSTRCILPREPARE, inlineODSTRCILPREPARE)(
      CHFp_FRA(OData)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ODSTRCILPREPARE");
 FORTRAN_NAME( ODSTRCILPREPARE ,odstrcilprepare )(
      CHFt_FRA(OData)
      ,CHFt_BOX(box) );
}
#endif  // GUARDODSTRCILPREPARE 

#ifndef GUARDODSTRCILINIT 
#define GUARDODSTRCILINIT 
// Prototype for Fortran procedure ODSTRCILINIT ...
//
void FORTRAN_NAME( ODSTRCILINIT ,odstrcilinit )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_FRA(OData)
      ,CHFp_BOX(box) );

#define FORT_ODSTRCILINIT FORTRAN_NAME( inlineODSTRCILINIT, inlineODSTRCILINIT)
#define FORTNT_ODSTRCILINIT FORTRAN_NAME( ODSTRCILINIT, odstrcilinit)

inline void FORTRAN_NAME(inlineODSTRCILINIT, inlineODSTRCILINIT)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_FRA(OData)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ODSTRCILINIT");
 FORTRAN_NAME( ODSTRCILINIT ,odstrcilinit )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_FRA(OData)
      ,CHFt_BOX(box) );
}
#endif  // GUARDODSTRCILINIT 

#ifndef GUARDODSTRCILGSSPHERICAL 
#define GUARDODSTRCILGSSPHERICAL 
// Prototype for Fortran procedure ODSTRCILGSSPHERICAL ...
//
void FORTRAN_NAME( ODSTRCILGSSPHERICAL ,odstrcilgsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_FRA(OData)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_ODSTRCILGSSPHERICAL FORTRAN_NAME( inlineODSTRCILGSSPHERICAL, inlineODSTRCILGSSPHERICAL)
#define FORTNT_ODSTRCILGSSPHERICAL FORTRAN_NAME( ODSTRCILGSSPHERICAL, odstrcilgsspherical)

inline void FORTRAN_NAME(inlineODSTRCILGSSPHERICAL, inlineODSTRCILGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_FRA(OData)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ODSTRCILGSSPHERICAL");
 FORTRAN_NAME( ODSTRCILGSSPHERICAL ,odstrcilgsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_FRA(OData)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDODSTRCILGSSPHERICAL 

#ifndef GUARDHELIOGMIRINITSPHERICAL 
#define GUARDHELIOGMIRINITSPHERICAL 
// Prototype for Fortran procedure HELIOGMIRINITSPHERICAL ...
//
void FORTRAN_NAME( HELIOGMIRINITSPHERICAL ,heliogmirinitspherical )(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_HELIOGMIRINITSPHERICAL FORTRAN_NAME( inlineHELIOGMIRINITSPHERICAL, inlineHELIOGMIRINITSPHERICAL)
#define FORTNT_HELIOGMIRINITSPHERICAL FORTRAN_NAME( HELIOGMIRINITSPHERICAL, heliogmirinitspherical)

inline void FORTRAN_NAME(inlineHELIOGMIRINITSPHERICAL, inlineHELIOGMIRINITSPHERICAL)(
      CHFp_FRA(U)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HELIOGMIRINITSPHERICAL");
 FORTRAN_NAME( HELIOGMIRINITSPHERICAL ,heliogmirinitspherical )(
      CHFt_FRA(U)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHELIOGMIRINITSPHERICAL 

#ifndef GUARDGMIRGSSPHERICAL 
#define GUARDGMIRGSSPHERICAL 
// Prototype for Fortran procedure GMIRGSSPHERICAL ...
//
void FORTRAN_NAME( GMIRGSSPHERICAL ,gmirgsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(startGMIR)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_N,ia_Nhi0)
      ,CHFp_CONST_R1D(a_V,ia_Vhi0)
      ,CHFp_CONST_R1D(a_T,ia_Thi0)
      ,CHFp_CONST_R1D(a_B,ia_Bhi0)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) );

#define FORT_GMIRGSSPHERICAL FORTRAN_NAME( inlineGMIRGSSPHERICAL, inlineGMIRGSSPHERICAL)
#define FORTNT_GMIRGSSPHERICAL FORTRAN_NAME( GMIRGSSPHERICAL, gmirgsspherical)

inline void FORTRAN_NAME(inlineGMIRGSSPHERICAL, inlineGMIRGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(U)
      ,CHFp_CONST_REAL(startGMIR)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_N,ia_Nhi0)
      ,CHFp_CONST_R1D(a_V,ia_Vhi0)
      ,CHFp_CONST_R1D(a_T,ia_Thi0)
      ,CHFp_CONST_R1D(a_B,ia_Bhi0)
      ,CHFp_CONST_INT(iRhoN)
      ,CHFp_CONST_INT(fluids)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_GMIRGSSPHERICAL");
 FORTRAN_NAME( GMIRGSSPHERICAL ,gmirgsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(U)
      ,CHFt_CONST_REAL(startGMIR)
      ,CHFt_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFt_CONST_R1D(a_N,ia_Nhi0)
      ,CHFt_CONST_R1D(a_V,ia_Vhi0)
      ,CHFt_CONST_R1D(a_T,ia_Thi0)
      ,CHFt_CONST_R1D(a_B,ia_Bhi0)
      ,CHFt_CONST_INT(iRhoN)
      ,CHFt_CONST_INT(fluids)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(ltime)
      ,CHFt_BOX(box) );
}
#endif  // GUARDGMIRGSSPHERICAL 

#ifndef GUARDOMNIGSSPHERICAL 
#define GUARDOMNIGSSPHERICAL 
// Prototype for Fortran procedure OMNIGSSPHERICAL ...
//
void FORTRAN_NAME( OMNIGSSPHERICAL ,omnigsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_N,ia_Nhi0)
      ,CHFp_CONST_R1D(a_V,ia_Vhi0)
      ,CHFp_CONST_R1D(a_T,ia_Thi0)
      ,CHFp_CONST_R1D(a_Br,ia_Brhi0)
      ,CHFp_CONST_R1D(a_Bp,ia_Bphi0)
      ,CHFp_CONST_R1D(a_Bt,ia_Bthi0)
      ,CHFp_CONST_R1D(a_B,ia_Bhi0)
      ,CHFp_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iGMIR)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) );

#define FORT_OMNIGSSPHERICAL FORTRAN_NAME( inlineOMNIGSSPHERICAL, inlineOMNIGSSPHERICAL)
#define FORTNT_OMNIGSSPHERICAL FORTRAN_NAME( OMNIGSSPHERICAL, omnigsspherical)

inline void FORTRAN_NAME(inlineOMNIGSSPHERICAL, inlineOMNIGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_N,ia_Nhi0)
      ,CHFp_CONST_R1D(a_V,ia_Vhi0)
      ,CHFp_CONST_R1D(a_T,ia_Thi0)
      ,CHFp_CONST_R1D(a_Br,ia_Brhi0)
      ,CHFp_CONST_R1D(a_Bp,ia_Bphi0)
      ,CHFp_CONST_R1D(a_Bt,ia_Bthi0)
      ,CHFp_CONST_R1D(a_B,ia_Bhi0)
      ,CHFp_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iGMIR)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_OMNIGSSPHERICAL");
 FORTRAN_NAME( OMNIGSSPHERICAL ,omnigsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_REAL(startBC)
      ,CHFt_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFt_CONST_R1D(a_N,ia_Nhi0)
      ,CHFt_CONST_R1D(a_V,ia_Vhi0)
      ,CHFt_CONST_R1D(a_T,ia_Thi0)
      ,CHFt_CONST_R1D(a_Br,ia_Brhi0)
      ,CHFt_CONST_R1D(a_Bp,ia_Bphi0)
      ,CHFt_CONST_R1D(a_Bt,ia_Bthi0)
      ,CHFt_CONST_R1D(a_B,ia_Bhi0)
      ,CHFt_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iGMIR)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(ltime)
      ,CHFt_BOX(box) );
}
#endif  // GUARDOMNIGSSPHERICAL 

#ifndef GUARDWSOGSSPHERICAL 
#define GUARDWSOGSSPHERICAL 
// Prototype for Fortran procedure WSOGSSPHERICAL ...
//
void FORTRAN_NAME( WSOGSSPHERICAL ,wsogsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_CONST_REAL(startBCPhys)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) );

#define FORT_WSOGSSPHERICAL FORTRAN_NAME( inlineWSOGSSPHERICAL, inlineWSOGSSPHERICAL)
#define FORTNT_WSOGSSPHERICAL FORTRAN_NAME( WSOGSSPHERICAL, wsogsspherical)

inline void FORTRAN_NAME(inlineWSOGSSPHERICAL, inlineWSOGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_REAL(startBC)
      ,CHFp_CONST_REAL(startBCPhys)
      ,CHFp_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFp_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(ltime)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_WSOGSSPHERICAL");
 FORTRAN_NAME( WSOGSSPHERICAL ,wsogsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_REAL(startBC)
      ,CHFt_CONST_REAL(startBCPhys)
      ,CHFt_CONST_R1D(a_Time,ia_Timehi0)
      ,CHFt_CONST_R1D(a_Tilt,ia_Tilthi0)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(ltime)
      ,CHFt_BOX(box) );
}
#endif  // GUARDWSOGSSPHERICAL 

#ifndef GUARDTHETAINTERPOLATEDATA 
#define GUARDTHETAINTERPOLATEDATA 
// Prototype for Fortran procedure thetainterpolatedata ...
//
void FORTRAN_NAME( THETAINTERPOLATEDATA ,thetainterpolatedata )(
      CHFp_R1D(tdomain,itdomainhi0)
      ,CHFp_CONST_R1D(tread,itreadhi0)
      ,CHFp_CONST_FRA1(thetac) );

#define FORT_THETAINTERPOLATEDATA FORTRAN_NAME( inlineTHETAINTERPOLATEDATA, inlineTHETAINTERPOLATEDATA)
#define FORTNT_THETAINTERPOLATEDATA FORTRAN_NAME( THETAINTERPOLATEDATA, thetainterpolatedata)

inline void FORTRAN_NAME(inlineTHETAINTERPOLATEDATA, inlineTHETAINTERPOLATEDATA)(
      CHFp_R1D(tdomain,itdomainhi0)
      ,CHFp_CONST_R1D(tread,itreadhi0)
      ,CHFp_CONST_FRA1(thetac) )
{
 CH_TIMELEAF("FORT_THETAINTERPOLATEDATA");
 FORTRAN_NAME( THETAINTERPOLATEDATA ,thetainterpolatedata )(
      CHFt_R1D(tdomain,itdomainhi0)
      ,CHFt_CONST_R1D(tread,itreadhi0)
      ,CHFt_CONST_FRA1(thetac) );
}
#endif  // GUARDTHETAINTERPOLATEDATA 

#ifndef GUARDSUESSTHETABC 
#define GUARDSUESSTHETABC 
// Prototype for Fortran procedure SUESSTHETABC ...
//
void FORTRAN_NAME( SUESSTHETABC ,suessthetabc )(
      CHFp_FRA(W)
      ,CHFp_CONST_R1D(N1,iN1hi0)
      ,CHFp_CONST_R1D(V1,iV1hi0)
      ,CHFp_CONST_R1D(T1,iT1hi0)
      ,CHFp_CONST_REAL(time1)
      ,CHFp_CONST_R1D(N2,iN2hi0)
      ,CHFp_CONST_R1D(V2,iV2hi0)
      ,CHFp_CONST_R1D(T2,iT2hi0)
      ,CHFp_CONST_REAL(time2)
      ,CHFp_CONST_REAL(timeCycle)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(tilt)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) );

#define FORT_SUESSTHETABC FORTRAN_NAME( inlineSUESSTHETABC, inlineSUESSTHETABC)
#define FORTNT_SUESSTHETABC FORTRAN_NAME( SUESSTHETABC, suessthetabc)

inline void FORTRAN_NAME(inlineSUESSTHETABC, inlineSUESSTHETABC)(
      CHFp_FRA(W)
      ,CHFp_CONST_R1D(N1,iN1hi0)
      ,CHFp_CONST_R1D(V1,iV1hi0)
      ,CHFp_CONST_R1D(T1,iT1hi0)
      ,CHFp_CONST_REAL(time1)
      ,CHFp_CONST_R1D(N2,iN2hi0)
      ,CHFp_CONST_R1D(V2,iV2hi0)
      ,CHFp_CONST_R1D(T2,iT2hi0)
      ,CHFp_CONST_REAL(time2)
      ,CHFp_CONST_REAL(timeCycle)
      ,CHFp_CONST_REAL(time)
      ,CHFp_CONST_REAL(tilt)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_SUESSTHETABC");
 FORTRAN_NAME( SUESSTHETABC ,suessthetabc )(
      CHFt_FRA(W)
      ,CHFt_CONST_R1D(N1,iN1hi0)
      ,CHFt_CONST_R1D(V1,iV1hi0)
      ,CHFt_CONST_R1D(T1,iT1hi0)
      ,CHFt_CONST_REAL(time1)
      ,CHFt_CONST_R1D(N2,iN2hi0)
      ,CHFt_CONST_R1D(V2,iV2hi0)
      ,CHFt_CONST_R1D(T2,iT2hi0)
      ,CHFt_CONST_REAL(time2)
      ,CHFt_CONST_REAL(timeCycle)
      ,CHFt_CONST_REAL(time)
      ,CHFt_CONST_REAL(tilt)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_BOX(box) );
}
#endif  // GUARDSUESSTHETABC 

#ifndef GUARDOMNIINTERPOLATETILT 
#define GUARDOMNIINTERPOLATETILT 
// Prototype for Fortran procedure omniinterpolatetilt ...
//
void FORTRAN_NAME( OMNIINTERPOLATETILT ,omniinterpolatetilt )(
      CHFp_R1D(OMNITilt,iOMNITilthi0)
      ,CHFp_CONST_R1D(OMNITime,iOMNITimehi0)
      ,CHFp_CONST_R1D(WSOTilt,iWSOTilthi0)
      ,CHFp_CONST_R1D(WSOTime,iWSOTimehi0) );

#define FORT_OMNIINTERPOLATETILT FORTRAN_NAME( inlineOMNIINTERPOLATETILT, inlineOMNIINTERPOLATETILT)
#define FORTNT_OMNIINTERPOLATETILT FORTRAN_NAME( OMNIINTERPOLATETILT, omniinterpolatetilt)

inline void FORTRAN_NAME(inlineOMNIINTERPOLATETILT, inlineOMNIINTERPOLATETILT)(
      CHFp_R1D(OMNITilt,iOMNITilthi0)
      ,CHFp_CONST_R1D(OMNITime,iOMNITimehi0)
      ,CHFp_CONST_R1D(WSOTilt,iWSOTilthi0)
      ,CHFp_CONST_R1D(WSOTime,iWSOTimehi0) )
{
 CH_TIMELEAF("FORT_OMNIINTERPOLATETILT");
 FORTRAN_NAME( OMNIINTERPOLATETILT ,omniinterpolatetilt )(
      CHFt_R1D(OMNITilt,iOMNITilthi0)
      ,CHFt_CONST_R1D(OMNITime,iOMNITimehi0)
      ,CHFt_CONST_R1D(WSOTilt,iWSOTilthi0)
      ,CHFt_CONST_R1D(WSOTime,iWSOTimehi0) );
}
#endif  // GUARDOMNIINTERPOLATETILT 

#ifndef GUARDCIRGSSPHERICAL 
#define GUARDCIRGSSPHERICAL 
// Prototype for Fortran procedure CIRGSSPHERICAL ...
//
void FORTRAN_NAME( CIRGSSPHERICAL ,cirgsspherical )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(CIRData)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_CIRGSSPHERICAL FORTRAN_NAME( inlineCIRGSSPHERICAL, inlineCIRGSSPHERICAL)
#define FORTNT_CIRGSSPHERICAL FORTRAN_NAME( CIRGSSPHERICAL, cirgsspherical)

inline void FORTRAN_NAME(inlineCIRGSSPHERICAL, inlineCIRGSSPHERICAL)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(CIRData)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_CIRGSSPHERICAL");
 FORTRAN_NAME( CIRGSSPHERICAL ,cirgsspherical )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(CIRData)
      ,CHFt_CONST_REAL(time0)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCIRGSSPHERICAL 

#ifndef GUARDHDF5INPUTGCROTATING 
#define GUARDHDF5INPUTGCROTATING 
// Prototype for Fortran procedure HDF5INPUTGCROTATING ...
//
void FORTRAN_NAME( HDF5INPUTGCROTATING ,hdf5inputgcrotating )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(gcdata)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_HDF5INPUTGCROTATING FORTRAN_NAME( inlineHDF5INPUTGCROTATING, inlineHDF5INPUTGCROTATING)
#define FORTNT_HDF5INPUTGCROTATING FORTRAN_NAME( HDF5INPUTGCROTATING, hdf5inputgcrotating)

inline void FORTRAN_NAME(inlineHDF5INPUTGCROTATING, inlineHDF5INPUTGCROTATING)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(gcdata)
      ,CHFp_CONST_REAL(time0)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HDF5INPUTGCROTATING");
 FORTRAN_NAME( HDF5INPUTGCROTATING ,hdf5inputgcrotating )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(gcdata)
      ,CHFt_CONST_REAL(time0)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHDF5INPUTGCROTATING 

#ifndef GUARDHDF5INPUTGC 
#define GUARDHDF5INPUTGC 
// Prototype for Fortran procedure HDF5INPUTGC ...
//
void FORTRAN_NAME( HDF5INPUTGC ,hdf5inputgc )(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(gcdata)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) );

#define FORT_HDF5INPUTGC FORTRAN_NAME( inlineHDF5INPUTGC, inlineHDF5INPUTGC)
#define FORTNT_HDF5INPUTGC FORTRAN_NAME( HDF5INPUTGC, hdf5inputgc)

inline void FORTRAN_NAME(inlineHDF5INPUTGC, inlineHDF5INPUTGC)(
      CHFp_FRA(W)
      ,CHFp_CONST_FRA(gcdata)
      ,CHFp_CONST_REAL(R0)
      ,CHFp_CONST_INT(iHCS)
      ,CHFp_CONST_INT(iHCSb)
      ,CHFp_CONST_INT(level)
      ,CHFp_CONST_REAL(time)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_HDF5INPUTGC");
 FORTRAN_NAME( HDF5INPUTGC ,hdf5inputgc )(
      CHFt_FRA(W)
      ,CHFt_CONST_FRA(gcdata)
      ,CHFt_CONST_REAL(R0)
      ,CHFt_CONST_INT(iHCS)
      ,CHFt_CONST_INT(iHCSb)
      ,CHFt_CONST_INT(level)
      ,CHFt_CONST_REAL(time)
      ,CHFt_BOX(box) );
}
#endif  // GUARDHDF5INPUTGC 

}

#endif
