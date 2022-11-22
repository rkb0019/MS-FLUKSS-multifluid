#ifndef _CSHSPHERICALF_F_H_
#define _CSHSPHERICALF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDSETLEVELGRIDSPACING 
#define GUARDSETLEVELGRIDSPACING 
// Prototype for Fortran procedure SETLEVELGRIDSPACING ...
//
void FORTRAN_NAME( SETLEVELGRIDSPACING ,setlevelgridspacing )(
      CHFp_FRA1(dudvdw)
      ,CHFp_FRA1(uvwc)
      ,CHFp_FRA1(uvwn)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) );

#define FORT_SETLEVELGRIDSPACING FORTRAN_NAME( inlineSETLEVELGRIDSPACING, inlineSETLEVELGRIDSPACING)
#define FORTNT_SETLEVELGRIDSPACING FORTRAN_NAME( SETLEVELGRIDSPACING, setlevelgridspacing)

inline void FORTRAN_NAME(inlineSETLEVELGRIDSPACING, inlineSETLEVELGRIDSPACING)(
      CHFp_FRA1(dudvdw)
      ,CHFp_FRA1(uvwc)
      ,CHFp_FRA1(uvwn)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_SETLEVELGRIDSPACING");
 FORTRAN_NAME( SETLEVELGRIDSPACING ,setlevelgridspacing )(
      CHFt_FRA1(dudvdw)
      ,CHFt_FRA1(uvwc)
      ,CHFt_FRA1(uvwn)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDSETLEVELGRIDSPACING 

#ifndef GUARDVOLUMESPHERICAL 
#define GUARDVOLUMESPHERICAL 
// Prototype for Fortran procedure VOLUMESPHERICAL ...
//
void FORTRAN_NAME( VOLUMESPHERICAL ,volumespherical )(
      CHFp_FRA1(vol)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_VOLUMESPHERICAL FORTRAN_NAME( inlineVOLUMESPHERICAL, inlineVOLUMESPHERICAL)
#define FORTNT_VOLUMESPHERICAL FORTRAN_NAME( VOLUMESPHERICAL, volumespherical)

inline void FORTRAN_NAME(inlineVOLUMESPHERICAL, inlineVOLUMESPHERICAL)(
      CHFp_FRA1(vol)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_VOLUMESPHERICAL");
 FORTRAN_NAME( VOLUMESPHERICAL ,volumespherical )(
      CHFt_FRA1(vol)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDVOLUMESPHERICAL 

#ifndef GUARDGETNODECOORDSSPHERICAL_CS 
#define GUARDGETNODECOORDSSPHERICAL_CS 
// Prototype for Fortran procedure GETNODECOORDSSPHERICAL_CS ...
//
void FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS ,getnodecoordsspherical_cs )(
      CHFp_REALVECT(coords)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_INT(level) );

#define FORT_GETNODECOORDSSPHERICAL_CS FORTRAN_NAME( inlineGETNODECOORDSSPHERICAL_CS, inlineGETNODECOORDSSPHERICAL_CS)
#define FORTNT_GETNODECOORDSSPHERICAL_CS FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS, getnodecoordsspherical_cs)

inline void FORTRAN_NAME(inlineGETNODECOORDSSPHERICAL_CS, inlineGETNODECOORDSSPHERICAL_CS)(
      CHFp_REALVECT(coords)
      ,CHFp_CONST_INTVECT(iv)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_GETNODECOORDSSPHERICAL_CS");
 FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS ,getnodecoordsspherical_cs )(
      CHFt_REALVECT(coords)
      ,CHFt_CONST_INTVECT(iv)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDGETNODECOORDSSPHERICAL_CS 

#ifndef GUARDGETNODECOORDSSPHERICAL_CS_FAB 
#define GUARDGETNODECOORDSSPHERICAL_CS_FAB 
// Prototype for Fortran procedure GETNODECOORDSSPHERICAL_CS_FAB ...
//
void FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS_FAB ,getnodecoordsspherical_cs_fab )(
      CHFp_FRA(coords)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_GETNODECOORDSSPHERICAL_CS_FAB FORTRAN_NAME( inlineGETNODECOORDSSPHERICAL_CS_FAB, inlineGETNODECOORDSSPHERICAL_CS_FAB)
#define FORTNT_GETNODECOORDSSPHERICAL_CS_FAB FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS_FAB, getnodecoordsspherical_cs_fab)

inline void FORTRAN_NAME(inlineGETNODECOORDSSPHERICAL_CS_FAB, inlineGETNODECOORDSSPHERICAL_CS_FAB)(
      CHFp_FRA(coords)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_GETNODECOORDSSPHERICAL_CS_FAB");
 FORTRAN_NAME( GETNODECOORDSSPHERICAL_CS_FAB ,getnodecoordsspherical_cs_fab )(
      CHFt_FRA(coords)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDGETNODECOORDSSPHERICAL_CS_FAB 

#ifndef GUARDAREASPOLAR 
#define GUARDAREASPOLAR 
// Prototype for Fortran procedure AREASPOLAR ...
//
void FORTRAN_NAME( AREASPOLAR ,areaspolar )(
      CHFp_FRA1(area)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) );

#define FORT_AREASPOLAR FORTRAN_NAME( inlineAREASPOLAR, inlineAREASPOLAR)
#define FORTNT_AREASPOLAR FORTRAN_NAME( AREASPOLAR, areaspolar)

inline void FORTRAN_NAME(inlineAREASPOLAR, inlineAREASPOLAR)(
      CHFp_FRA1(area)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_AREASPOLAR");
 FORTRAN_NAME( AREASPOLAR ,areaspolar )(
      CHFt_FRA1(area)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDAREASPOLAR 

#ifndef GUARDAREASSPHERICAL 
#define GUARDAREASSPHERICAL 
// Prototype for Fortran procedure AREASSPHERICAL ...
//
void FORTRAN_NAME( AREASSPHERICAL ,areasspherical )(
      CHFp_FRA1(area)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) );

#define FORT_AREASSPHERICAL FORTRAN_NAME( inlineAREASSPHERICAL, inlineAREASSPHERICAL)
#define FORTNT_AREASSPHERICAL FORTRAN_NAME( AREASSPHERICAL, areasspherical)

inline void FORTRAN_NAME(inlineAREASSPHERICAL, inlineAREASSPHERICAL)(
      CHFp_FRA1(area)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_AREASSPHERICAL");
 FORTRAN_NAME( AREASSPHERICAL ,areasspherical )(
      CHFt_FRA1(area)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDAREASSPHERICAL 

#ifndef GUARDTRANSCARTESIANTOPOLAR 
#define GUARDTRANSCARTESIANTOPOLAR 
// Prototype for Fortran procedure TRANSCARTESIANTOPOLAR ...
//
void FORTRAN_NAME( TRANSCARTESIANTOPOLAR ,transcartesiantopolar )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_TRANSCARTESIANTOPOLAR FORTRAN_NAME( inlineTRANSCARTESIANTOPOLAR, inlineTRANSCARTESIANTOPOLAR)
#define FORTNT_TRANSCARTESIANTOPOLAR FORTRAN_NAME( TRANSCARTESIANTOPOLAR, transcartesiantopolar)

inline void FORTRAN_NAME(inlineTRANSCARTESIANTOPOLAR, inlineTRANSCARTESIANTOPOLAR)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_TRANSCARTESIANTOPOLAR");
 FORTRAN_NAME( TRANSCARTESIANTOPOLAR ,transcartesiantopolar )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDTRANSCARTESIANTOPOLAR 

#ifndef GUARDTRANSCARTESIANTOSPHERICAL 
#define GUARDTRANSCARTESIANTOSPHERICAL 
// Prototype for Fortran procedure TRANSCARTESIANTOSPHERICAL ...
//
void FORTRAN_NAME( TRANSCARTESIANTOSPHERICAL ,transcartesiantospherical )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_TRANSCARTESIANTOSPHERICAL FORTRAN_NAME( inlineTRANSCARTESIANTOSPHERICAL, inlineTRANSCARTESIANTOSPHERICAL)
#define FORTNT_TRANSCARTESIANTOSPHERICAL FORTRAN_NAME( TRANSCARTESIANTOSPHERICAL, transcartesiantospherical)

inline void FORTRAN_NAME(inlineTRANSCARTESIANTOSPHERICAL, inlineTRANSCARTESIANTOSPHERICAL)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_TRANSCARTESIANTOSPHERICAL");
 FORTRAN_NAME( TRANSCARTESIANTOSPHERICAL ,transcartesiantospherical )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDTRANSCARTESIANTOSPHERICAL 

#ifndef GUARDTRANSCARTESIANTOSPHERICALZ 
#define GUARDTRANSCARTESIANTOSPHERICALZ 
// Prototype for Fortran procedure TRANSCARTESIANTOSPHERICALZ ...
//
void FORTRAN_NAME( TRANSCARTESIANTOSPHERICALZ ,transcartesiantosphericalz )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_TRANSCARTESIANTOSPHERICALZ FORTRAN_NAME( inlineTRANSCARTESIANTOSPHERICALZ, inlineTRANSCARTESIANTOSPHERICALZ)
#define FORTNT_TRANSCARTESIANTOSPHERICALZ FORTRAN_NAME( TRANSCARTESIANTOSPHERICALZ, transcartesiantosphericalz)

inline void FORTRAN_NAME(inlineTRANSCARTESIANTOSPHERICALZ, inlineTRANSCARTESIANTOSPHERICALZ)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_TRANSCARTESIANTOSPHERICALZ");
 FORTRAN_NAME( TRANSCARTESIANTOSPHERICALZ ,transcartesiantosphericalz )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDTRANSCARTESIANTOSPHERICALZ 

#ifndef GUARDTRANSPOLARTOCARTESIAN 
#define GUARDTRANSPOLARTOCARTESIAN 
// Prototype for Fortran procedure TRANSPOLARTOCARTESIAN ...
//
void FORTRAN_NAME( TRANSPOLARTOCARTESIAN ,transpolartocartesian )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(iType)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_TRANSPOLARTOCARTESIAN FORTRAN_NAME( inlineTRANSPOLARTOCARTESIAN, inlineTRANSPOLARTOCARTESIAN)
#define FORTNT_TRANSPOLARTOCARTESIAN FORTRAN_NAME( TRANSPOLARTOCARTESIAN, transpolartocartesian)

inline void FORTRAN_NAME(inlineTRANSPOLARTOCARTESIAN, inlineTRANSPOLARTOCARTESIAN)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(iType)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_TRANSPOLARTOCARTESIAN");
 FORTRAN_NAME( TRANSPOLARTOCARTESIAN ,transpolartocartesian )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_CONST_INT(iType)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDTRANSPOLARTOCARTESIAN 

#ifndef GUARDTRANSSPHERICALTOCARTESIAN 
#define GUARDTRANSSPHERICALTOCARTESIAN 
// Prototype for Fortran procedure TRANSSPHERICALTOCARTESIAN ...
//
void FORTRAN_NAME( TRANSSPHERICALTOCARTESIAN ,transsphericaltocartesian )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(iType)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_TRANSSPHERICALTOCARTESIAN FORTRAN_NAME( inlineTRANSSPHERICALTOCARTESIAN, inlineTRANSSPHERICALTOCARTESIAN)
#define FORTNT_TRANSSPHERICALTOCARTESIAN FORTRAN_NAME( TRANSSPHERICALTOCARTESIAN, transsphericaltocartesian)

inline void FORTRAN_NAME(inlineTRANSSPHERICALTOCARTESIAN, inlineTRANSSPHERICALTOCARTESIAN)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(iType)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_TRANSSPHERICALTOCARTESIAN");
 FORTRAN_NAME( TRANSSPHERICALTOCARTESIAN ,transsphericaltocartesian )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_CONST_INT(iType)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDTRANSSPHERICALTOCARTESIAN 

#ifndef GUARDCOMPUTE_VELOCITIES 
#define GUARDCOMPUTE_VELOCITIES 
// Prototype for Fortran procedure COMPUTE_VELOCITIES ...
//
void FORTRAN_NAME( COMPUTE_VELOCITIES ,compute_velocities )(
      CHFp_FRA(U)
      ,CHFp_CONST_I1D(velvec,ivelvechi0)
      ,CHFp_BOX(box) );

#define FORT_COMPUTE_VELOCITIES FORTRAN_NAME( inlineCOMPUTE_VELOCITIES, inlineCOMPUTE_VELOCITIES)
#define FORTNT_COMPUTE_VELOCITIES FORTRAN_NAME( COMPUTE_VELOCITIES, compute_velocities)

inline void FORTRAN_NAME(inlineCOMPUTE_VELOCITIES, inlineCOMPUTE_VELOCITIES)(
      CHFp_FRA(U)
      ,CHFp_CONST_I1D(velvec,ivelvechi0)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPUTE_VELOCITIES");
 FORTRAN_NAME( COMPUTE_VELOCITIES ,compute_velocities )(
      CHFt_FRA(U)
      ,CHFt_CONST_I1D(velvec,ivelvechi0)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPUTE_VELOCITIES 

#ifndef GUARDCOMPUTE_MOMENTUM 
#define GUARDCOMPUTE_MOMENTUM 
// Prototype for Fortran procedure COMPUTE_MOMENTUM ...
//
void FORTRAN_NAME( COMPUTE_MOMENTUM ,compute_momentum )(
      CHFp_FRA(U)
      ,CHFp_CONST_I1D(velvec,ivelvechi0)
      ,CHFp_BOX(box) );

#define FORT_COMPUTE_MOMENTUM FORTRAN_NAME( inlineCOMPUTE_MOMENTUM, inlineCOMPUTE_MOMENTUM)
#define FORTNT_COMPUTE_MOMENTUM FORTRAN_NAME( COMPUTE_MOMENTUM, compute_momentum)

inline void FORTRAN_NAME(inlineCOMPUTE_MOMENTUM, inlineCOMPUTE_MOMENTUM)(
      CHFp_FRA(U)
      ,CHFp_CONST_I1D(velvec,ivelvechi0)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_COMPUTE_MOMENTUM");
 FORTRAN_NAME( COMPUTE_MOMENTUM ,compute_momentum )(
      CHFt_FRA(U)
      ,CHFt_CONST_I1D(velvec,ivelvechi0)
      ,CHFt_BOX(box) );
}
#endif  // GUARDCOMPUTE_MOMENTUM 

#ifndef GUARDFR_INCREMENTFINE_POLAR 
#define GUARDFR_INCREMENTFINE_POLAR 
// Prototype for Fortran procedure FR_INCREMENTFINE_POLAR ...
//
void FORTRAN_NAME( FR_INCREMENTFINE_POLAR ,fr_incrementfine_polar )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) );

#define FORT_FR_INCREMENTFINE_POLAR FORTRAN_NAME( inlineFR_INCREMENTFINE_POLAR, inlineFR_INCREMENTFINE_POLAR)
#define FORTNT_FR_INCREMENTFINE_POLAR FORTRAN_NAME( FR_INCREMENTFINE_POLAR, fr_incrementfine_polar)

inline void FORTRAN_NAME(inlineFR_INCREMENTFINE_POLAR, inlineFR_INCREMENTFINE_POLAR)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_FR_INCREMENTFINE_POLAR");
 FORTRAN_NAME( FR_INCREMENTFINE_POLAR ,fr_incrementfine_polar )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDFR_INCREMENTFINE_POLAR 

#ifndef GUARDFR_INCREMENTFINE_SPHERICAL 
#define GUARDFR_INCREMENTFINE_SPHERICAL 
// Prototype for Fortran procedure FR_INCREMENTFINE_SPHERICAL ...
//
void FORTRAN_NAME( FR_INCREMENTFINE_SPHERICAL ,fr_incrementfine_spherical )(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) );

#define FORT_FR_INCREMENTFINE_SPHERICAL FORTRAN_NAME( inlineFR_INCREMENTFINE_SPHERICAL, inlineFR_INCREMENTFINE_SPHERICAL)
#define FORTNT_FR_INCREMENTFINE_SPHERICAL FORTRAN_NAME( FR_INCREMENTFINE_SPHERICAL, fr_incrementfine_spherical)

inline void FORTRAN_NAME(inlineFR_INCREMENTFINE_SPHERICAL, inlineFR_INCREMENTFINE_SPHERICAL)(
      CHFp_FRA(a_U)
      ,CHFp_CONST_I1D(vectors,ivectorshi0)
      ,CHFp_CONST_INT(refRatio)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(dir)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_FR_INCREMENTFINE_SPHERICAL");
 FORTRAN_NAME( FR_INCREMENTFINE_SPHERICAL ,fr_incrementfine_spherical )(
      CHFt_FRA(a_U)
      ,CHFt_CONST_I1D(vectors,ivectorshi0)
      ,CHFt_CONST_INT(refRatio)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(dir)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDFR_INCREMENTFINE_SPHERICAL 

#ifndef GUARDMETRICTENSOR_R 
#define GUARDMETRICTENSOR_R 
// Prototype for Fortran procedure METRICTENSOR_R ...
//
void FORTRAN_NAME( METRICTENSOR_R ,metrictensor_r )(
      CHFp_FRA1(g)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_METRICTENSOR_R FORTRAN_NAME( inlineMETRICTENSOR_R, inlineMETRICTENSOR_R)
#define FORTNT_METRICTENSOR_R FORTRAN_NAME( METRICTENSOR_R, metrictensor_r)

inline void FORTRAN_NAME(inlineMETRICTENSOR_R, inlineMETRICTENSOR_R)(
      CHFp_FRA1(g)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_METRICTENSOR_R");
 FORTRAN_NAME( METRICTENSOR_R ,metrictensor_r )(
      CHFt_FRA1(g)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDMETRICTENSOR_R 

#ifndef GUARDMETRICTENSOR11_SPHERICAL 
#define GUARDMETRICTENSOR11_SPHERICAL 
// Prototype for Fortran procedure METRICTENSOR11_SPHERICAL ...
//
void FORTRAN_NAME( METRICTENSOR11_SPHERICAL ,metrictensor11_spherical )(
      CHFp_FRA1(g)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) );

#define FORT_METRICTENSOR11_SPHERICAL FORTRAN_NAME( inlineMETRICTENSOR11_SPHERICAL, inlineMETRICTENSOR11_SPHERICAL)
#define FORTNT_METRICTENSOR11_SPHERICAL FORTRAN_NAME( METRICTENSOR11_SPHERICAL, metrictensor11_spherical)

inline void FORTRAN_NAME(inlineMETRICTENSOR11_SPHERICAL, inlineMETRICTENSOR11_SPHERICAL)(
      CHFp_FRA1(g)
      ,CHFp_BOX(box)
      ,CHFp_CONST_INT(level) )
{
 CH_TIMELEAF("FORT_METRICTENSOR11_SPHERICAL");
 FORTRAN_NAME( METRICTENSOR11_SPHERICAL ,metrictensor11_spherical )(
      CHFt_FRA1(g)
      ,CHFt_BOX(box)
      ,CHFt_CONST_INT(level) );
}
#endif  // GUARDMETRICTENSOR11_SPHERICAL 

#ifndef GUARDZAXISCOPY 
#define GUARDZAXISCOPY 
// Prototype for Fortran procedure ZAXISCOPY ...
//
void FORTRAN_NAME( ZAXISCOPY ,zaxiscopy )(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(a_GS)
      ,CHFp_CONST_INT(sign) );

#define FORT_ZAXISCOPY FORTRAN_NAME( inlineZAXISCOPY, inlineZAXISCOPY)
#define FORTNT_ZAXISCOPY FORTRAN_NAME( ZAXISCOPY, zaxiscopy)

inline void FORTRAN_NAME(inlineZAXISCOPY, inlineZAXISCOPY)(
      CHFp_CONST_FRA(a_from)
      ,CHFp_FRA(a_to)
      ,CHFp_CONST_INT(a_GS)
      ,CHFp_CONST_INT(sign) )
{
 CH_TIMELEAF("FORT_ZAXISCOPY");
 FORTRAN_NAME( ZAXISCOPY ,zaxiscopy )(
      CHFt_CONST_FRA(a_from)
      ,CHFt_FRA(a_to)
      ,CHFt_CONST_INT(a_GS)
      ,CHFt_CONST_INT(sign) );
}
#endif  // GUARDZAXISCOPY 

#ifndef GUARDZAXISFIRSTORDER 
#define GUARDZAXISFIRSTORDER 
// Prototype for Fortran procedure ZAXISFIRSTORDER ...
//
void FORTRAN_NAME( ZAXISFIRSTORDER ,zaxisfirstorder )(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(kbnd)
      ,CHFp_CONST_INT(nGS)
      ,CHFp_BOX(box) );

#define FORT_ZAXISFIRSTORDER FORTRAN_NAME( inlineZAXISFIRSTORDER, inlineZAXISFIRSTORDER)
#define FORTNT_ZAXISFIRSTORDER FORTRAN_NAME( ZAXISFIRSTORDER, zaxisfirstorder)

inline void FORTRAN_NAME(inlineZAXISFIRSTORDER, inlineZAXISFIRSTORDER)(
      CHFp_FRA(W)
      ,CHFp_CONST_INT(lohisign)
      ,CHFp_CONST_INT(kbnd)
      ,CHFp_CONST_INT(nGS)
      ,CHFp_BOX(box) )
{
 CH_TIMELEAF("FORT_ZAXISFIRSTORDER");
 FORTRAN_NAME( ZAXISFIRSTORDER ,zaxisfirstorder )(
      CHFt_FRA(W)
      ,CHFt_CONST_INT(lohisign)
      ,CHFt_CONST_INT(kbnd)
      ,CHFt_CONST_INT(nGS)
      ,CHFt_BOX(box) );
}
#endif  // GUARDZAXISFIRSTORDER 

}

#endif
