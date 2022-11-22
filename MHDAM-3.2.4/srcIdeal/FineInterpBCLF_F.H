#ifndef _FINEINTERPBCLF_F_H_
#define _FINEINTERPBCLF_F_H_

#include "FORT_PROTO.H"
#include "CH_Timer.H"
#include "REAL.H"

extern "C"
{

#ifndef GUARDINTERPLIMIT_BCL 
#define GUARDINTERPLIMIT_BCL 
// Prototype for Fortran procedure interplimit_BCL ...
//
void FORTRAN_NAME( INTERPLIMIT_BCL ,interplimit_bcl )(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) );

#define FORT_INTERPLIMIT_BCL FORTRAN_NAME( inlineINTERPLIMIT_BCL, inlineINTERPLIMIT_BCL)
#define FORTNT_INTERPLIMIT_BCL FORTRAN_NAME( INTERPLIMIT_BCL, interplimit_bcl)

inline void FORTRAN_NAME(inlineINTERPLIMIT_BCL, inlineINTERPLIMIT_BCL)(
      CHFp_FRA(islope)
      ,CHFp_FRA(jslope)
      ,CHFp_FRA(kslope)
      ,CHFp_CONST_FRA(state)
      ,CHFp_BOX(bcoarse)
      ,CHFp_BOX(bn)
      ,CHFp_BOX(physdomain) )
{
 CH_TIMELEAF("FORT_INTERPLIMIT_BCL");
 FORTRAN_NAME( INTERPLIMIT_BCL ,interplimit_bcl )(
      CHFt_FRA(islope)
      ,CHFt_FRA(jslope)
      ,CHFt_FRA(kslope)
      ,CHFt_CONST_FRA(state)
      ,CHFt_BOX(bcoarse)
      ,CHFt_BOX(bn)
      ,CHFt_BOX(physdomain) );
}
#endif  // GUARDINTERPLIMIT_BCL 

}

#endif
