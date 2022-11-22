#if CH_LANG_CC
// Constants used between C++ and Fortran to index the conserved and
// primitive variables and to know the number of each.  Also, the number
// of slopes and fluxes are defined here.
#endif

#ifndef CH_LGINTEGRATOR_H
#define CH_LGINTEGRATOR_H

#define WRHO   0
#define WVELX  1
#define WVELY  2
#define WVELZ  3
#define WPRES  4
#define WBX    5
#define WBY    6
#define WBZ    7

#define WVELR  1
#define WVELP  2
#define WVELT  3
#define WBR    5
#define WBP    6
#define WBT    7


#define URHO   0
#define UMOMX  1
#define UMOMY  2
#define UMOMZ  3
#define UENG   4
#define UBX    5
#define UBY    6
#define UBZ    7


#define UNUM   8
#define WNUM   8
#define WNUM_E 5
#define UNUM_E 5

#define WSLOPE 8
#define WSLOPE 8
#define WSLO_E 5
#define UFLU_E 5
#define UFLUX  8

#endif
