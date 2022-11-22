#if CH_LANG_CC
// Constants used between C++ and Fortran to index the conserved and
// source terms and to know the number of each.  
#endif

#ifndef CH_KS1DINTEGRATOR_H
#define CH_KS1DINTEGRATOR_H

#if CH_LANG_CC
// Variables needed by source calculator
#endif

#define SWRHO   0 
#define SWVELX  1 
#define SWVELY  2 
#define SWTEMP  3 

#if CH_LANG_CC
// Variables that are stored in AMRLevelIdealMHD::m_KineticData
#endif

#define SMOMX  0
#define SMOMY  1
#define SENG   2
#define NRHO   3 
#define NVELX  4 
#define NVELY  5 
#define NTEMP  6 

#if CH_LANG_CC
//  NRHO  -  neutral density
//  NVEL  -  neutral velocity
//  NTEMP -  neutral temperature
#endif



#endif
