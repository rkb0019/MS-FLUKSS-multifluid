#if CH_LANG_CC
// Constants used between C++ and Fortran to index the conserved and
// source terms and to know the number of each.  
#endif

#ifndef CH_KS2DINTEGRATOR_H
#define CH_KS2DINTEGRATOR_H

#if CH_LANG_CC
//  Variables needed by source calculator
// 'K' means kinetic plasma data
#endif

#if CH_SPACEDIM == 2
  #define KCOMP  5

  #define KRHO   0 
  #define KVELX  1 
  #define KVELY  2 
  #define KVELZ  KVELY
  #define KTEMP  3
  #define KREG   4

  #if CH_LANG_CC
  // Variables that are stored in AMRLevelIdealMHD::m_SCData
  //  NRHO  -  neutral density
  //  NVEL  -  neutral velocity
  //  NTEMP -  neutral temperature
  #endif

  #define SCOMP    4
  
  #define SRHO     0
  #define SMOMX    1
  #define SMOMY    2
  #define SMOMZ    SMOMY
  #define SENG     3
  
  #define NCOMP    17
  
  #define NCHEX    4
  #define NRHO     5 
  #define NVELX    6 
  #define NVELY    7 
  #define NVELZ    NVELY
  #define NTEMP    8 
  #define NRHO1    9 
  #define NVELX1   10 
  #define NVELY1   11 
  #define NVELZ1   NVELY1   
  #define NTEMP1   12 
  #define NRHO2    13 
  #define NVELX2   14 
  #define NVELY2   15 
  #define NVELZ2   NVELY2   
  #define NTEMP2   16 
  #define NRHO3    17 
  #define NVELX3   18 
  #define NVELY3   19 
  #define NVELZ3   NVELY3   
  #define NTEMP3   20 
#endif

#if CH_SPACEDIM == 3

  #define KCOMP  6
  
  #define KRHO   0 
  #define KVELX  1 
  #define KVELY  2 
  #define KVELZ  3 
  #define KTEMP  4
  #define KREG   5

  #if CH_LANG_CC
  // Variables that are stored in AMRLevelIdealMHD::m_SCData
  //  NRHO  -  neutral density
  //  NVEL  -  neutral velocity
  //  NTEMP -  neutral temperature
  #endif

  #define SCOMP    5
  
  #define SRHO     0
  #define SMOMX    1
  #define SMOMY    2
  #define SMOMZ    3
  #define SENG     4
  
  #define NCOMP    6
  
  #define NCHEX    5
  #define NRHO     6 
  #define NVELX    7 
  #define NVELY    8 
  #define NVELZ    9 
  #define NTEMP    10 
  
  #if CH_LANG_CC
  // Currently not used. Declared for compilation of 2d code  
  #endif

  #define NRHO1    11 
  #define NVELX1   12 
  #define NVELY1   13 
  #define NVELZ1   14
  #define NTEMP1   15 
  #define NRHO2    16 
  #define NVELX2   17 
  #define NVELY2   18 
  #define NVELZ2   19   
  #define NTEMP2   20 
  #define NRHO3    21 
  #define NVELX3   22 
  #define NVELY3   23 
  #define NVELZ3   24   
  #define NTEMP3   25  
#endif




#endif
