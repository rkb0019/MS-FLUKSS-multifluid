#ifndef _REF_CRITERIA_H_
#define _REF_CRITERIA_H_

#include "REAL.H"
#include <limits.h>

#define REF_RHO 0
#define REF_B   1

#define REF_NUM 2


struct RefineCriterion
{
  bool m_PerformRef;
  Real m_Threshold;
  int  m_maxlevel; // The finest level of the refinement where criterion is apllied to 
                   // Thus, m_maxlevel+1 is the finest level which can be created using this criteria
  
  RefineCriterion()
  {
    m_PerformRef = false;
    m_Threshold  = 0.0;    
    m_maxlevel   = INT_MAX;
  }
};


#endif
