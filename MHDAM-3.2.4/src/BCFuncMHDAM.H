#ifndef _BCFUNCMHDAM_H_
#define _BCFUNCMHDAM_H_

#include <BCFunc.H>

class PhysProblem;

class BCFuncMHDAM : public BCFunction
{
public:
  BCFuncMHDAM(PhysProblem * a_physProblem){m_physProblem = a_physProblem;}

  ~BCFuncMHDAM(){}

  virtual void operator()(FArrayBox&           a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous);

protected:
  PhysProblem * m_physProblem;
};

#endif