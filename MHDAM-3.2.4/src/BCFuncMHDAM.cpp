#include "BCFuncMHDAM.H"
#include "PhysProblem.H"


void BCFuncMHDAM::operator()(FArrayBox&        a_state,
                          const Box&           a_valid,
                          const ProblemDomain& a_domain,
                          Real                 a_dx,
                          bool                 a_homogeneous)
{
  m_physProblem->projectBbc(a_state, a_valid, a_domain, a_dx, a_homogeneous);
}

