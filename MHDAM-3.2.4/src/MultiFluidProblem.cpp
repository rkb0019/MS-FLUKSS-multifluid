#include "MultiFluidProblem.H"

                                    // Indicate that define() hasn't been called
MultiFluidProblem::MultiFluidProblem()
  : PhysProblem()
{
  m_isDefined = false;
  m_physModel = PP_Undefined;
}

                                                            // Define the object
void MultiFluidProblem::define( const ProblemDomain & a_domain,
                                const int             a_level  )
{
  m_domain = a_domain;  
  m_level  = a_level;

  m_isDefined = true;
}

void MultiFluidProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  this->PhysProblem::copy_PhysProblem(a_PP);
  
  const MultiFluidProblem* PP = dynamic_cast<const MultiFluidProblem*>(a_PP);
  
  this->m_physModel = PP->m_physModel;
  
}

  // Default implementation of artificial viscosity at the boundary does nothing
void MultiFluidProblem::artViscBC(       FArrayBox & a_F,
                                   const FArrayBox & a_U,
                                   const FArrayBox & a_divVel,
                                   const int       & a_dir,
                                   const Real      & a_time    )
{
}

                      // Default implementation of problem specific source terms
void MultiFluidProblem::explicitSource(       FArrayBox & a_U,
                                              FArrayBox & a_S,
                                        const FArrayBox & a_W,
                                           BaseFab<int> & a_REG,
                                        const Real      & a_dt,
                                        const Box       & a_box       )
{
}

                                              // Problem specific postprocessing
void MultiFluidProblem::postprocessing(       FArrayBox & a_U,
                                        const FArrayBox & a_W,
                                        const Real      & a_dt,
                                        const Real      & a_time,                                        
                                        const Box       & a_box       )
{
}

                             // Creates tagged cells for dynamic mesh refinement
void MultiFluidProblem::tagCells( const FArrayBox&  a_U,
                                  const Box&        a_box,
                                        IntVectSet& a_tags )
{
}
                                                      // Return number of fluids
int MultiFluidProblem::nFluids()
{
  switch (m_physModel) 
  {
    case PP_MHDPM:    return 1;break;
    case PP_2FluidPM: return 2;break; 
    case PP_3FluidPM: return 3;break;
    case PP_4FluidPM: return 4;break; 
    case PP_5FluidPM: return 5;break; 
  }
  
  MayDay::Error("MultiFluidProblem::nFluids: undefined number of fluids");
  
  return -1;
}

//                            Return boundary condition flags for all boundaries
void MultiFluidProblem::getBCFlags( eBoundaryConditions leftBC,
                                    eBoundaryConditions rightBC,
                                    eBoundaryConditions bottomBC,
                                    eBoundaryConditions topBC,
                                    eBoundaryConditions frontBC,
                                    eBoundaryConditions behindBC )
{
  leftBC   = BC_Undefined;
  rightBC  = BC_Undefined;
  bottomBC = BC_Undefined;
  topBC    = BC_Undefined;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions MultiFluidProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  return BC_Undefined;
}
