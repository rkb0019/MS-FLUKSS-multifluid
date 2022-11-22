#include "PhysProblem.H"
#include "CSHandler.H"
#include "EquationSystem.H"
#include "LGintegrator.H"
#include "RiemannSolver.H"
#include "DebugF_F.H"

                                    // Indicate that define() hasn't been called
PhysProblem::PhysProblem()
{
  m_isDefined = false;

  m_pSourceC  = NULL;  
  m_RS        = NULL;  
  m_RSGD      = NULL;  
  
  m_csh       = NULL;
  m_eqSys     = NULL;
  
  m_verbosity = 0;  
  m_level     = -1;
  
  m_numGhostCells  = 0;
}

                                                            // Define the object
void PhysProblem::define( const ProblemDomain & a_domain,                          
                          const int             a_level  )
{
  m_domain = a_domain;  
  m_level  = a_level;

  m_isDefined = true;
}

// Copy method
void PhysProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  this->setSourceCalculator( a_PP->getSourceCalculator() );
  this->setRiemannSolver   ( a_PP->getRiemannSolver()  );
  this->setRiemannSolverGD ( a_PP->getRiemannSolverGD());  
  
  this->setCoordinateSystem ( a_PP->coordinateSystem() );
  this->setEquationSystem ( a_PP->equationSystem() );
  
  this->m_verbosity      = a_PP->m_verbosity;
  this->m_numGhostCells  = a_PP->m_numGhostCells;
}

// Default implementation of artificial viscosity at the boundary does nothing
void PhysProblem::artViscBC(       FArrayBox & a_F,
                             const FArrayBox & a_U,
                             const FArrayBox & a_divVel,
                             const int       & a_dir,
                             const Real      & a_time    )
{
}

                      // Default implementation of problem specific source terms
void PhysProblem::explicitSource(       FArrayBox & a_U,
                                        FArrayBox & a_S,
                                  const FArrayBox & a_W,
                                     BaseFab<int> & a_REG,
                                  const Real      & a_dt,                                  
                                  const Box       & a_box       )
{
}

                                              // Problem specific postprocessing
void PhysProblem::postprocessing(       FArrayBox & a_U,
                                  const FArrayBox & a_W,
                                  const Real      & a_dt,
                                  const Real      & a_time,                                  
                                  const Box       & a_box       )
{
}

                            // Shock boundary conditions
void PhysProblem::shockBC( const FArrayBox & a_W,
                                 FArrayBox & a_U,
                              BaseFab<int> & a_REG,
			         const int   a_level)
{
}

                             // Creates tagged cells for dynamic mesh refinement
void PhysProblem::tagCells( const FArrayBox&  a_U,
                            const Box&        a_box,
                                  IntVectSet& a_tags )
{
}

                            // Check geometrical/problem limitations for grid adaptation
void PhysProblem::lockedCellsRegrid( BaseFab<int> & a_flag,
                               const FArrayBox&  a_U,
                               const Box&     a_box)
{ 
}

void PhysProblem::defineMesh(const ProblemDomain & a_prob_domain,
                             const Vector<Real>  & a_domainBox)
{
}

/// Set up initial conditions
  /**
   */
void PhysProblem::initialize(LevelData<FArrayBox>& a_U,
                                        Interval & a_comp)
{
  MayDay::Error("This problem does not support this type of initialization [initialize(LevelData<FArrayBox>& a_U,Interval & a_comp]");
}                                        

//                                             Set the source calculator pointer
void PhysProblem::setSourceCalculator( SourceCalculator * a_pSC )
{
  if( a_pSC != NULL )
  {
    m_pSourceC = a_pSC;
  }
}
//                                             Get the source calculator pointer
SourceCalculator * PhysProblem::getSourceCalculator( void ) const
{
  return m_pSourceC;
}


void PhysProblem::setRiemannSolver( RiemannSolver * a_RS )
{   
  if( a_RS != NULL )
  {
    m_RS = a_RS;
  }    
}

//                                  Return the current MHD Riemann solver object
RiemannSolver * PhysProblem::getRiemannSolver( void ) const
{  
  return m_RS;
}

void PhysProblem::setRiemannSolverGD( RiemannSolver * a_RSGD )
{   
  if( a_RSGD != NULL )
  {
    m_RSGD = a_RSGD;
  }    
}

//                                  
RiemannSolver * PhysProblem::getRiemannSolverGD( void ) const
{  
  return m_RSGD;
}

/// Set the number of ghost cells
  /**
   */
void PhysProblem::setNumGhostCells(int a_numGhostCells)
{
  m_numGhostCells = a_numGhostCells;
}


// Number additional variables for writing to plot file
int PhysProblem::numPlotVars()
{
  return 0;
}

//               Generate default names for the primitive variables, "variable#"
Vector<std::string> PhysProblem::plotNames()
{
  Vector<std::string> retval;

  int cnum = numPlotVars();

  for (int ivar = 0; ivar < cnum; ivar++)
  {
    char varNameChar[80];
    sprintf(varNameChar,"variable%d",ivar);
    retval.push_back(std::string(varNameChar));
  }

  return retval; 
}

// Calculates variables for plotting using primitive variables  
void PhysProblem::calcPlotVars(FArrayBox&      a_vars,
                         int              a_comp,
                         const FArrayBox& a_W,
                         const Box&       a_box)
{                         
}

void PhysProblem::primForPlot(    FArrayBox& a_W,
                            const Box&       a_box)
{

}                           
                                               // Set the Coordinate System flag
void PhysProblem::setCoordinateSystem( CoordinateSystemHandler* a_csh )
{
  m_csh = a_csh;
}

                                               // Get the Coordinate System flag
CoordinateSystemHandler * PhysProblem::coordinateSystem( void ) const
{
  return m_csh;
}

                                               // Set the Equation System flag
void PhysProblem::setEquationSystem( EquationSystem* a_eqSys )
{
  m_eqSys = a_eqSys;
}

                                               // Get the Equation System flag
EquationSystem * PhysProblem::equationSystem( void ) const
{
  return m_eqSys;
}


//                            Return boundary condition flags for all boundaries
void PhysProblem::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions PhysProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  return BC_Undefined;
}

void PhysProblem::projectBbc(FArrayBox &a_state, const Box &a_valid, const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous)
{
  Box valid = a_valid;
  for(int dir=0; dir<CH_SPACEDIM; ++dir)
  {
    // don't do anything if periodic
    if (!a_domain.isPeriodic(dir)) continue;
    
    SideIterator sit;
    for (sit.begin();sit.ok();sit.next())
    {    
      int isign = sign(sit());
      Box ghostBox = adjCellBox(valid, dir, sit(), 1) & a_state.box();
      
      if(!a_domain.domainBox().contains(ghostBox))
      {
        eBoundaryConditions ebc = getBCFlags(dir, sit());
        if (ebc == BC_Fixed)
        {
          Box & toRegion = ghostBox;        

          Box   fromRegion = toRegion;
          fromRegion.shift(dir, -isign);

          a_state.copy(a_state, fromRegion, 0, toRegion, 0, a_state.nComp());          
        } else
        if (ebc == BC_Continuous)
        {          
          a_state.setVal(0.0, ghostBox, 0, a_state.nComp());          
        } else a_state.setVal(0.0, ghostBox, 0, a_state.nComp());          
      }
    }
    
  }

}


//                         Things to do before advancing one level by one time step.
void PhysProblem::preTimeStep(LevelData<FArrayBox>&       a_U, Real a_time) 
{
}

// Things to do after advancing this level by one time step.  
void PhysProblem::postTimeStep(LevelData<FArrayBox>&      a_U)
{
}

/// Check geometrical/problem limitations for grid adaptation
Real PhysProblem::computeDt( const FArrayBox& a_U,
                             const FArrayBox& a_dt,
                             const Box&     a_box,
                             IntVect&       a_minDtCell)
{
  a_minDtCell = a_dt.minIndex();
  Real dt     = a_dt.get(a_minDtCell,0);
  
  if (m_verbosity>=4)
  { 
    pout() << "local min dt =  " << dt << " at ";a_minDtCell.p();pout() << endl;
  }  
  
#ifndef NDEBUG
  if (dt <= 0.0)
    FORT_VIEWBOXDATACONST(CHF_CONST_FRA(a_U));
#endif  

  CH_assert(dt > 0.0);

  return dt;
}

Real PhysProblem::getPhysTime(Real a_time)
{
  return a_time;
}

// For a tracking surface 'a_s' returns first component of the velocity field for the level set method
int PhysProblem::lsIndexField(int a_s)
{
  return WVELX;
}
  
// Returns curvilinear coordinates for 1D probes
void PhysProblem::probeFilesHeader(int & a_numVars, std::string& a_header)
{
}

Real PhysProblem::probeFilesDimensiolessDt(Real a_dt)
{
  return a_dt;
}


void PhysProblem::conclude(const LevelData<FArrayBox>& a_U,
                        Real a_time,
                        int  a_cur_step,                        
                        const std::string & a_inputfile)
{
}


/// Creates "DATASETAUXDATA" fields for tecplot
void PhysProblem::auxDataTecplot(std::string & a_str,
                              Real          a_time,
                              int           a_type)
{
}
