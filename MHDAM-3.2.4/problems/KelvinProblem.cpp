#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "KelvinProblem.H"
#include "KelvinF_F.H"
#include "EqSysMHDMF.H"
#include "LoHiCenter.H"

// Null constructor
KelvinProblem::KelvinProblem()
{
    m_isFortranCommonSet = false;
}

// Input parameters
void KelvinProblem::input( ParmParse & parser, int verbosity )
{
  m_gamma      = 1.6666666667;

  parser.query( "gamma",   m_gamma   );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Kelvin-Helmholtz instability problem input:"      << endl;
    pout() << "gamma     = " << m_gamma     << endl;
  }

  setFortranCommon( m_gamma );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void KelvinProblem::setFortranCommon( const Real&     a_gamma )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETKELVIN( CHF_CONST_REAL( a_gamma     ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void KelvinProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* KelvinProblem::new_PhysProblem()
{
  KelvinProblem* retval = new KelvinProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma  = this->m_gamma;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void KelvinProblem::fluxBC(       FArrayBox&      a_F,
                                  FArrayBox&      a_Bn,
                            const FArrayBox&      a_W,
                            const FArrayBox&      a_Wextrap,
                            const int&            a_dir,
                            const Side::LoHiSide& a_side,
                            const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything
}

// Set up initial conditions
void KelvinProblem::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = m_csh->dx(0,m_level);
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
  int iCP = eqSys->correctionPotentialIndex();
  

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition in this grid
    FORT_KELVININIT(CHF_CONST_FRA(U),
                    CHF_CONST_INT(iCP),
                    CHF_CONST_REAL(dx),
                    CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void KelvinProblem::fillGhostCells(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything

}

//                            Return boundary condition flags for all boundaries
void KelvinProblem::getBCFlags( eBoundaryConditions leftBC,
                                eBoundaryConditions rightBC,
                                eBoundaryConditions bottomBC,
                                eBoundaryConditions topBC,
                                eBoundaryConditions frontBC,
                                eBoundaryConditions behindBC )
{
  leftBC   = BC_Periodic;
  rightBC  = BC_Periodic;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions KelvinProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 2 )
    return BC_Undefined;

  return BC_Periodic;
}
