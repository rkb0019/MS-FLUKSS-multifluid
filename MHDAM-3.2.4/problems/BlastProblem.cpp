#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "BlastProblem.H"
#include "BlastF_F.H"
#include "EqSysMHDMF.H"
#include "LoHiCenter.H"

// Null constructor
BlastProblem::BlastProblem()
{
  m_isFortranCommonSet = false;
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void BlastProblem::setFortranCommon( const Real&     a_gamma,
                                     const Real&     a_Blastp0,
                                     const Real&     a_Blastp1,
                                     const Real&     a_BlastXc,
                                     const Real&     a_BlastYc,  
                                     const Real&     a_BlastZc,             
                                     const Real&     a_BlastBx,
                                     const Real&     a_BlastBy,
                                     const Real&     a_Blastr0 )
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_SETBLAST( CHF_CONST_REAL( a_gamma    ),
                 CHF_CONST_REAL( a_Blastp0  ),
                 CHF_CONST_REAL( a_Blastp1  ),
                 CHF_CONST_REAL( a_BlastXc  ),
                 CHF_CONST_REAL( a_BlastYc  ),
                 CHF_CONST_REAL( a_BlastZc  ),
                 CHF_CONST_REAL( a_BlastBx  ),
                 CHF_CONST_REAL( a_BlastBy  ),
                 CHF_CONST_REAL( a_Blastr0  ) );

  m_isFortranCommonSet = true;
}

// Input parameters
void BlastProblem::input( ParmParse & parser, int verbosity )
{
  m_Blastp0  = 1.0;
  m_Blastp1  = 1000.0;
  m_BlastXc  = 0.5;
  m_BlastYc  = 0.5;
  m_BlastZc  = 0.0;      
  m_BlastBx  = 30.0;
  m_BlastBy  = 0.0;
  m_Blastr0  = 0.4;

  m_gamma    = 1.6666666667;

  parser.query( "gamma",   m_gamma   );
  parser.query( "Blastp0", m_Blastp0 );
  parser.query( "Blastp1", m_Blastp1 );
  parser.query( "BlastXc", m_BlastXc );
  parser.query( "BlastYc", m_BlastYc );
  parser.query( "BlastZc", m_BlastZc );
  parser.query( "BlastBx", m_BlastBx );
  parser.query( "BlastBy", m_BlastBy );
  parser.query( "Blastr0", m_Blastr0 );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Blast problem input:"      << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "Blastp0   = " << m_Blastp0   << endl;
    pout() << "Blastp1   = " << m_Blastp1   << endl;
    pout() << "BlastXc   = " << m_BlastXc   << endl;
    pout() << "BlastYc   = " << m_BlastYc   << endl;
    pout() << "BlastZc   = " << m_BlastZc   << endl;
    pout() << "BlastBx   = " << m_BlastBx   << endl;
    pout() << "BlastBy   = " << m_BlastBy   << endl;
    pout() << "Blastr0   = " << m_Blastr0   << endl;
  }

  setFortranCommon( m_gamma,
                    m_Blastp0, m_Blastp1,
                    m_BlastXc, m_BlastYc, m_BlastZc,             
                    m_BlastBx, m_BlastBy,
                    m_Blastr0 );
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void BlastProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* BlastProblem::new_PhysProblem()
{
  BlastProblem* retval = new BlastProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma    = this->m_gamma;
    retval->m_Blastp0  = this->m_Blastp0;
    retval->m_Blastp1  = this->m_Blastp1;
    retval->m_BlastXc  = this->m_BlastXc;
    retval->m_BlastYc  = this->m_BlastYc;
    retval->m_BlastZc  = this->m_BlastZc;      
    retval->m_BlastBx  = this->m_BlastBx;
    retval->m_BlastBy  = this->m_BlastBy;
    retval->m_Blastr0  = this->m_Blastr0;

    retval->setFortranCommonSet();
  }
  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void BlastProblem::fluxBC( FArrayBox&            a_F,
                           FArrayBox&            a_Bn,
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
void BlastProblem::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
  int iCP = eqSys->correctionPotentialIndex();

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    
    // Set up initial condition in this grid
    FORT_BLASTINIT( CHF_CONST_FRA(U),
                    CHF_CONST_INT(iCP),
                    CHF_CONST_REAL(dx),
                    CHF_BOX(uBox));
  }
}

                                                             // Fill ghost cells
void BlastProblem::fillGhostCells(       FArrayBox&      a_W,
                                   const FArrayBox&      a_U,
                                   const int&            a_dir,
                                   const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything

}

//                            Return boundary condition flags for all boundaries
void BlastProblem::getBCFlags( eBoundaryConditions leftBC,
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
  frontBC  = BC_Periodic;
  behindBC = BC_Periodic;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions BlastProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  return BC_Periodic;
}
