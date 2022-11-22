#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "RotorProblem.H"
#include "rotorF_F.H"
#include "EqSysMHDMF.H"
#include "LoHiCenter.H"

// Null constructor
RotorProblem::RotorProblem()
{
    m_isFortranCommonSet = false;
}

// Input parameters
void RotorProblem::input( ParmParse & parser, int verbosity )
{
  m_rotorU     = 2.0;
  m_rotorP     = 1.0;
  m_rotorR     = 10.0;
  m_rotorB     = 5.0;
  m_rotorXC    = 0.5;
  m_rotorYC    = 0.5;
  m_rotorR0    = 0.1;
  m_rotorR1    = 0.115;
  m_gamma      = 1.6666666667;

  parser.query( "gamma",   m_gamma   );
  parser.query( "rotorU",  m_rotorU  );
  parser.query( "rotorP",  m_rotorP  );
  parser.query( "rotorR",  m_rotorR  );
  parser.query( "rotorB",  m_rotorB  );
  parser.query( "rotorXC", m_rotorXC );
  parser.query( "rotorYC", m_rotorYC );
  parser.query( "rotorR0", m_rotorR0 );
  parser.query( "rotorR1", m_rotorR1 );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Rotor problem input:"    << endl;
    pout() << "gamma    = " << m_gamma    << endl;
    pout() << "rotorU   = " << m_rotorU   << endl;
    pout() << "rotorP   = " << m_rotorP   << endl;
    pout() << "rotorR   = " << m_rotorR   << endl;
    pout() << "rotorB   = " << m_rotorB   << endl;
    pout() << "rotorXC  = " << m_rotorXC  << endl;
    pout() << "rotorYC  = " << m_rotorYC  << endl;
    pout() << "rotorR0  = " << m_rotorR0  << endl;
    pout() << "rotorR1  = " << m_rotorR1  << endl;
  }

  setFortranCommon( m_gamma,
                    m_rotorU,  m_rotorP,  m_rotorR,  m_rotorB,
                    m_rotorXC, m_rotorYC, m_rotorR0, m_rotorR1 );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RotorProblem::setFortranCommon( const Real&     a_gamma,
                                     const Real&     a_rotorU,
                                     const Real&     a_rotorP,
                                     const Real&     a_rotorR,
                                     const Real&     a_rotorB,
                                     const Real&     a_rotorXC,
                                     const Real&     a_rotorYC,
                                     const Real&     a_rotorR0,
                                     const Real&     a_rotorR1 )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETROTOR( CHF_CONST_REAL( a_gamma   ),
                   CHF_CONST_REAL( a_rotorU  ),
                   CHF_CONST_REAL( a_rotorP  ),
                   CHF_CONST_REAL( a_rotorR  ),
                   CHF_CONST_REAL( a_rotorB  ),
                   CHF_CONST_REAL( a_rotorXC ),
                   CHF_CONST_REAL( a_rotorYC ),
                   CHF_CONST_REAL( a_rotorR0 ),
                   CHF_CONST_REAL( a_rotorR1 ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RotorProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* RotorProblem::new_PhysProblem()
{
  RotorProblem* retval = new RotorProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma    = this->m_gamma;
    retval->m_rotorU   = this->m_rotorU;
    retval->m_rotorP   = this->m_rotorP;
    retval->m_rotorR   = this->m_rotorR;
    retval->m_rotorB   = this->m_rotorB;
    retval->m_rotorXC  = this->m_rotorXC;
    retval->m_rotorYC  = this->m_rotorYC;
    retval->m_rotorR0  = this->m_rotorR0;
    retval->m_rotorR1  = this->m_rotorR1;

    retval->setFortranCommonSet();
  }
  
 
  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void RotorProblem::fluxBC(       FArrayBox&      a_F,
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
void RotorProblem::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_ROTORINIT(CHF_CONST_FRA(U),
                    CHF_CONST_INT(iCP),
                    CHF_CONST_REAL(dx),
                    CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void RotorProblem::fillGhostCells(       FArrayBox&      a_W,
                                   const FArrayBox&      a_U,
                                   const int&            a_dir,
                                   const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything

}

//                            Return boundary condition flags for all boundaries
void RotorProblem::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions RotorProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 2 )
    return BC_Undefined;

  return BC_Periodic;
}
