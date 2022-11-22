#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "RampProblem.H"
#include "RampF_F.H"

#include "LoHiCenter.H"

// Null constructor
RampProblem::RampProblem()
{
    m_isFortranCommonSet = false;
}

// Input parameters
void RampProblem::input( ParmParse & parser, int verbosity )
{
  m_alpha      = 2.0;
  m_M          = 1.0;
  m_xcorner    = 0.115;
  m_gamma      = 1.4;

  parser.query( "gamma",      m_gamma   );
  parser.query( "alpha_deg",  m_alpha   );
  parser.query( "shock_mach", m_M       );
  parser.query( "xcorner",    m_xcorner );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Ramp problem input:"    << endl;
    pout() << "gamma    = " << m_gamma    << endl;
    pout() << "alpha    = " << m_alpha    << endl;
    pout() << "M        = " << m_M        << endl;
    pout() << "xcorner  = " << m_xcorner  << endl;
  }

  setFortranCommon( m_gamma, m_alpha, m_M, m_xcorner );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RampProblem::setFortranCommon( const Real&     a_gamma,
                                    const Real&     a_alpha,
                                    const Real&     a_M,
                                    const Real&     a_xcorner )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETRAMP( CHF_CONST_REAL( a_gamma   ),
                  CHF_CONST_REAL( a_alpha  ),
                  CHF_CONST_REAL( a_M  ),
                  CHF_CONST_REAL( a_xcorner ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RampProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* RampProblem::new_PhysProblem()
{
  RampProblem* retval = new RampProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma    = this->m_gamma;
    retval->m_alpha    = this->m_alpha;
    retval->m_M        = this->m_M;
    retval->m_xcorner  = this->m_xcorner;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void RampProblem::fluxBC(       FArrayBox&      a_F,
                                FArrayBox&      a_Bn,
                          const FArrayBox&      a_W,
                          const FArrayBox&      a_Wextrap,
                          const int&            a_dir,
                          const Side::LoHiSide& a_side,
                          const Real&           a_time )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  int sign;
  Box FBox = a_F.box();
  Box tmp = FBox;
  
  Real dx = m_csh->dx(0,m_level);

  // Determine which side and thus shifting directions
  if (a_side == Side::Lo)
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }

  tmp.shiftHalf(a_dir,sign);

  // Is there a domain boundary next to this grid
  if (!m_domain.contains(tmp))
  {
    tmp &= m_domain;

    Box boundaryBox;

    // Find the strip of cells next to the domain boundary
    if (a_side == Side::Lo)
    {
      boundaryBox = bdryLo(tmp,a_dir);
    }
    else
    {
      boundaryBox = bdryHi(tmp,a_dir);
    }

    // Shift things to all line up correctly
    boundaryBox.shiftHalf(a_dir,-sign);
    a_F.shiftHalf(a_dir,-sign);

    // Set the solid wall boundary fluxes
    FORT_RAMPSOLIDBC(CHF_FRA(a_F),
                     CHF_CONST_FRA(a_Wextrap),
                     CHF_CONST_INT(sign),
                     CHF_CONST_REAL(dx),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(boundaryBox));

    // Set the time varying boundary fluxes
    FORT_RAMPBC(CHF_FRA(a_F),
                CHF_CONST_FRA(a_W),
                CHF_CONST_REAL(a_time),
                CHF_CONST_INT(sign),
                CHF_CONST_REAL(dx),
                CHF_CONST_INT(a_dir),
                CHF_BOX(boundaryBox));

    // Shift returned fluxes to be face centered
    a_F.shiftHalf(a_dir,sign);
  }
}

// Set up initial conditions
void RampProblem::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = m_csh->dx(0,m_level);

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition in this grid
    FORT_RAMPINIT( CHF_CONST_FRA(U),
                   CHF_CONST_REAL(dx),
                   CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void RampProblem::fillGhostCells(       FArrayBox&      a_W,
                                  const FArrayBox&      a_U,
                                  const int&            a_dir,
                                  const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();

  Box tmp  = WBox;
  tmp     &= m_domain;
  
  Real dx = m_csh->dx(0,m_level);

                         // See if this chops off the high side of the input box
  int indW = WBox.bigEnd( a_dir );
  int indD =  tmp.bigEnd( a_dir );

  if( indW != indD )
  {
    int sign         = 1;
    Box boundaryBox  = WBox;
    boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                         // Fill the ghost cells
    FORT_RAMPGS( CHF_FRA(a_W),
                 CHF_CONST_REAL(a_time),
                 CHF_CONST_REAL(dx),
                 CHF_CONST_INT(sign),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(boundaryBox) );
  }

  indW     = WBox.smallEnd( a_dir );
  indD     =  tmp.smallEnd( a_dir );

  if( indW != indD )
  {
    int sign         =-1;
    Box boundaryBox  = WBox;
    boundaryBox.setRange( a_dir, indW, indD - indW );

                                                         // Fill the ghost cells
    FORT_RAMPGS( CHF_FRA(a_W),
                 CHF_CONST_REAL(a_time),
                 CHF_CONST_REAL(dx),
                 CHF_CONST_INT(sign),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(boundaryBox) );
  }
}

//                            Return boundary condition flags for all boundaries
void RampProblem::getBCFlags( eBoundaryConditions leftBC,
                              eBoundaryConditions rightBC,
                              eBoundaryConditions bottomBC,
                              eBoundaryConditions topBC,
                              eBoundaryConditions frontBC,
                              eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Continuous;
  bottomBC = BC_Axis;
  topBC    = BC_Mixed;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions RampProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  case 1  : return (a_sd == Side::Lo) ? BC_Axis  : BC_Mixed;
  default : return BC_Undefined;
  }
}
