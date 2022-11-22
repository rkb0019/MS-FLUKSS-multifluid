#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "RiemannProblem.H"
#include "RiemannF_F.H"

#include "LoHiCenter.H"

// Null constructor
RiemannProblem::RiemannProblem()
{
    m_isFortranCommonSet = false;
}

// Input parameters
void RiemannProblem::input( ParmParse & parser, int verbosity )
{
  m_densityL   = 1.0;
  m_densityR   = 0.125;
  m_pressureL  = 1.0;
  m_pressureR  = 0.1;
  m_velxL      = 0.0;
  m_velxR      = 0.0;
  m_velyL      = 0.0;
  m_velyR      = 0.0;
  m_velzL      = 0.0;
  m_velzR      = 0.0;
  m_BxL        = 0.0;
  m_BxR        = 0.0;
  m_ByL        = 0.0;
  m_ByR        = 0.0;
  m_BzL        = 0.0;
  m_BzR        = 0.0;
  m_startX     = 0.5;
  m_gamma      = 1.6666666667;

  parser.query( "gamma",     m_gamma     );
  parser.query( "densityL",  m_densityL  );
  parser.query( "densityR",  m_densityR  );
  parser.query( "pressureL", m_pressureL );
  parser.query( "pressureR", m_pressureR );
  parser.query( "velxL",     m_velxL     );
  parser.query( "velxR",     m_velxR     );
  parser.query( "velyL",     m_velyL     );
  parser.query( "velyR",     m_velyR     );
  parser.query( "velzL",     m_velzL     );
  parser.query( "velzR",     m_velzR     );
  parser.query( "BxL",       m_BxL       );
  parser.query( "BxR",       m_BxR       );
  parser.query( "ByL",       m_ByL       );
  parser.query( "ByR",       m_ByR       );
  parser.query( "BzL",       m_BzL       );
  parser.query( "BzR",       m_BzR       );
  parser.query( "startX",    m_startX    );

                                                         // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "1D Riemann problem input:"   << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "densityL  = " << m_densityL  << endl;
    pout() << "densityR  = " << m_densityR  << endl;
    pout() << "pressureL = " << m_pressureL << endl;
    pout() << "pressureR = " << m_pressureR << endl;
    pout() << "velxL     = " << m_velxL     << endl;
    pout() << "velxR     = " << m_velxR     << endl;
    pout() << "velyL     = " << m_velyL     << endl;
    pout() << "velyR     = " << m_velyR     << endl;
    pout() << "velzL     = " << m_velzL     << endl;
    pout() << "velzR     = " << m_velzR     << endl;
    pout() << "BxL       = " << m_BxL       << endl;
    pout() << "BxR       = " << m_BxR       << endl;
    pout() << "ByL       = " << m_ByL       << endl;
    pout() << "ByR       = " << m_ByR       << endl;
    pout() << "BzL       = " << m_BzL       << endl;
    pout() << "BzR       = " << m_BzR       << endl;
    pout() << "startX    = " << m_startX    << endl;
  }

  setFortranCommon( m_gamma,
                    m_densityL,   m_densityR,
                    m_pressureL,  m_pressureR,
                    m_velxL,      m_velxR,
                    m_velyL,      m_velyR,
                    m_velzL,      m_velzR,
                    m_BxL,        m_BxR,
                    m_ByL,        m_ByR,
                    m_BzL,        m_BzR,       m_startX    );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RiemannProblem::setFortranCommon( const Real&     a_gamma,
                                       const Real&     a_densityL,   const Real&     a_densityR,
                                       const Real&     a_pressureL,  const Real&     a_pressureR,
                                       const Real&     a_velxL,      const Real&     a_velxR,
                                       const Real&     a_velyL,      const Real&     a_velyR,
                                       const Real&     a_velzL,      const Real&     a_velzR,
                                       const Real&     a_BxL,        const Real&     a_BxR,
                                       const Real&     a_ByL,        const Real&     a_ByR,
                                       const Real&     a_BzL,        const Real&     a_BzR,
                                       const Real&     a_startX                                  )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETRIEMANN( CHF_CONST_REAL( a_gamma     ),
                     CHF_CONST_REAL( a_densityL  ),
                     CHF_CONST_REAL( a_densityR  ),
                     CHF_CONST_REAL( a_pressureL ),
                     CHF_CONST_REAL( a_pressureR ),
                     CHF_CONST_REAL( a_velxL     ),
                     CHF_CONST_REAL( a_velxR     ),
                     CHF_CONST_REAL( a_velyL     ),
                     CHF_CONST_REAL( a_velyR     ),
                     CHF_CONST_REAL( a_velzL     ),
                     CHF_CONST_REAL( a_velzR     ),
                     CHF_CONST_REAL( a_BxL       ),
                     CHF_CONST_REAL( a_BxR       ),
                     CHF_CONST_REAL( a_ByL       ),
                     CHF_CONST_REAL( a_ByR       ),
                     CHF_CONST_REAL( a_BzL       ),
                     CHF_CONST_REAL( a_BzR       ),
                     CHF_CONST_REAL( a_startX    ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RiemannProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* RiemannProblem::new_PhysProblem()
{
  RiemannProblem* retval = new RiemannProblem();

  if (m_isFortranCommonSet == true)
  {
    retval->m_gamma      = this->m_gamma;
    retval->m_densityL   = this->m_densityL;
    retval->m_densityR   = this->m_densityR;
    retval->m_pressureL  = this->m_pressureL;
    retval->m_pressureR  = this->m_pressureR;
    retval->m_velxL      = this->m_velxL;
    retval->m_velxR      = this->m_velxR;
    retval->m_velyL      = this->m_velyL;
    retval->m_velyR      = this->m_velyR;
    retval->m_velzL      = this->m_velzL;
    retval->m_velzR      = this->m_velzR;
    retval->m_BxL        = this->m_BxL;
    retval->m_BxR        = this->m_BxR;
    retval->m_ByL        = this->m_ByL;
    retval->m_ByR        = this->m_ByR;
    retval->m_BzL        = this->m_BzL;
    retval->m_BzR        = this->m_BzR;
    retval->m_startX     = this->m_startX;

    retval->setFortranCommonSet();
  }
  
  retval->copy_PhysProblem(this);

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void RiemannProblem::fluxBC(       FArrayBox&      a_F,
                                   FArrayBox&      a_Bn,
                             const FArrayBox&      a_W,
                             const FArrayBox&      a_Wextrap,
                             const int&            a_dir,
                             const Side::LoHiSide& a_side,
                             const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);
    
  int iRhoN  = 0;
  int fluids = 1;

  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    int sign;
    Box FBox = a_F.box();
    Box tmp = FBox;

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
      boundaryBox.shiftHalf( a_dir,-sign );
      a_F.shiftHalf( a_dir,-sign );

      
                                                      // Set the boundary fluxes
      FORT_FLUXBC( CHF_FRA(a_F),
                   CHF_FRA1(a_Bn,0),
                   CHF_CONST_FRA(a_W),
                   CHF_CONST_INT(sign),
                   CHF_CONST_REAL(dx),
                   CHF_CONST_INT(a_dir),
                   CHF_CONST_INT(iRhoN),
                   CHF_CONST_INT(fluids),
                   CHF_BOX(boundaryBox) );

                                    // Shift returned fluxes to be face centered
      a_F.shiftHalf( a_dir, sign );
    }
  }
}

// Set up initial conditions
void RiemannProblem::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_RIEMANNINIT(CHF_CONST_FRA(U),
                     CHF_CONST_REAL(dx),
                     CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void RiemannProblem::fillGhostCells(       FArrayBox&      a_W,
                                     const FArrayBox&      a_U,
                                     const int&            a_dir,
                                     const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int iRhoN  = 0;
  int fluids = 1;

                                   // In periodic case, this doesn't do anything
  if( !m_domain.isPeriodic(a_dir) )
  {
    Box WBox = a_W.box();

                         // See if this chops off the high side of the input box
    Box tmp  = WBox;
    tmp     &= m_domain;

    int fluids = 1;

    int indW = WBox.bigEnd( a_dir );
    int indD =  tmp.bigEnd( a_dir );

    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                         // Fill the ghost cells
      FORT_RIEMANNGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
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
      FORT_RIEMANNGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void RiemannProblem::getBCFlags( eBoundaryConditions leftBC,
                                 eBoundaryConditions rightBC,
                                 eBoundaryConditions bottomBC,
                                 eBoundaryConditions topBC,
                                 eBoundaryConditions frontBC,
                                 eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Fixed;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions RiemannProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return BC_Fixed;
  case 1  : return BC_Periodic;
  default : return BC_Undefined;
  }
}
