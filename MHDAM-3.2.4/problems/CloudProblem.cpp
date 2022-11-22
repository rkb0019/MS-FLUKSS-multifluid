#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "CloudProblem.H"
#include "cloudF_F.H"
#include "EqSysMHDMF.H"
#include "LoHiCenter.H"
#include "RiemannSolver.H"
#include "LGintegrator.H"

// Null constructor
CloudProblem::CloudProblem()
{
    m_isFortranCommonSet = false;
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void CloudProblem::setFortranCommon( const Real&     a_gamma,
                                     const Real&     a_cloudRL,
                                     const Real&     a_cloudUL,
                                     const Real&     a_cloudVL,
                                     const Real&     a_cloudWL,
                                     const Real&     a_cloudPL,
                                     const Real&     a_cloudBXL,
                                     const Real&     a_cloudBYL,
                                     const Real&     a_cloudBZL,
                                     const Real&     a_cloudRR,
                                     const Real&     a_cloudUR,
                                     const Real&     a_cloudPR,
                                     const Real&     a_cloudBYR,
                                     const Real&     a_cloudBZR,
                                     const Real&     a_cloudXS,
                                     const Real&     a_cloudXC,
                                     const Real&     a_cloudYC,
                                     const Real&     a_cloudZC,
                                     const Real&     a_cloudR0,
                                     const Real&     a_cloudRho )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETCLOUD( CHF_CONST_REAL( a_gamma    ),
                   CHF_CONST_REAL( a_cloudRL  ),
                   CHF_CONST_REAL( a_cloudUL  ),
                   CHF_CONST_REAL( a_cloudVL  ),
                   CHF_CONST_REAL( a_cloudWL  ),
                   CHF_CONST_REAL( a_cloudPL  ),
                   CHF_CONST_REAL( a_cloudBXL ),
                   CHF_CONST_REAL( a_cloudBYL ),
                   CHF_CONST_REAL( a_cloudBZL ),
                   CHF_CONST_REAL( a_cloudRR  ),
                   CHF_CONST_REAL( a_cloudUR  ),
                   CHF_CONST_REAL( a_cloudPR  ),
                   CHF_CONST_REAL( a_cloudBYR ),
                   CHF_CONST_REAL( a_cloudBZR ),
                   CHF_CONST_REAL( a_cloudXS  ),
                   CHF_CONST_REAL( a_cloudXC  ),
                   CHF_CONST_REAL( a_cloudYC  ),
                   CHF_CONST_REAL( a_cloudZC  ),
                   CHF_CONST_REAL( a_cloudR0  ),
                   CHF_CONST_REAL( a_cloudRho ) );

    m_isFortranCommonSet = true;
}

// Input parameters
void CloudProblem::input( ParmParse & parser, int verbosity )
{
  m_densityL   = 1.0;
  m_densityR   = 0.125;
  m_pressureL  = 1.0;
  m_pressureR  = 0.1;
  m_velxL      = 0.0;
  m_velxR      = 0.0;
  m_velyL      = 0.0;
  m_velzL      = 0.0;
  m_BxL        = 0.0;
  m_ByL        = 0.0;
  m_ByR        = 0.0;
  m_BzL        = 0.0;
  m_BzR        = 0.0;
  m_startX     = 0.5;

  m_cloudRho   = 10.0;
  m_XC         = 0.5;
  m_YC         = 0.5;
  m_ZC         = 0.0;
  m_R0         = 0.1;

  m_gamma      = 1.6666666667;

  parser.query( "gamma",     m_gamma     );

  parser.query( "densityL",  m_densityL  );
  parser.query( "velxL",     m_velxL     );
  parser.query( "velyL",     m_velyL     );
  parser.query( "velzL",     m_velzL     );
  parser.query( "pressureL", m_pressureL );
  parser.query( "BxL",       m_BxL       );
  parser.query( "ByL",       m_ByL       );
  parser.query( "BzL",       m_BzL       );

  parser.query( "densityR",  m_densityR  );
  parser.query( "pressureR", m_pressureR );
  parser.query( "velxR",     m_velxR     );
  parser.query( "ByR",       m_ByR       );
  parser.query( "BzR",       m_BzR       );

  parser.query( "startX",    m_startX    );
  parser.query( "XC",        m_XC        );
  parser.query( "YC",        m_YC        );
  parser.query( "ZC",        m_ZC        );
  parser.query( "R0",        m_R0        );
  parser.query( "densityC",  m_cloudRho  );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Cloud problem input:"      << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "densityL  = " << m_densityL  << endl;
    pout() << "densityR  = " << m_densityR  << endl;
    pout() << "pressureL = " << m_pressureL << endl;
    pout() << "pressureR = " << m_pressureR << endl;
    pout() << "velxL     = " << m_velxL     << endl;
    pout() << "velyL     = " << m_velyL     << endl;
    pout() << "velzL     = " << m_velzL     << endl;
    pout() << "velxR     = " << m_velxR     << endl;
    pout() << "BxL       = " << m_BxL       << endl;
    pout() << "ByL       = " << m_ByL       << endl;
    pout() << "ByR       = " << m_ByR       << endl;
    pout() << "BzL       = " << m_BzL       << endl;
    pout() << "BzR       = " << m_BzR       << endl;
    pout() << "startX    = " << m_startX    << endl;
    pout() << "XC        = " << m_XC        << endl;
    pout() << "YC        = " << m_YC        << endl;
    pout() << "ZC        = " << m_ZC        << endl;
    pout() << "R0        = " << m_R0        << endl;
    pout() << "densityC  = " << m_cloudRho  << endl;
  }

  setFortranCommon( m_gamma,
                    m_densityL, m_velxL, m_velyL, m_velzL, m_pressureL, m_BxL, m_ByL, m_BzL,
                    m_densityR, m_velxR, m_pressureR, m_ByR, m_BzR,
                    m_startX, m_XC, m_YC, m_ZC, m_R0, m_cloudRho );
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void CloudProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* CloudProblem::new_PhysProblem()
{
  CloudProblem* retval = new CloudProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma      = this->m_gamma;
    retval->m_densityL   = this->m_densityL;
    retval->m_densityR   = this->m_densityR;
    retval->m_pressureL  = this->m_pressureL;
    retval->m_pressureR  = this->m_pressureR;
    retval->m_velxL      = this->m_velxL;
    retval->m_velxR      = this->m_velxR;
    retval->m_velyL      = this->m_velyL;
    retval->m_velzL      = this->m_velzL;
    retval->m_BxL        = this->m_BxL;
    retval->m_ByL        = this->m_ByL;
    retval->m_ByR        = this->m_ByR;
    retval->m_BzL        = this->m_BzL;
    retval->m_BzR        = this->m_BzR;
    retval->m_startX     = this->m_startX;
    retval->m_cloudRho   = this->m_cloudRho;
    retval->m_XC         = this->m_XC;
    retval->m_YC         = this->m_YC;
    retval->m_ZC         = this->m_ZC;
    retval->m_R0         = this->m_R0;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void CloudProblem::fluxBC(       FArrayBox&      a_F,
                                 FArrayBox&      a_Bn,
                           const FArrayBox&      a_WMinus,
                           const FArrayBox&      a_WPlus,
                           const int&            a_dir,
                           const Side::LoHiSide& a_side,
                           const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

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
      
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );

/*
      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-sign);
      a_F.shiftHalf(a_dir,-sign);
      
      Real dx = m_csh->dx(0,m_level);

      // Set the boundary fluxes
      FORT_CLOUDBC( CHF_FRA(a_F),
                    CHF_FRA1(a_Bn,0),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox) );

      // Shift returned fluxes to be face centered
      a_F.shiftHalf(a_dir,sign);*/
    }
  }
}

// Set up initial conditions
void CloudProblem::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_CLOUDINIT( CHF_CONST_FRA(U),
                    CHF_CONST_INT(iCP),
                    CHF_CONST_REAL(dx),
                    CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void CloudProblem::fillGhostCells(       FArrayBox&      a_W,
                                   const FArrayBox&      a_U,
                                   const int&            a_dir,
                                   const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
  int iCP = eqSys->correctionPotentialIndex();

                                   // In periodic case, this doesn't do anything
  if( !m_domain.isPeriodic(a_dir) )
  {
    Box WBox = a_W.box();

                         // See if this chops off the high side of the input box
    Box tmp  = WBox;
    tmp     &= m_domain;

    int indW = WBox.bigEnd( a_dir );
    int indD =  tmp.bigEnd( a_dir );

    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                         // Fill the ghost cells
      FORT_CLOUDGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iCP),
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
      FORT_CLOUDGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iCP),
                    CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void CloudProblem::getBCFlags( eBoundaryConditions leftBC,
                               eBoundaryConditions rightBC,
                               eBoundaryConditions bottomBC,
                               eBoundaryConditions topBC,
                               eBoundaryConditions frontBC,
                               eBoundaryConditions behindBC )
{
  leftBC   = BC_Continuous;
  rightBC  = BC_Fixed;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Periodic;
  behindBC = BC_Periodic;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions CloudProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    
    return (a_sd == Side::Lo ) ? BC_Continuous : BC_Fixed;
  } else {
    return BC_Periodic;
  }
}
