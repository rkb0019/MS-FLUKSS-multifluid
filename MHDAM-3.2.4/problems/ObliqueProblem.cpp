#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "ObliqueProblem.H"
#include "ObliqueF_F.H"

#include "LoHiCenter.H"

// Null constructor
ObliqueProblem::ObliqueProblem()
{
  m_isFortranCommonSet = false;
}

// Input parameters
void ObliqueProblem::input( ParmParse & parser, int verbosity )
{
  m_densityL   = 1.0;
  m_velxL      = 1.0;
  m_velyL      = 0.0;
  m_angleSW    = 29.0;
  m_startX     = 0.5;
  m_gamma      = 1.4;
  m_pressureL  = 1.0/m_gamma;
  m_obliqueAe  = 0.0;
  m_obliqueAv  = 0.0;
  m_obliquePsi = 0.0;
  m_obliqueK   = 1.0;

  parser.query( "gamma",      m_gamma      );
  parser.query( "densityL",   m_densityL   );
  parser.query( "pressureL",  m_pressureL  );
  parser.query( "velxL",      m_velxL      );
  parser.query( "velyL",      m_velyL      );
  parser.query( "angleSW",    m_angleSW    );
  parser.query( "startX",     m_startX     );
  parser.query( "obliqueAe",  m_obliqueAe  );
  parser.query( "obliqueAv",  m_obliqueAv  );
  parser.query( "obliquePsi", m_obliquePsi );
  parser.query( "obliqueK",   m_obliqueK   );

                                                         // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "1D Riemann problem input:"     << endl;
    pout() << "gamma      = " << m_gamma      << endl;
    pout() << "densityL   = " << m_densityL   << endl;
    pout() << "pressureL  = " << m_pressureL  << endl;
    pout() << "velxL      = " << m_velxL      << endl;
    pout() << "velyL      = " << m_velyL      << endl;
    pout() << "angleSW    = " << m_angleSW    << endl;
    pout() << "startX     = " << m_startX     << endl;
    pout() << "obliqueAe  = " << m_obliqueAe  << endl;
    pout() << "obliqueAv  = " << m_obliqueAv  << endl;
    pout() << "obliquePsi = " << m_obliquePsi << endl;
    pout() << "obliqueK   = " << m_obliqueK   << endl;
  }

  setFortranCommon( m_gamma,
                    m_densityL,     m_pressureL,
                    m_velxL,        m_velyL,
                    m_angleSW,      m_startX,
                    m_obliqueAe,    m_obliqueAv,
                    m_obliquePsi,   m_obliqueK );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void ObliqueProblem::setFortranCommon( const Real&     a_gamma,
                                       const Real&     a_densityL,
                                       const Real&     a_pressureL,
                                       const Real&     a_velxL,
                                       const Real&     a_velyL,
                                       const Real&     a_angleSW,
                                       const Real&     a_startX,
                                       const Real&     a_obliqueAe,
                                       const Real&     a_obliqueAv,
                                       const Real&     a_obliquePsi,
                                       const Real&     a_obliqueK   )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETOBLIQUE( CHF_CONST_REAL( a_gamma      ),
                     CHF_CONST_REAL( a_densityL   ),
                     CHF_CONST_REAL( a_pressureL  ),
                     CHF_CONST_REAL( a_velxL      ),
                     CHF_CONST_REAL( a_velyL      ),
                     CHF_CONST_REAL( a_angleSW    ),
                     CHF_CONST_REAL( a_startX     ),
                     CHF_CONST_REAL( a_obliqueAe  ),
                     CHF_CONST_REAL( a_obliqueAv  ),
                     CHF_CONST_REAL( a_obliquePsi ),
                     CHF_CONST_REAL( a_obliqueK   ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void ObliqueProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* ObliqueProblem::new_PhysProblem()
{
  ObliqueProblem* retval = new ObliqueProblem();
  
  retval->copy_PhysProblem(this);

  if (m_isFortranCommonSet == true)
  {
    retval->m_gamma      = this->m_gamma;
    retval->m_densityL   = this->m_densityL;
    retval->m_pressureL  = this->m_pressureL;
    retval->m_velxL      = this->m_velxL;
    retval->m_velyL      = this->m_velyL;
    retval->m_angleSW    = this->m_angleSW;
    retval->m_startX     = this->m_startX;
    retval->m_obliqueAe  = this->m_obliqueAe;
    retval->m_obliqueAv  = this->m_obliqueAv;
    retval->m_obliquePsi = this->m_obliquePsi;
    retval->m_obliqueK   = this->m_obliqueK;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void ObliqueProblem::fluxBC(       FArrayBox&      a_F,
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
      
      Real dx = m_csh->dx(0,m_level);
                                                      // Set the boundary fluxes
      FORT_FLUXBCOS( CHF_FRA(a_F),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(sign),
                     CHF_CONST_REAL(dx),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(boundaryBox) );

                                    // Shift returned fluxes to be face centered
      a_F.shiftHalf( a_dir, sign );
    }
  }
}

// Set up initial conditions
void ObliqueProblem::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_OBLIQUEINIT( CHF_CONST_FRA(U),
                      CHF_CONST_REAL(dx),
                      CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void ObliqueProblem::fillGhostCells(       FArrayBox&      a_W,
                                     const FArrayBox&      a_U,
                                     const int&            a_dir,
                                     const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);

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
      FORT_OBLIQUEGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(a_time),
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
      FORT_OBLIQUEGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(a_time),
                      CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void ObliqueProblem::getBCFlags( eBoundaryConditions leftBC,
                                 eBoundaryConditions rightBC,
                                 eBoundaryConditions bottomBC,
                                 eBoundaryConditions topBC,
                                 eBoundaryConditions frontBC,
                                 eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Continuous;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Periodic;
  behindBC = BC_Periodic;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions ObliqueProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  } else {
    return BC_Periodic;
  }
}
