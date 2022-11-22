#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "FastShock.H"
#include "FastShockF_F.H"

#include "LoHiCenter.H"

// Null constructor
FastShock::FastShock()
{
  m_isFortranCommonSet = false;
}

// Input parameters
void FastShock::input( ParmParse & parser, int verbosity )
{
  m_VelModuleL   = 1.0;
  m_VelAngleL    = 0.125;
  m_MagModuleL   = 1.0;
  m_MagAngleL    = 0.1;
  m_PressRatioL  = 0.0;
  m_FastMachL    = 0.0;
  m_startX       = 0.5;
  m_gamma        = 1.6666666667;
  m_AlfvenMag    = 0.0;
  m_AlfvenK      = 2;

  parser.query( "gamma",      m_gamma       );
  parser.query( "VelModule",  m_VelModuleL  );
  parser.query( "VelAngle",   m_VelAngleL   );
  parser.query( "MagModule",  m_MagModuleL  );
  parser.query( "MagAngle",   m_MagAngleL   );
  parser.query( "PressRatio", m_PressRatioL );
  parser.query( "FastMach",   m_FastMachL   );
  parser.query( "startX",     m_startX      );
  parser.query( "AlfvenMag",  m_AlfvenMag   );
  parser.query( "AlfvenK",    m_AlfvenK     );

                                                         // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Fast MHD Shock problem input:"  << endl;
    pout() << "gamma      = " << m_gamma       << endl;
    pout() << "VelModule  = " << m_VelModuleL  << endl;
    pout() << "VelAngle   = " << m_VelAngleL   << endl;
    pout() << "MagModule  = " << m_MagModuleL  << endl;
    pout() << "MagAngle   = " << m_MagAngleL   << endl;
    pout() << "PressRatio = " << m_PressRatioL << endl;
    pout() << "FastMach   = " << m_FastMachL   << endl;
    pout() << "startX     = " << m_startX      << endl;
    pout() << "AlfvenMag  = " << m_AlfvenMag   << endl;
    pout() << "AlfvenK    = " << m_AlfvenK     << endl;
  }

  setFortranCommon( m_gamma,
                    m_VelModuleL,  m_VelAngleL,
                    m_MagModuleL,  m_MagAngleL,
                    m_PressRatioL, m_FastMachL, m_startX, m_AlfvenMag, m_AlfvenK );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void FastShock::setFortranCommon( const Real&     a_gamma,
                                  const Real&     a_VelModuleL,  const Real&     a_VelAngleL,
                                  const Real&     a_MagModuleL,  const Real&     a_MagAngleL,
                                  const Real&     a_PressRatioL, const Real&     a_FastMachL,
                                  const Real&     a_startX,      const Real&     a_AlfvenMag,
                                  const Real&     a_AlfvenK                                   )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETFASTSHOCK( CHF_CONST_REAL( a_gamma       ),
                       CHF_CONST_REAL( a_VelModuleL  ),
                       CHF_CONST_REAL( a_VelAngleL   ),
                       CHF_CONST_REAL( a_MagModuleL  ),
                       CHF_CONST_REAL( a_MagAngleL   ),
                       CHF_CONST_REAL( a_PressRatioL ),
                       CHF_CONST_REAL( a_FastMachL   ),
                       CHF_CONST_REAL( a_startX      ),
                       CHF_CONST_REAL( a_AlfvenMag   ),
                       CHF_CONST_REAL( a_AlfvenK     ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void FastShock::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* FastShock::new_PhysProblem()
{
  FastShock* retval = new FastShock();
  
  retval->copy_PhysProblem(this);

  if (m_isFortranCommonSet == true)
  {
    retval->m_gamma        = this->m_gamma;
    retval->m_VelModuleL   = this->m_VelModuleL;
    retval->m_VelAngleL    = this->m_VelAngleL;
    retval->m_MagModuleL   = this->m_MagModuleL;
    retval->m_MagAngleL    = this->m_MagAngleL;
    retval->m_PressRatioL  = this->m_PressRatioL;
    retval->m_FastMachL    = this->m_FastMachL;
    retval->m_startX       = this->m_startX;
    retval->m_AlfvenMag    = this->m_AlfvenMag;
    retval->m_AlfvenK      = this->m_AlfvenK;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void FastShock::fluxBC(       FArrayBox&      a_F,
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
      FORT_FASTSHOCKBC( CHF_FRA(a_F),
                        CHF_FRA1(a_Bn,0),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_INT(sign),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(a_time),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(boundaryBox) );

                                    // Shift returned fluxes to be face centered
      a_F.shiftHalf( a_dir, sign );
    }
  }
}

// Set up initial conditions
void FastShock::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_FASTSHOCKINIT( CHF_CONST_FRA(U),
                        CHF_CONST_REAL(dx),
                        CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void FastShock::fillGhostCells(       FArrayBox&      a_W,
                                const FArrayBox&      a_U,
                                const int&            a_dir,
                                const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

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

      Real dx = m_csh->dx(0,m_level);
                                                         // Fill the ghost cells
      FORT_FASTSHOCKGS( CHF_FRA(a_W),
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

      Real dx = m_csh->dx(0,m_level);
                                                         // Fill the ghost cells
      FORT_FASTSHOCKGS( CHF_FRA(a_W),
                        CHF_CONST_INT(sign),
                        CHF_CONST_INT(a_dir),
                        CHF_CONST_REAL(dx),
                        CHF_CONST_REAL(a_time),
                        CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void FastShock::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions FastShock::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return BC_Fixed;
  case 1  : return BC_Periodic;
  default : return BC_Undefined;
  }
}
