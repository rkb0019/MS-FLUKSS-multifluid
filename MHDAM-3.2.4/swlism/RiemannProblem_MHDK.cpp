#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "RiemannProblem_MHDK.H"
#include "RiemannProblem_MHDKF_F.H"

#include "LoHiCenter.H"

// Null constructor
RiemannProblem_MHDK::RiemannProblem_MHDK()
{
  m_isFortranCommonSet = false;
}

// Input parameters
void RiemannProblem_MHDK::input( ParmParse & parser, int verbosity )
{
  m_numdenL      = 1.0;
  m_temperatureL = 5000.0;
  m_velxL        = 25000.0;
  m_velyL        = 25000.0;
  m_BxL          = 0.0;
  m_ByL          = 0.0;

  m_numdenR      = 1.0;
  m_temperatureR = 5000.0;
  m_velxR        = 25000.0;
  m_velyR        = 25000.0;
  m_BxR          = 0.0;
  m_ByR          = 0.0;

  m_startX       = 0.5;
  m_netnumL      = 1.0;
  m_gamma        = 1.6666666667;

  parser.query( "gamma",        m_gamma        );
  parser.query( "netnumL",      m_netnumL      );
  parser.query( "startX",       m_startX       );

  parser.query( "numdenL",      m_numdenL      );
  parser.query( "temperatureL", m_temperatureL );
  parser.query( "velxL",        m_velxL        );
  parser.query( "velyL",        m_velyL        );
  parser.query( "BxL",          m_BxL          );
  parser.query( "ByL",          m_ByL          );

  parser.query( "numdenR",      m_numdenR      );
  parser.query( "temperatureR", m_temperatureR );
  parser.query( "velxR",        m_velxR        );
  parser.query( "velyR",        m_velyR        );
  parser.query( "BxR",          m_BxR          );
  parser.query( "ByR",          m_ByR          );

                                                         // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "1D Riemann problem with kinetic source terms:"      << endl;
    pout() << "gamma        = " << m_gamma        << endl;
    pout() << "netnumL      = " << m_netnumL      << endl;
    pout() << "startX       = " << m_startX       << endl << endl;

    pout() << "numdenL      = " << m_numdenL      << endl;
    pout() << "temperatureL = " << m_temperatureL << endl;
    pout() << "velxL        = " << m_velxL        << endl;
    pout() << "velyL        = " << m_velyL        << endl;
    pout() << "BxL          = " << m_BxL          << endl;
    pout() << "ByL          = " << m_ByL          << endl;

    pout() << "numdenR      = " << m_numdenR      << endl;
    pout() << "temperatureR = " << m_temperatureR << endl;
    pout() << "velxR        = " << m_velxR        << endl;
    pout() << "velyR        = " << m_velyR        << endl;
    pout() << "BxR          = " << m_BxR          << endl;
    pout() << "ByR          = " << m_ByR          << endl;
  }

  setFortranCommon( m_gamma,   m_startX,       m_netnumL,
                    m_numdenL, m_temperatureL, m_velxL, m_velyL, m_BxL, m_ByL,
                    m_numdenR, m_temperatureR, m_velxR, m_velyR, m_BxR, m_ByR );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RiemannProblem_MHDK::setFortranCommon( const Real& a_gamma,   const Real& a_startX, const Real& a_netnumL,
                         const Real& a_numdenL, const Real& a_temperatureL,
                         const Real& a_velxL,   const Real& a_velyL,
                         const Real& a_BxL,     const Real& a_ByL,
                         const Real& a_numdenR, const Real& a_temperatureR,
                         const Real& a_velxR,   const Real& a_velyR,
                         const Real& a_BxR,     const Real& a_ByR     )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETRIEMANN_MHDK(
      CHF_CONST_REAL(a_gamma),   CHF_CONST_REAL(a_startX), CHF_CONST_REAL(a_netnumL),
      CHF_CONST_REAL(a_numdenL), CHF_CONST_REAL(a_temperatureL),
      CHF_CONST_REAL(a_velxL),   CHF_CONST_REAL(a_velyL), CHF_CONST_REAL(a_BxL), CHF_CONST_REAL(a_ByL),
      CHF_CONST_REAL(a_numdenR), CHF_CONST_REAL(a_temperatureR),
      CHF_CONST_REAL(a_velxR),   CHF_CONST_REAL(a_velyR), CHF_CONST_REAL(a_BxR), CHF_CONST_REAL(a_ByR));

  //Real scaleLen  = 1.5e+13;
  //Real scaleNDen = a_numdenL;

  //FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( a_velxL   ),
  //                       CHF_CONST_REAL( scaleLen  ),
  //                       CHF_CONST_REAL( scaleNDen ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RiemannProblem_MHDK::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* RiemannProblem_MHDK::new_PhysProblem()
{
  RiemannProblem_MHDK* retval = new RiemannProblem_MHDK();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma        = this->m_gamma;
    retval->m_numdenL      = this->m_numdenL;
    retval->m_temperatureL = this->m_temperatureL;
    retval->m_velxL        = this->m_velxL;
    retval->m_velyL        = this->m_velyL;
    retval->m_BxL          = this->m_BxL;
    retval->m_ByL          = this->m_ByL;

    retval->m_numdenR      = this->m_numdenR;
    retval->m_temperatureR = this->m_temperatureR;
    retval->m_velxR        = this->m_velxR;
    retval->m_velyR        = this->m_velyR;
    retval->m_BxR          = this->m_BxR;
    retval->m_ByR          = this->m_ByR;

    retval->m_startX       = this->m_startX;
    retval->m_netnumL      = this->m_netnumL;

    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void RiemannProblem_MHDK::fluxBC(     FArrayBox&      a_F,
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
      boundaryBox.shiftHalf(a_dir,-sign);
      a_F.shiftHalf(a_dir,-sign);

      // Set the boundary fluxes
      FORT_FLUXBC_MHDK( CHF_FRA(a_F),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(boundaryBox) );

      // Shift returned fluxes to be face centered
      a_F.shiftHalf(a_dir,sign);
    }
  }
}

// Set up initial conditions
void RiemannProblem_MHDK::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_RIEMANNINIT_MHDK( CHF_CONST_FRA(U),
                         CHF_CONST_REAL(dx),
                         CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void RiemannProblem_MHDK::fillGhostCells(       FArrayBox&      a_W,
                                          const FArrayBox&      a_U,
                                          const int&            a_dir,
                                          const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
}

//                            Return boundary condition flags for all boundaries
void RiemannProblem_MHDK::getBCFlags( eBoundaryConditions leftBC,
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
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions RiemannProblem_MHDK::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  case 1  : return BC_Periodic;
  default : return BC_Undefined;
  }
}
