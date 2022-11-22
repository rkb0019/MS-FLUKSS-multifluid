#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "CurrentSheet.H"
#include "CurrentSheetF_F.H"
#include "EqSysMHDMF.H"
#include "LoHiCenter.H"

// Null constructor
CurrentSheet::CurrentSheet()
{
  m_isFortranCommonSet = false;
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void CurrentSheet::setFortranCommon( const Real&     a_gamma,
                                     const Real&     a_CSheetU0,
                                     const Real&     a_CSheetB0,
                                     const Real&     a_CSheetP0 )
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_SETCSHEET( CHF_CONST_REAL( a_gamma    ),
                  CHF_CONST_REAL( a_CSheetU0 ),
                  CHF_CONST_REAL( a_CSheetB0 ),
                  CHF_CONST_REAL( a_CSheetP0 ) );

  m_isFortranCommonSet = true;
}

// Input parameters
void CurrentSheet::input( ParmParse & parser, int verbosity )
{
  m_CSheetU0  = 0.1;
  m_CSheetB0  = 3.5449;
  m_CSheetP0  = 0.1;

  m_gamma    = 1.6666666667;

  parser.query( "gamma",    m_gamma    );
  parser.query( "CSheetU0", m_CSheetU0 );
  parser.query( "CSheetB0", m_CSheetB0 );
  parser.query( "CSheetP0", m_CSheetP0 );

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Current sheet problem input:" << endl;
    pout() << "gamma     = " << m_gamma      << endl;
    pout() << "CSheetU0  = " << m_CSheetU0   << endl;
    pout() << "CSheetB0  = " << m_CSheetB0   << endl;
    pout() << "CSheetP0  = " << m_CSheetP0   << endl;
  }

  setFortranCommon( m_gamma,
                    m_CSheetU0,
                    m_CSheetB0,
                    m_CSheetP0 );
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void CurrentSheet::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* CurrentSheet::new_PhysProblem()
{
  CurrentSheet* retval = new CurrentSheet();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma    = this->m_gamma;
    retval->m_CSheetU0 = this->m_CSheetU0;
    retval->m_CSheetB0 = this->m_CSheetB0;
    retval->m_CSheetP0 = this->m_CSheetP0;

    retval->setFortranCommonSet();
  }
  

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void CurrentSheet::fluxBC( FArrayBox&            a_F,
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
void CurrentSheet::initialize(LevelData<FArrayBox>& a_U)
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
    FORT_CSHEETINIT( CHF_CONST_FRA(U),
                     CHF_CONST_INT(iCP),
                     CHF_CONST_REAL(dx),
                     CHF_BOX(uBox));
  }
}

                                                             // Fill ghost cells
void CurrentSheet::fillGhostCells(       FArrayBox&      a_W,
                                   const FArrayBox&      a_U,
                                   const int&            a_dir,
                                   const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  // In periodic case, this doesn't do anything

}

//                            Return boundary condition flags for all boundaries
void CurrentSheet::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions CurrentSheet::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  return BC_Periodic;
}
