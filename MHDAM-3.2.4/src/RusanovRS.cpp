#include "RusanovRS.H"

#include "RiemannSolverF_F.H"
#include "RusanovF_F.H"

RusanovRS::RusanovRS( void )
 : m_dSmallB(1.0e-12), m_iAveraging(0)
{
  Real dEntropyFix   = 0.0;
  int  iLaxFriedrix  = 1;
  FORT_SETRSCONST( CHF_CONST_REAL( dEntropyFix  ),
                   CHF_CONST_REAL( m_dSmallB    ),
                   CHF_CONST_INT ( m_iAveraging ),
                   CHF_CONST_INT ( iLaxFriedrix ) );
}

RiemannSolver * RusanovRS::new_RiemannSolver()
{
  RusanovRS * retval = new RusanovRS();

  retval->setParameters( m_dSmallB, m_iAveraging );

  return static_cast<RiemannSolver*>(retval);
}

void RusanovRS::setParameters( double dSmallB, int iAver )
{
  m_dSmallB      = dSmallB;
  m_iAveraging   = iAver;

  Real dEntropyFix   = 0.0;
  int  iLaxFriedrix  = 1;

  FORT_SETRSCONST( CHF_CONST_REAL( dEntropyFix  ),
                   CHF_CONST_REAL( m_dSmallB    ),
                   CHF_CONST_INT ( m_iAveraging ),
                   CHF_CONST_INT ( iLaxFriedrix ) );
}

void RusanovRS::fluxes(       FArrayBox & a_F,
                        const FArrayBox & a_WPlus,
                        const FArrayBox & a_WMinus,
                        const int &       a_dir,
                        const int &       a_iRho,
                        const Box &       a_box )
{
  CH_assert( a_F.box().contains(a_box) );

                                        // Get the numbers of relevant variables
  CH_assert(a_F     .nComp() >= a_iRho + 8);
  CH_assert(a_WPlus .nComp() >= a_iRho + 8);
  CH_assert(a_WMinus.nComp() >= a_iRho + 8);

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
  FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);
  
  CH_assert( shiftWLeft.box().contains(a_box) );
  CH_assert( shiftWRight.box().contains(a_box) );

                                   // Riemann solver computes a_F for all edges
                                      // that are not on the physical boundary.
  FORT_RUSANOVF( CHF_FRA(a_F),
                 CHF_CONST_FRA(shiftWLeft),
                 CHF_CONST_FRA(shiftWRight),
                 CHF_CONST_INT(a_dir),
                 CHF_CONST_INT(a_iRho),
                 CHF_BOX(a_box));
 
                            // Shift the left and right primitive variable boxes
                            // back to their original position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}
