#include "Roe8WavesRS.H"

#include "RiemannSolverF_F.H"
#include "Roe8WaveF_F.H"

Roe8WavesRS::Roe8WavesRS( void )
 : m_dEntropyFix(0.1), m_dSmallB(1.0e-12), m_iAveraging(0), m_iLaxFriedrix(0)
{
  FORT_SETRSCONST( CHF_CONST_REAL( m_dEntropyFix  ),
                   CHF_CONST_REAL( m_dSmallB      ),
                   CHF_CONST_INT ( m_iAveraging   ),
                   CHF_CONST_INT ( m_iLaxFriedrix ) );
}

RiemannSolver * Roe8WavesRS::new_RiemannSolver()
{
  Roe8WavesRS * retval = new Roe8WavesRS();

  retval->setParameters( m_dEntropyFix, m_dSmallB, m_iAveraging, m_iLaxFriedrix );

  return static_cast<RiemannSolver*>(retval);
}

void Roe8WavesRS::setParameters( double EntFix, double dSmallB,
                                 int    iAver,  int    iLF )
{
  m_dEntropyFix  = EntFix;
  m_dSmallB      = dSmallB;
  m_iAveraging   = iAver;
  m_iLaxFriedrix = iLF;

  FORT_SETRSCONST( CHF_CONST_REAL( m_dEntropyFix  ),
                   CHF_CONST_REAL( m_dSmallB      ),
                   CHF_CONST_INT ( m_iAveraging   ),
                   CHF_CONST_INT ( m_iLaxFriedrix ) );
}

void Roe8WavesRS::fluxes(       FArrayBox & a_F,
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
  FORT_RIEMANNF(CHF_FRA(a_F),
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
