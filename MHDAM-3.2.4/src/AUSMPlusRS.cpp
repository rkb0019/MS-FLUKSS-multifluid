#include "AUSMPlusRS.H"

#include "AUSMPlusRSF_F.H"

AUSMPlusRS::AUSMPlusRS( void )
{
}

RiemannSolver * AUSMPlusRS::new_RiemannSolver()
{
  AUSMPlusRS * retval = new AUSMPlusRS();

  return static_cast<RiemannSolver*>(retval);
}

void AUSMPlusRS::fluxes(       FArrayBox & a_F,
                         const FArrayBox & a_WLeft,
                         const FArrayBox & a_WRight,
                         const int &       a_dir,
                         const int &       a_iRho,
                         const Box &       a_box )
{
  CH_assert( a_F.box().contains(a_box) );

                                        // Get the numbers of relevant variables
  CH_assert(a_F     .nComp() >= a_iRho + 5);
  CH_assert(a_WLeft .nComp() >= a_iRho + 5);
  CH_assert(a_WRight.nComp() >= a_iRho + 5);

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WLeft;
  FArrayBox& shiftWRight = (FArrayBox&)a_WRight;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);
  
  CH_assert( shiftWLeft.box().contains(a_box) );
  CH_assert( shiftWRight.box().contains(a_box) );

                                   // Riemann solver computes a_F for all edges
                                      // that are not on the physical boundary.
  FORT_AUSMPLUSF(CHF_FRA(a_F),
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
