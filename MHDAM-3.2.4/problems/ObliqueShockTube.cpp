#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "ObliqueShockTube.H"
#include "ObliqueShockTubeF_F.H"
#include "EqSysMHDMF.H"
#include "MHDAMDefs.H"
#include "LoHiCenter.H"
#include "LGintegrator.H"
#include "RiemannSolver.H"
#include "PatchIdealMHDF_F.H"

// Null constructor
ObliqueShockTube::ObliqueShockTube()
{
    m_isFortranCommonSet = false;
}

// Input parameters
void ObliqueShockTube::input( ParmParse & parser, int verbosity )
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
  m_startY     = 0.5;
  m_tanangle      = 1;
  
  m_gamma      = 1.6666666667;
  
  Real VParL = 0.0,VParR = 0.0,VPerL = 0.0,VPerR = 0.0;
  Real BParL = 0.0,BParR = 0.0,BPerL = 0.0,BPerR = 0.0;

  parser.query( "gamma",     m_gamma     );
  parser.query( "densityL",  m_densityL  );
  parser.query( "densityR",  m_densityR  );
  parser.query( "pressureL", m_pressureL );
  parser.query( "pressureR", m_pressureR );
  parser.query( "velxL",   VParL     );
  parser.query( "velxR",   VParR     );
  parser.query( "velyL",   VPerL     );
  parser.query( "velyR",   VPerR     );
  parser.query( "velzL",     m_velzL     );
  parser.query( "velzR",     m_velzR     );
  parser.query( "BxL",     BParL       );
  parser.query( "BxR",     BParR       );
  parser.query( "ByL",     BPerL       );
  parser.query( "ByR",     BPerR       );
  parser.query( "BzL",       m_BzL       );
  parser.query( "BzR",       m_BzR       );
  parser.query( "startX",    m_startX    );  
  parser.query( "tanangle",     m_tanangle    );
  
  if ((m_tanangle<=0) || (m_tanangle>4))
  {
    MayDay::Error("tanangle should be 1/n, 0, 1, 2, 3 or 4");
  }
    
  
  Real domainLength = 1.0;
  parser.query("domain_length",domainLength);
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim); 
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  parser.queryarr("num_cells",numCells,0,SpaceDim);
      
  Real dx = domainLength/numCells[0];    
  
  m_startY = 0.5*dx*numCells[1]+0.1*dx;
  
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
    pout() << "startY    = " << m_startY    << endl;
    pout() << "tanangle  = " << m_tanangle     << endl;    
  }
  
  Real angle = atan(m_tanangle);
  
  m_velxL      = VParL*sin(angle) + VPerL*cos(angle);
  m_velxR      = VParR*sin(angle) + VPerR*cos(angle);
  m_velyL      = VParL*cos(angle) - VPerL*sin(angle);
  m_velyR      = VParR*cos(angle) - VPerR*sin(angle);
  
  m_BxL      = BParL*sin(angle) + BPerL*cos(angle);
  m_BxR      = BParR*sin(angle) + BPerR*cos(angle);
  m_ByL      = BParL*cos(angle) - BPerL*sin(angle);
  m_ByR      = BParR*cos(angle) - BPerR*sin(angle);


  setFortranCommon( );
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void ObliqueShockTube::setFortranCommon( )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETOBLIQUESHOCKTUBE( CHF_CONST_REAL( m_gamma     ),
                     CHF_CONST_REAL( m_densityL  ),
                     CHF_CONST_REAL( m_densityR  ),
                     CHF_CONST_REAL( m_pressureL ),
                     CHF_CONST_REAL( m_pressureR ),
                     CHF_CONST_REAL( m_velxL     ),
                     CHF_CONST_REAL( m_velxR     ),
                     CHF_CONST_REAL( m_velyL     ),
                     CHF_CONST_REAL( m_velyR     ),
                     CHF_CONST_REAL( m_velzL     ),
                     CHF_CONST_REAL( m_velzR     ),
                     CHF_CONST_REAL( m_BxL       ),
                     CHF_CONST_REAL( m_BxR       ),
                     CHF_CONST_REAL( m_ByL       ),
                     CHF_CONST_REAL( m_ByR       ),
                     CHF_CONST_REAL( m_BzL       ),
                     CHF_CONST_REAL( m_BzR       ),
                     CHF_CONST_REAL( m_startX    ),
                     CHF_CONST_REAL( m_startY    ) );

    m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void ObliqueShockTube::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* ObliqueShockTube::new_PhysProblem()
{
  ObliqueShockTube* retval = new ObliqueShockTube();

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
    retval->m_startY     = this->m_startY;
    retval->m_tanangle   = this->m_tanangle;    

    retval->setFortranCommonSet();
  }
  
  retval->copy_PhysProblem(this);

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void ObliqueShockTube::fluxBC(       FArrayBox&      a_F,
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
      
      
            // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
    FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);

                                        // Computes a_Bn for all edges
                                        // that are not on the physical boundary
    CH_assert( shiftWLeft.box().contains(boundaryBox) );
    CH_assert( shiftWRight.box().contains(boundaryBox) );

    FORT_COMPUTEBN( CHF_FRA1(a_Bn,0),
                    CHF_CONST_FRA(shiftWLeft),
                    CHF_CONST_FRA(shiftWRight),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(boundaryBox));
                    
    // Shift the left and right primitive variable boxes
    // back to their original position
    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);


    }
  }
}

// Set up initial conditions
void ObliqueShockTube::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = m_csh->dx(0,m_level);
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
  int iCP = eqSys->correctionPotentialIndex();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition in this grid
    FORT_OBLIQUESHOCKTUBEINIT(CHF_CONST_FRA(U),
                     CHF_CONST_REAL(m_tanangle),
                     CHF_CONST_INT(iCP),
                     CHF_CONST_REAL(dx),
                     CHF_BOX(uBox));
  }
}
                                                             // Fill ghost cells
void ObliqueShockTube::fillGhostCells(       FArrayBox&      a_W,
                                     const FArrayBox&      a_U,
                                     const int&            a_dir,
                                     const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);
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
      FORT_OBLIQUESHOCKTUBEGS( CHF_FRA(a_W),
                      CHF_CONST_REAL(m_tanangle),
                      CHF_CONST_INT(iCP),
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
      FORT_OBLIQUESHOCKTUBEGS( CHF_FRA(a_W),
                      CHF_CONST_REAL(m_tanangle),
                      CHF_CONST_INT(iCP),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),                      
                      CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void ObliqueShockTube::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions ObliqueShockTube::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return BC_Fixed;
  case 1  : return BC_Periodic;
  default : return BC_Undefined;
  }
}

void ObliqueShockTube::primForPlot(      FArrayBox& a_W,
                           const Box&       a_box)
{
  //return;
  
  Real vx,vy,vper,vpar,sinA = sin(atan(m_tanangle)),cosA = cos(atan(m_tanangle));
  
  BoxIterator bit(a_box);
  
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    
    vx = a_W(iv,WVELX);
    vy = a_W(iv,WVELY);    
    vpar = vx*sinA+vy*cosA;
    vper = vx*cosA-vy*sinA;    
    a_W(iv,WVELX) = vpar;
    a_W(iv,WVELY) = vper;
    
    vx = a_W(iv,WBX);
    vy = a_W(iv,WBY);
    vpar = vx*sinA+vy*cosA;
    vper = vx*cosA-vy*sinA;    
    a_W(iv,WBX) = vpar;
    a_W(iv,WBY) = vper;    
  }  
  
}

/// Number additional variables for writing to plot file
int ObliqueShockTube::numPlotVars()
{
  return 2;
}
  
/// Names of the additional variables for writing to plot file
Vector<std::string> ObliqueShockTube::plotNames()
{
  Vector<string> retval;  

  retval.push_back("BxError");
  retval.push_back("ByError");
  
  return retval;
}
  
/// Calculates variables for plotting using primitive variables
void ObliqueShockTube::calcPlotVars(FArrayBox&      a_vars,
                           int              a_comp,
                           const FArrayBox& a_W,
                           const Box&       a_box)
{
  Real dx = m_csh->dx(0,m_level);
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
  int iCP = eqSys->correctionPotentialIndex();
  
  FArrayBox Uexact(a_box, m_eqSys->numStates());
  FArrayBox Wexact(a_box, m_eqSys->numPrimitives());
      
  FORT_OBLIQUESHOCKTUBEINIT(CHF_CONST_FRA(Uexact),
                     CHF_CONST_REAL(m_tanangle),
                     CHF_CONST_INT(iCP),
                     CHF_CONST_REAL(dx),
                     CHF_BOX(a_box));
                     
  m_eqSys->stateToPrim(Wexact, Uexact, a_box);
  
  Real vx,vy,vper,vpar,sinA = sin(atan(m_tanangle)),cosA = cos(atan(m_tanangle));
  
  BoxIterator bit(a_box);
  
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
            
    vx = a_W(iv,WBX)-Wexact(iv,WBX);
    vy = a_W(iv,WBY)-Wexact(iv,WBY);
    vpar = vx*sinA+vy*cosA;
    vper = vx*cosA-vy*sinA;    
    a_vars(iv,a_comp)   = fabs(vpar);
    a_vars(iv,a_comp+1) = fabs(vper);    
  }  
  
}
