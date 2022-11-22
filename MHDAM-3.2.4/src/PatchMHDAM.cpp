#include <cstdio>
#include <string>
#include <limits>
#include "DebugOut.H"

#include "PatchMHDAM.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"

#include "printArrayF_F.H"
#include "DivCleaningF_F.H"

#include "LGintegrator.H"
#include "PatchMHDAMF_F.H"
#include "LevelSetF_F.H"
#include "EquationSystem.H"
#include "MHDAMDefs.H"

#include "PatchIdealMHDF_F.H" // temporarily!!!!

//                                         Flag everything as not defined or set
PatchMHDAM::PatchMHDAM()
{
  m_isDefined          = false;
  m_isPPSet            = false;
  m_isRecSet           = false;
  m_isCurrentTimeSet   = false;  
  m_isRSSet            = false;
  m_isRSSetGD          = false;
  m_bLSonly            = false;

  m_PhPr               = NULL;
  m_RS                 = NULL;
  m_RSGD               = NULL;
  m_Reconstruction     = NULL;

  m_TimeMethod         = TA_FirstOrder;
  m_iDivBMethod        = 0;
  m_level              = -1;
  
  m_csh                = NULL;
  m_eqSys              = NULL;
    
}

PatchMHDAM::~PatchMHDAM()
{
                  // Delete the initial/boundary condition object - if it exists
  if( m_PhPr != NULL )
  {
    delete m_PhPr;
  }
                              // Delete the Riemann solver object - if it exists
  if( m_RS != NULL )
  {
    delete m_RS;
  }
  if( m_RSGD != NULL )
  {
    delete m_RSGD;
  }
}

void PatchMHDAM::input( ParmParse & parser, int verbosity )
{
}

//                          Define this object and the boundary condition object
void PatchMHDAM::define(       ProblemDomain& a_domain,                         
                         const int            a_level,
                         const int            a_verbosity)
{
  CH_assert(m_isPPSet);

  // Store the domain and grid spacing
  m_domain = a_domain;  
  m_level  = a_level;

  // Set the domain and grid spacing in the physical problem object
  m_PhPr->setNumGhostCells(numGhostCells());
  m_PhPr->define(m_domain,a_level);
  
  m_verbosity = a_verbosity;

  m_isDefined = true;
}

//                                      Copy internal data to another PatchMHDAM
void PatchMHDAM::copyTo( PatchMHDAM * pPatch ) const
{
  pPatch->setPhysProblem( getPhysProblem() );

  pPatch->setRiemannSolver( getRiemannSolver() );

  if( m_isRSSetGD == 1 )
  {
    pPatch->setRiemannSolverGD( getRiemannSolverGD() );
  }

  pPatch->setReconstruction( getReconstruction() );

  pPatch->setTimeApproximation( timeApproximation() );

  pPatch->setDivBMethod( getDivBMethod() );
  pPatch->setDivergenceCleaning( getDivergenceCleaning() );

  pPatch->setCurrentTime( m_currentTime );  
  
  pPatch->m_b8waveUse = m_b8waveUse;
  pPatch->m_bLSonly   = m_bLSonly;
    

  for( int i = 0; i < REF_NUM; i++ )
    pPatch->RefineParams[i]  = RefineParams[i];
}

//                                             Set the boundary condition object
void PatchMHDAM::setPhysProblem( PhysProblem* a_PhPr )
{
  // Delete old boundary condition object - if any
  if (m_PhPr != NULL)
  {
    delete m_PhPr;
  }

  // Store new boundary condition object
  m_PhPr    = a_PhPr->new_PhysProblem();
  m_isPPSet = true;
  
  
  m_csh   = m_PhPr->coordinateSystem();
  m_eqSys = m_PhPr->equationSystem();
}

//                                  Return the current boundary condition object
PhysProblem* PatchMHDAM::getPhysProblem() const
{
  CH_assert(m_isPPSet);

  return m_PhPr;
}

//                                             Set the MHD Riemann solver object
void PatchMHDAM::setRiemannSolver( RiemannSolver * a_RS )
{
                                    // Delete old Riemann solver object - if any
  if( m_RS != NULL )
  {
    delete m_RS;
  }

                                          // Store new MHD Riemann solver object
  m_RS      = a_RS->new_RiemannSolver();
  m_isRSSet = true;
}

//                                  Return the current MHD Riemann solver object
RiemannSolver * PatchMHDAM::getRiemannSolver( void ) const
{
  CH_assert(m_isRSSet);

  return m_RS;
}

//                                    Set the Gas Dynamics Riemann solver object
void PatchMHDAM::setRiemannSolverGD( RiemannSolver * a_RS )
{
                                    // Delete old Riemann solver object - if any
  if( m_RSGD != NULL )
  {
    delete m_RSGD;
  }

                                 // Store new Gas Dynamics Riemann solver object
  m_RSGD    = a_RS->new_RiemannSolver();
  m_isRSSetGD = true;
}

//                         Return the current Gas Dynamics Riemann solver object
RiemannSolver * PatchMHDAM::getRiemannSolverGD( void ) const
{
  CH_assert(m_isRSSetGD);

  return m_RSGD;
}

//                                                 Set the Reconstruction object
void PatchMHDAM::setReconstruction( Reconstruction * a_Rec )
{
                                    // Delete old Reconstruction object - if any
  if( m_Reconstruction != NULL )
  {
    delete m_Reconstruction;
  }

                                              // Store new Reconstruction object
  m_Reconstruction = a_Rec->new_Reconstruction();
  m_isRecSet       = true;
}

//                                      Return the current Reconstruction object
Reconstruction * PatchMHDAM::getReconstruction( void ) const
{
  CH_assert(m_isRecSet);

  return m_Reconstruction;
}

//                                             Get the source calculator pointer
SourceCalculator * PatchMHDAM::getSourceCalculator( void ) const
{
  if( m_isPPSet == false ) return NULL;
  if( m_PhPr    == NULL  ) return NULL;

  return m_PhPr->getSourceCalculator();
}

                                              // Set the Time Approximation flag
void PatchMHDAM::setTimeApproximation( PatchMHDAM::eTimeApproximation method )
{
  m_TimeMethod  = method;
}

                                              // Get the Time Approximation flag
PatchMHDAM::eTimeApproximation PatchMHDAM::timeApproximation( void ) const
{
  return m_TimeMethod;
}

                                            // Set the Divergence Cleaning flag
void PatchMHDAM::setDivergenceCleaning( eDivergenceCleaning method )
{
  m_DivergenceCleaning = method;
}

                                             // Get the Divergence Cleaning flag
eDivergenceCleaning PatchMHDAM::getDivergenceCleaning( void ) const
{
  return m_DivergenceCleaning;
}

                                       // Set the div(B) calculation method flag
void PatchMHDAM::setDivBMethod( int method )
{
  m_iDivBMethod  = method;
}

                                       // Get the div(B) calculation method flag
int PatchMHDAM::getDivBMethod( void ) const
{
  return m_iDivBMethod;
}


void PatchMHDAM::set8waveFlag( int a_8waveUse )
{
  m_b8waveUse = (a_8waveUse == 1);
}

// Set m_bLSonly flag  
void PatchMHDAM::setLSonlyFlag( int a_LSonly )
{
  m_bLSonly = a_LSonly;
}

bool PatchMHDAM::getLSonlyFlag( void ) const
{
  return m_bLSonly;
}

//   Set the current physical time - used for time dependent boundary conditions
void PatchMHDAM::setCurrentTime(const Real& a_currentTime)
{
  m_currentTime = a_currentTime;
  m_isCurrentTimeSet = true;
}


// Determing the number of ghost cells necessary here
int PatchMHDAM::numGhostCells()
{
  if (m_TimeMethod == TA_RK2)     return 4;
  if (m_TimeMethod == TA_Hancock) return 3;
  
  return 2;    
}

void PatchMHDAM::computeBn(       FArrayBox & a_Bn,
                            const FArrayBox & a_W,
                            const FArrayBox & a_WMinus,
                            const FArrayBox & a_WPlus ,
                            const int       & a_method,
                            const int       & a_dir,
                            const Box       & a_box)
{

}

void PatchMHDAM::computeUn(       FArrayBox & a_Un,
                            const FArrayBox & a_W,
                            const FArrayBox & a_WMinus,
                            const FArrayBox & a_WPlus ,
                            const int       & a_method,
                            const int       & a_dir,
                            const Box       & a_box)
{

}

//                                       Compute dt and returns a cell with the minimum dt. 
Real PatchMHDAM::computeDtLevelSet( const FArrayBox& a_U,
                                 const Box&     a_box,
                                 IntVect&       a_minDtCell)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  
  FArrayBox dtFab(a_box,1);  
  dtFab.setVal(numeric_limits<Real>::max());
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
  {    
    int ivf = m_PhPr->lsIndexField(i);
            
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))  
    {
      Real dx = m_csh->dx(0,m_level);
      
      FORT_MINDT_LSCARTESIAN(
        CHF_FRA1(dtFab,0),
        CHF_CONST_FRA(a_U),               
        CHF_CONST_INT(ivf), 
        CHF_CONST_REAL(dx),    
        CHF_BOX(a_box));
    }
    
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Spherical) )
    {
      FArrayBox USph(a_box,a_U.nComp());
      USph.copy(a_U);
        
      m_csh->transCartesianVectToCurv(USph,a_box,m_level);    
      FORT_MINDT_LSSPHERICAL(
        CHF_FRA1(dtFab,0),
        CHF_CONST_FRA(USph),               
        CHF_CONST_INT(ivf), 
        CHF_CONST_INT(m_level),    
        CHF_BOX(a_box));
    }    
  }    
    
  Real dt = m_PhPr->computeDt(a_U, dtFab, a_box, a_minDtCell);
  
  return dt;  
}

void PatchMHDAM::levelSetFluxes(       FArrayBox & a_F,
                         const FArrayBox & a_WPlus,
                         const FArrayBox & a_WMinus,
                         const int &       a_dir,
                         const Box &       a_box)
{
  
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

                                                                   
  
  for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
  {
    int ils  = m_eqSys->lsIndexCons(i);
    int ivf = m_PhPr->lsIndexField(i);
    
    FORT_LEVELSETFLUX(CHF_FRA(a_F),                
                CHF_CONST_FRA(shiftWLeft),
                CHF_CONST_FRA(shiftWRight),
                CHF_CONST_INT(a_dir),
                CHF_CONST_INT(ils),
                CHF_CONST_INT(ivf),
                CHF_BOX(a_box));        
  }  
  
                            // Shift the left and right primitive variable boxes
                            // back to their original position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
                             
}                         
  

              // Update the conserved variable using fluxes and a scaling factor
void PatchMHDAM::updateCons(       FArrayBox& a_U,
                                  const FArrayBox& a_F,
                                        Real       a_dt,
                                  const FArrayBox& a_scale,
                                  const int&       a_dir,
                                  const Box&       a_box)
{
  CH_assert(isDefined());      
  
  int iBgn = m_eqSys->consInterval().begin();
  int iEnd = m_eqSys->consInterval().end();
        
  FORT_UPDATECONS(CHF_FRA(a_U),
                  CHF_CONST_FRA(a_F),
                  CHF_CONST_INT(iBgn),
                  CHF_CONST_INT(iEnd),
                  CHF_CONST_REAL(a_dt),
                  CHF_CONST_FRA1(a_scale,0),
                  CHF_CONST_INT(a_dir),
                  CHF_BOX(a_box));
}

              // Update the conserved variable using fluxes and a scaling factor
void PatchMHDAM::updateLevelSet(        FArrayBox& a_U,
                                  const FArrayBox& a_W,
                                  const FArrayBox& a_F,
                                        Real       a_dt,
                                  const int&       a_dir,
                                  const Box&       a_box)
{
  CH_assert(isDefined());      
  
  if (m_eqSys->numTrackingSurfaces()==0) return;
  
  FArrayBox dx,g(a_box,1);
  
  IntVect iv_off(IntVect::Zero);
  iv_off[a_dir]=1;

  Box dxBox(a_box.smallEnd()*iv_off, 
            a_box.bigEnd()  *iv_off);            

  dx.define(dxBox,1);                     
  m_csh->dx(dx,dxBox,a_dir,m_level);          
      
  m_csh->sqrtMetricTensorCoeff(g, a_box, a_dir, m_level);
  
  for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
  {
    int ils = m_eqSys->lsIndexCons(i);
    int ivf = m_PhPr->lsIndexField(i);  
    
    FORT_UPDATELEVELSET(CHF_FRA(a_U),  
                    CHF_CONST_FRA(a_W),                  
                    CHF_CONST_FRA(a_F),
                    CHF_CONST_FRA1(dx,0),
                    CHF_CONST_FRA1(g,0),
                    CHF_CONST_REAL(a_dt),
                    CHF_CONST_INT(a_dir),                      
                    CHF_CONST_INT(ils),
                    CHF_CONST_INT(ivf),
                    CHF_BOX(a_box));      
  }

}

              // Update the conserved variable using fluxes and a scaling factor
void PatchMHDAM::finalUpdate(      FArrayBox&       a_U,                      
                                   const FArrayBox& a_W,
                                   const FArrayBox& a_F,                                   
                                         Real       a_dt,
                                   const FArrayBox& a_scale,
                                   const int&       a_dir,
                                   const Box&       a_box)
{
  FArrayBox Fnew(a_F.box(),a_F.nComp());
  Fnew.copy(a_F);
      
  m_csh->transFluxesForUpdateState(Fnew, Fnew.box(), a_dir, m_level);
    
  
  if (m_bLSonly==false)
  updateCons    (a_U,      Fnew, a_dt, a_scale, a_dir, a_box);      
     
  updateLevelSet(a_U, a_W, Fnew, a_dt,          a_dir, a_box);       
}

void PatchMHDAM::EulerStep(         FArrayBox& a_U,
                                    FArrayBox  a_F[CH_SPACEDIM],
                                    EdgeBox  & a_E,
                                    FArrayBox& a_divB,
                              const FArrayBox& a_scale,
                              const FArrayBox  a_areas[CH_SPACEDIM],
                              const FArrayBox& a_S,
                              const Real&      a_dt,
                              const Box&       a_box)
{    
  // Get the number of various variables
  int numFlux  = m_eqSys->numFluxes();
  int numPrim  = m_eqSys->numPrimitives();
  int numSlope = m_eqSys->numSlopes();
  int numCons  = a_U.nComp();
  
  int dir;
  
  Box UBox     = a_U.box();
  Box WBox     = UBox;
  Box WBoxB    = UBox;
      WBoxB   &= m_domain;
  
  FArrayBox Uold( UBox, numCons );
  Uold.copy( a_U );
                                     // Primitive variables    
  FArrayBox W( WBox, numPrim );
                                     // Source terms arrays
  FArrayBox S( UBox, numCons );
      
  Box valueBox = a_box;  
  valueBox &= m_domain;
  
  Box slopeBox = valueBox;
  slopeBox.grow( 1 );
  
  // Calculate the primitive variables from the conserved variables
  m_eqSys->stateToPrim( W, a_U, WBoxB );
  m_csh->transCartesianVectToCurv(W, WBoxB, m_level);

  ///////#pragma omp for private(dir) num_threads(CH_SPACEDIM)
  for( dir = 0; dir < SpaceDim; dir++ )
  {
    m_PhPr->fillGhostCells( W, a_U, dir, m_currentTime );
  }
  
  BaseFab<int> Region;

  preprocessing( W, S, Region, slopeBox );
    
  FArrayBox WMinus[SpaceDim];  
  FArrayBox WPlus [SpaceDim];
  FluxBox   Bn(valueBox,1);  //FArrayBox     Bn[SpaceDim];  
    

////////////////////////////////////////////////////////////////  Reconstruction

  //#pragma omp parallel for private(dir) num_threads(CH_SPACEDIM) 
  for( dir = 0; dir < SpaceDim; dir++ )
  {
    WMinus[dir].clear();
    WPlus [dir].clear();

                    // Size the intermediate, extrapolated primitive variables
    WMinus[dir].resize( slopeBox, numPrim  );
    WPlus [dir].resize( slopeBox, numPrim  );

                                       // Compute primitive variables on faces
    m_Reconstruction->faceValues( W, WMinus[dir], WPlus[dir], m_level, dir, slopeBox, m_eqSys, m_csh );
    
              
    m_Reconstruction->checkPositivity( WMinus[dir], WPlus[dir], W, slopeBox );
  }


//////////////////////////////////////////////////////  Final fluxes calculation

  // Boxes for face centered state - used for MUSCL approach
  Box faceBox[SpaceDim];  
  // Boxes for face centered fluxes
  Box fluxBox[SpaceDim];

  for( dir = 0; dir < SpaceDim; ++dir )
  {
     //  This box is face centered in direction "dir", is one bigger than the
     // input box in all directions except "dir" and stays one cell away from
     // the domain boundary in "dir"
    faceBox[dir] = valueBox;
    faceBox[dir].grow(dir,1);
    faceBox[dir] &= m_domain;
    faceBox[dir].grow(dir,-1);
    faceBox[dir].surroundingNodes(dir);

               // This box is face centered in direction "dir", is one bigger
               // than the input box in all directions except "dir"
    fluxBox[dir] = valueBox;    
    fluxBox[dir].surroundingNodes(dir);

    FArrayBox& F = a_F[dir];
               F.resize( fluxBox[dir], numFlux );
                   
  }

  if (m_bLSonly == false)
  {
  
    
        //#pragma omp parallel for private(dir) num_threads(CH_SPACEDIM)         
        for( dir = 0; dir < SpaceDim; dir++ )
        {
          FArrayBox& F = a_F[dir];
          
          //#pragma omp task firstprivate(dir) untied
          {
          computeBn( Bn[dir], W, WMinus[dir], WPlus[dir], m_iDivBMethod, dir, faceBox[dir]);
          if (F.nComp()>WBX)  F.copy(Bn[dir],faceBox[dir],0,faceBox[dir],WBX+dir,1);
          }
          
          //#pragma omp task firstprivate(dir) untied
          {
          // Solve the Riemann problem and get fluxes, these final fluxes will be returned
          fluxesRP( F, WPlus[dir], WMinus[dir], dir, faceBox[dir] );    
          }
          
          //#pragma omp task firstprivate(dir) untied
          {
                    // Use the user supplied PhysBC object to obtain boundary fluxes
          m_PhPr->fluxBC( F, Bn[dir], WMinus[dir], WPlus [dir], dir, Side::Lo, m_currentTime );
          m_PhPr->fluxBC( F, Bn[dir], WMinus[dir], WPlus [dir], dir, Side::Hi, m_currentTime );
          }
              
          m_csh->multiplyFluxesByArea(F, m_eqSys->consInterval(), a_areas[dir], F.box(), dir, m_level);
        }
        
        //#pragma omp taskwait
      
      
      /*#pragma omp for num_threads(CH_SPACEDIM) private(dir)
      for( dir = 0; dir < SpaceDim; dir++ )
      {
        FArrayBox& F = a_F[dir];
        m_csh->multiplyFluxesByArea(F, m_eqSys->consInterval(), a_areas[dir], F.box(), dir, m_level);
      }*/
        
    
    
  }
  
  
  for( dir = 0; dir < SpaceDim; dir++ )
  {
    FArrayBox& F = a_F[dir];
    levelSetFluxes( F, WPlus[dir], WMinus[dir], dir, fluxBox[dir] );
  }
  
  
                           // Update conserved variables to the next time step
                                     // using the final flux in each direction
  //#pragma omp parallel for private(dir) num_threads(CH_SPACEDIM)                                    
  for( dir = 0; dir < SpaceDim; dir++ )
  {                  
    finalUpdate( a_U, W, a_F[dir], a_dt, a_scale, dir, valueBox );    
  }

/////////////////////////////////////////////////////// Source terms calculation

// Problem specific postprocessing
  //m_PhPr->postprocessing( a_U, W, a_dt, m_currentTime, valueBox );


//FFF take out for ugradp
  //with this section PUI pressure is written with pdivu source term
  FluxBox Un;
  Un.define(valueBox,1);
  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    computeUn( Un[dir], W, WMinus[dir], WPlus[dir], m_iDivBMethod, dir, fluxBox[dir]);
  }
//FFF end take out for ugradp



  if (m_bLSonly == false)
  addExplicitSources( a_U, W, a_S, S, Bn, a_divB, Un, Region, a_dt, a_scale, valueBox );
//FFF old below for Ugradp
//  addExplicitSources( a_U, W, a_S, S, Bn, a_divB, Region, a_dt, a_scale, valueBox );


///////////////////////////////////////////////////// Electric field calculation
 

  EdgeBox E;
  FArrayBox WOld( WBox, numPrim );  
  Box EBox = a_E.box();
  Interval EInt( 0, 0 );

  if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
  {            
    E.resize( valueBox );

    calculateElectricField( W, WOld, a_F, E, a_dt, valueBox );

    recalculateMagneticField( a_U, Uold, E, a_dt, valueBox );

    a_E.copy( EBox, EInt, E, EInt );
  }

////////////////////////////////////////////////////////////////  Postprocessing

                                            // Problem specific postprocessing
  m_PhPr->postprocessing( a_U, W, a_dt, m_currentTime, valueBox );

  postprocessing( a_U, Uold, a_dt, valueBox );

}                              

void PatchMHDAM::updateHancock(            FArrayBox & a_U,
                                     const FArrayBox & a_W,
                                     const FArrayBox & a_FMinus,
                                     const FArrayBox & a_FPlus,
                                     Real              a_dt,
                                     const FArrayBox & a_scale, 
                                     const FArrayBox   a_areas[CH_SPACEDIM],                                    
                                     const int       & a_dir,
                                     const Box       & a_box)
{
  if (m_eqSys->numTrackingSurfaces()>0)
  {
    FArrayBox dx,g(a_box,1);
  
    IntVect iv_off(IntVect::Zero);
    iv_off[a_dir]=1;

    Box dxBox(a_box.smallEnd()*iv_off, 
              a_box.bigEnd()  *iv_off);            

    dx.define(dxBox,1);                     
    m_csh->dx(dx,dxBox,a_dir,m_level);          
        
    m_csh->sqrtMetricTensorCoeff(g, a_box, a_dir, m_level);
  
    for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
    {
      int ils = m_eqSys->lsIndexCons(i);
      int ivf = m_PhPr->lsIndexField(i);  
      
      FORT_UPDATELEVELSETHANCOCK(CHF_FRA(a_U),
                    CHF_CONST_FRA(a_W),
                    CHF_CONST_FRA(a_FMinus),
                    CHF_CONST_FRA(a_FPlus),
                    CHF_CONST_FRA1(dx,0),
                    CHF_CONST_FRA1(g,0),
                    CHF_CONST_REAL(a_dt),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(ils),
                    CHF_CONST_INT(ivf),
                    CHF_BOX(a_box));     

    }
  }
  
  if (m_bLSonly == false)
  {
    int iBgn = m_eqSys->consInterval().begin();
    int iEnd = m_eqSys->consInterval().end();
            
           // Cast away "const" inputs so their boxes can be shifted left or right
           // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftFPlus  = (FArrayBox&)a_FPlus;
    FArrayBox& shiftFMinus = (FArrayBox&)a_FMinus;

                              // Shift the Plus and Minus primitive variable boxes
                              // 1/2 cell so they are face centered
    shiftFPlus .shiftHalf(a_dir, 1);
    shiftFMinus.shiftHalf(a_dir,-1);  
    
    //m_csh->multiplyFluxesByArea(shiftFPlus,  m_eqSys->consInterval(), a_areas[a_dir], shiftFPlus.box(),  a_dir, m_level);
    //m_csh->multiplyFluxesByArea(shiftFMinus, m_eqSys->consInterval(), a_areas[a_dir], shiftFMinus.box(), a_dir, m_level);
    
    m_csh->transFluxesForUpdateState(shiftFPlus,  shiftFPlus.box(),  a_dir, m_level);
    m_csh->transFluxesForUpdateState(shiftFMinus, shiftFMinus.box(), a_dir, m_level);
                              // Shift the Plus and Minus primitive variable boxes
                              // back to their original position
    shiftFPlus .shiftHalf(a_dir,-1);
    shiftFMinus.shiftHalf(a_dir, 1);
 

                                      
  
    FORT_UPDATECONSHANCOCK(CHF_FRA(a_U),
                    CHF_CONST_FRA(a_FMinus),
                    CHF_CONST_FRA(a_FPlus),
                    CHF_CONST_INT(iBgn),
                    CHF_CONST_INT(iEnd),
                    CHF_CONST_REAL(a_dt),
                    CHF_CONST_FRA1(a_scale,0),
                    CHF_CONST_FRA1(a_areas[a_dir],0),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(a_box));  
  }

//FF an alternative to UPDATELEVELSETHANCOCK used above
//    for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
//    {
//      int ils = m_eqSys->lsIndexCons(i);
//
//
//       FORT_UPDATELEVELSETHANCOCK2(CHF_FRA(a_U),
//                    CHF_CONST_FRA(a_W),
//                    CHF_CONST_FRA(a_FMinus),
//                    CHF_CONST_FRA(a_FPlus),
//                    CHF_CONST_INT(ils),
//                    CHF_CONST_REAL(a_dt),
//                    CHF_CONST_FRA1(a_scale,0),
//                    CHF_CONST_FRA1(a_areas[a_dir],0),
//                    CHF_CONST_INT(a_dir),
//                    CHF_BOX(a_box));
//
//    }



  
                  
                      
}

void PatchMHDAM::updateValuesOnFaces(       FArrayBox & a_WMinus,
                                                 FArrayBox & a_WPlus,
                                           const FArrayBox & a_W,
                                           const FArrayBox & a_WOld,
                                           const Box&        a_box )
{
  FORT_UPDATEVALUESONFACES( CHF_FRA(a_WMinus),
                            CHF_FRA(a_WPlus),
                            CHF_CONST_FRA(a_W),
                            CHF_CONST_FRA(a_WOld),
                            CHF_BOX(a_box) );
}


void PatchMHDAM::predictorHancock(  FArrayBox& a_U,    
                                    FArrayBox& a_WOld,
                                    FArrayBox& a_WPredictor,
                                    FArrayBox  a_WMinus[CH_SPACEDIM],
                                    FArrayBox  a_WPlus[CH_SPACEDIM],
                                    EdgeBox  & a_E,  
                                    FluxBox&   a_Bn,
                                    FArrayBox& a_divB,
                              const FArrayBox& a_scale,
                              const FArrayBox  a_areas[CH_SPACEDIM],
                              const FArrayBox& a_SOut,                              
                              const Real&      a_dt,
                              const Box&       a_box)
{    
  int numFlux  = m_eqSys->numFluxes();
  int numPrim  = m_eqSys->numPrimitives();
  int numSlope = m_eqSys->numSlopes();
  int numCons  = a_U.nComp();
  
  Box UBox     = a_U.box();
  Box WBox     = UBox;
  Box WBoxB    = UBox;
      WBoxB   &= m_domain;  
      
                                     // Primitive variables    
  FArrayBox W( WBox, numPrim );
                                     // Source terms arrays
  FArrayBox S( UBox, numCons );
      
  Box valueBox = a_box;
  valueBox.grow( 1 );    
  valueBox &= m_domain;
  
  Box slopeBox = valueBox;
  slopeBox.grow( 1 );
  
  // Calculate the primitive variables from the conserved variables
  m_eqSys->stateToPrim( W, a_U, WBoxB );
  m_csh->transCartesianVectToCurv(W, WBoxB, m_level);

  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    m_PhPr->fillGhostCells( W, a_U, dir, m_currentTime );
  }
  
  BaseFab<int> Region;

  preprocessing( W, S, Region, slopeBox );
  
  FArrayBox FPlus [SpaceDim];    
  FArrayBox FMinus[SpaceDim];    
  
  //FArrayBox Bn[SpaceDim];
  if (a_Bn.box().isEmpty()) a_Bn.define(valueBox,1);
    
////////////////////////////////////////////////////////////////  Reconstruction

  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    a_WMinus[dir].clear();
    a_WPlus [dir].clear();

                    // Size the intermediate, extrapolated primitive variables
    a_WMinus[dir].resize( slopeBox, numPrim  );
    a_WPlus [dir].resize( slopeBox, numPrim  );

                                       // Compute primitive variables on faces
    m_Reconstruction->faceValues( W, a_WMinus[dir], a_WPlus[dir], m_level, dir, slopeBox, m_eqSys, m_csh );

    m_Reconstruction->checkPositivity( a_WMinus[dir], a_WPlus[dir], W, slopeBox );
  }

//////////////////////////////////////////////////////  Final fluxes calculation
  
  // Boxes for face centered fluxes
  Box fluxBox[SpaceDim];

  for( int dir = 0; dir < SpaceDim; ++dir )
  {
               // This box is face centered in direction "dir", is one bigger
               // than the input box in all directions except "dir"
    fluxBox[dir] = valueBox;
    //fluxBox[dir].grow(1);
    //fluxBox[dir].grow(dir,-1);
    //fluxBox[dir] &= m_domain;
    fluxBox[dir].surroundingNodes(dir);
    
    FPlus [dir].resize( valueBox, numCons  );    
    FMinus[dir].resize( valueBox, numCons  );      

    //Bn[dir].clear();
    //Bn[dir].resize( fluxBox[dir], 1 );
  }
  
  //FArrayBox scale(valueBox,1);  
  //m_csh->scalingFactor(scale, valueBox, m_level);
  //scale *= a_dt;

  if (m_bLSonly == false)
  for( int dir = 0; dir < SpaceDim; dir++ )
  {        
    fluxesHancock( FMinus[dir], FPlus[dir], a_WMinus[dir], a_WPlus[dir],  dir, valueBox );
            
    computeBn( a_Bn[dir], W, a_WMinus[dir], a_WPlus[dir], m_iDivBMethod, dir, fluxBox[dir]);              
  }  
  
  for( int dir = 0; dir < SpaceDim; dir++ )
  {                
    updateHancock( a_U, W, FMinus[dir], FPlus[dir], a_dt, a_scale, a_areas, dir, valueBox );        
  }  
  
//take out for ugradp
//with this section PUI pressure is written with pdivu source term
  FluxBox Un;
  Un.define(valueBox,1);
  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    computeUn( Un[dir], W, a_WMinus[dir], a_WPlus[dir], m_iDivBMethod, dir, fluxBox[dir]);
  }
//end take out for ugradp

/////////////////////////////////////////////////////// Source terms calculation

  if (m_bLSonly == false)
  addExplicitSources( a_U, W, a_SOut, S, a_Bn, a_divB, Un, Region, a_dt, a_scale, valueBox );
//below is old, does not include ugradp
//  addExplicitSources( a_U, W, a_SOut, S, a_Bn, a_divB, Region, a_dt, a_scale, valueBox );

  a_WOld.copy(W);
  
  m_eqSys->stateToPrim( a_WPredictor, a_U, WBoxB );
  m_csh->transCartesianVectToCurv(a_WPredictor, WBoxB, m_level);
    

  for( int dir = 0; dir < SpaceDim; dir++ )
  {    
    m_PhPr->fillGhostCells( a_WPredictor, a_U, dir, m_currentTime + a_dt);
  }
  
  for (int dir = 0; dir < SpaceDim; dir++)
  {
    updateValuesOnFaces( a_WMinus[dir], a_WPlus[dir], a_WPredictor, a_WOld, slopeBox );
  }
  
  if (m_bLSonly == false)
  for( int dir = 0; dir < SpaceDim; dir++ )
  {               
    computeBn( a_Bn[dir], W, a_WMinus[dir], a_WPlus[dir], m_iDivBMethod, dir, fluxBox[dir]);              
  }
  
  if (m_DivergenceCleaning == DC_ProjHancock)
  {
    computeDivB(a_divB, a_Bn, W, a_box );
    a_divB*=(1.0/m_csh->dx(0,m_level));  
    Real maxdivB = MAX(fabs(a_divB.min()), fabs(a_divB.max()));
    Real sumdivB = a_divB.sumPow(a_divB.box(),2);
    //pout() << "max divB = " << maxdivB << " sum abs(divB) = " << sumdivB << endl;
  }



  //Shock boundary conditions for pickup ions
    PickupIons * pPickUp = m_eqSys->getPickupIons();
      if( pPickUp != NULL )
      {
        if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
        {
//          m_PhPr->shockBC(a_WPredictor, a_U, Region, m_level);
          FArrayBox Wtest( UBox, numCons );
          m_eqSys->stateToPrim( Wtest, a_U, WBoxB );
          m_csh->transCartesianVectToCurv(Wtest, WBoxB, m_level);

          BaseFab<int> region;

          preprocessing( Wtest, S, region, slopeBox );
          m_PhPr->shockBC(Wtest, a_U, region, m_level);
        }
      }
}

void PatchMHDAM::correctorHancock(  FArrayBox& a_U,
                              const FArrayBox& a_W,                   
                              const FArrayBox& a_WPredictor,                   
                                    FArrayBox  a_F[CH_SPACEDIM],
                                    FArrayBox  a_WMinus[CH_SPACEDIM],
                                    FArrayBox  a_WPlus[CH_SPACEDIM],
                                    EdgeBox  & a_E, 
                                    FluxBox  & a_Bn,
                                    FArrayBox& a_divB,
                              const FArrayBox& a_scale,
                              const FArrayBox  a_areas[CH_SPACEDIM],
                              const FArrayBox& a_SOut,                              
                              const Real&      a_dt,
                              const Box&       a_box)
{
  // Get the number of various variables
  int numFlux  = m_eqSys->numFluxes();
  int numPrim  = m_eqSys->numPrimitives();
  int numSlope = m_eqSys->numSlopes();
  int numCons  = a_U.nComp();
  
  Box UBox     = a_U.box();
  Box WBox     = UBox;
  Box WBoxB    = UBox;
      WBoxB   &= m_domain;
  
  FArrayBox Uold( UBox, numCons );
  Uold.copy( a_U );                                     
                                     // Source terms arrays
  FArrayBox S( UBox, numCons );
      
  Box valueBox = a_box;  
  valueBox &= m_domain;
      
  BaseFab<int> Region;

  preprocessing( a_WPredictor, S, Region, valueBox );
      
  //FArrayBox     Bn[SpaceDim];
    
//////////////////////////////////////////////////////  Final fluxes calculation

  // Boxes for face centered state - used for MUSCL approach
  Box faceBox[SpaceDim];  
  // Boxes for face centered fluxes
  Box fluxBox[SpaceDim];

  for( int dir = 0; dir < SpaceDim; ++dir )
  {
     //  This box is face centered in direction "dir", is one bigger than the
     // input box in all directions except "dir" and stays one cell away from
     // the domain boundary in "dir"
    faceBox[dir] = valueBox;
    faceBox[dir].grow(dir,1);//faceBox[dir].grow(1);
    faceBox[dir] &= m_domain;
    faceBox[dir].grow(dir,-1);
    faceBox[dir].surroundingNodes(dir);

               // This box is face centered in direction "dir", is one bigger
               // than the input box in all directions except "dir"
    fluxBox[dir] = valueBox;
    //fluxBox[dir].grow(1);
    //fluxBox[dir].grow(dir,-1);
    //fluxBox[dir] &= m_domain;
    fluxBox[dir].surroundingNodes(dir);

    FArrayBox& F = a_F[dir];
               F.resize( fluxBox[dir], numFlux );

    //Bn[dir].clear();
    //Bn[dir].resize( fluxBox[dir], 1 );    
    //Bn[dir].copy(a_Bn[dir]);
  }
  
  FArrayBox divB(a_box,1);divB.setVal(0.0);
  
  if (m_bLSonly == false)
  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    FArrayBox& F = a_F[dir];
    
    //computeBn( Bn[dir], a_WPredictor, a_WMinus[dir], a_WPlus[dir], m_iDivBMethod, dir, faceBox[dir]);
    //correctBn( Bn[dir], a_phi, dir, faceBox[dir]);
    if (F.nComp()>WBX)  F.copy(a_Bn[dir],fluxBox[dir],0,fluxBox[dir],WBX+dir,1);
        
    // Solve the Riemann problem and get fluxes, these final fluxes will be returned
    fluxesRP( F, a_WPlus[dir], a_WMinus[dir], dir, faceBox[dir] );    
    
              // Use the user supplied PhysBC object to obtain boundary fluxes
    m_PhPr->fluxBC( F, a_Bn[dir], a_WMinus[dir], a_WPlus [dir], dir, Side::Lo, m_currentTime );
    m_PhPr->fluxBC( F, a_Bn[dir], a_WMinus[dir], a_WPlus [dir], dir, Side::Hi, m_currentTime );
        
    m_csh->multiplyFluxesByArea(F, m_eqSys->consInterval(), a_areas[dir], F.box(), dir, m_level);
  }
  
  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    FArrayBox& F = a_F[dir];
    levelSetFluxes( F, a_WPlus[dir], a_WMinus[dir], dir, fluxBox[dir] );
  }
    

  //FArrayBox scale(valueBox,1);  
  //m_csh->scalingFactor(scale, valueBox, m_level);
  //scale *= a_dt;
                           // Update conserved variables to the next time step
                                     // using the final flux in each direction
  for( int dir = 0; dir < SpaceDim; dir++ )
  {                  
    finalUpdate( a_U, a_W, a_F[dir], a_dt, a_scale, dir, valueBox );
  }
/////////////////////////////////////////////////////// Source terms calculation

//take out for ugradp
//same as in predictor for PUI source term
  FluxBox Un;
  Un.define(valueBox,1);
  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    computeUn( Un[dir], a_WPredictor, a_WMinus[dir], a_WPlus[dir], m_iDivBMethod, dir, fluxBox[dir]);
  }
//end take out for ugradp
 
  if (m_bLSonly == false)
  addExplicitSources( a_U, a_WPredictor, a_SOut, S, a_Bn, a_divB, Un, Region, a_dt, a_scale, valueBox );
//old below, same as in predictor
//  addExplicitSources( a_U, a_WPredictor, a_SOut, S, a_Bn, a_divB, Region, a_dt, a_scale, valueBox );

///////////////////////////////////////////////////// Electric field calculation

  EdgeBox E;
  FArrayBox WOld( WBox, numPrim );  
  Box EBox = a_E.box();
  Interval EInt( 0, 0 );

  if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
  {            
    E.resize( valueBox );

    calculateElectricField( a_WPredictor, WOld, a_F, E, a_dt, valueBox );

    recalculateMagneticField( a_U, Uold, E, a_dt, valueBox );

    a_E.copy( EBox, EInt, E, EInt );
  }

////////////////////////////////////////////////////////////////  Postprocessing

                                            // Problem specific postprocessing
  m_PhPr->postprocessing( a_U, a_WPredictor, a_dt, m_currentTime, valueBox );

  postprocessing( a_U, Uold, a_dt, valueBox );
  //Shock boundary conditions for pickup ions
    PickupIons * pPickUp = m_eqSys->getPickupIons();
    if( pPickUp != NULL )
    {
      if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
      {
//        m_PhPr->shockBC(a_WPredictor, a_U, Region, m_level);

        FArrayBox Wtest( UBox, numCons );
        m_eqSys->stateToPrim( Wtest, a_U, WBoxB );
        m_csh->transCartesianVectToCurv(Wtest, WBoxB, m_level);

        BaseFab<int> region;

        Box b = a_box;
        b.grow(1);

        preprocessing( Wtest, S, region, b );
        m_PhPr->shockBC(Wtest, a_U, region, m_level);
      }
    }
}


//        Update the conserved variables, return the final fluxes used for this,
//                          and return the maximum wave speed on this patch/grid
void PatchMHDAM::updateState(       FArrayBox& a_U,
                                    FArrayBox  a_F[CH_SPACEDIM],
                                    EdgeBox  & a_E,
                                    Real&      a_newDt,
                              IntVect&         a_minDtCell,
                              const FArrayBox& a_S,
                              FArrayBox&       a_divB,
                              const Real&      a_dt,
                              const Box&       a_box             )
{
  if (m_verbosity >= 3) {
    pout() << "PatchMHDAM::updateState ";
    pout() << "(" << D_TERM( a_box.smallEnd()[0], << "," << a_box.smallEnd()[1], << "," << a_box.smallEnd()[2]) << ") ";
    pout() << "(" << D_TERM( a_box.bigEnd()[0],   << "," << a_box.bigEnd()[1],   << "," << a_box.bigEnd()[2])   << ")";
    pout() << endl;
  }  
  CH_assert(isDefined());  

#ifndef NDEBUG
  int LCorner[CH_SPACEDIM]={D_DECL(a_box.smallEnd()[0], a_box.smallEnd()[1], a_box.smallEnd()[2])};
  int UCorner[CH_SPACEDIM]={D_DECL(a_box.bigEnd()[0],   a_box.bigEnd()[1],   a_box.bigEnd()[2])};
#endif
  
  int dir;

  int numStates  = a_U.nComp();
  Box UBox = a_U.box();
  
  // Calculate geometrical info
  FArrayBox scale;                        // Inversion of cell volumes
  FArrayBox areas[CH_SPACEDIM];           // Cell areas
  Box geom_box = a_box;
  geom_box.grow(2); 
  scale.define(geom_box,1);
  m_csh->scalingFactor(scale, geom_box, m_level);
  Box areasBox;
  for (dir = 0; dir < SpaceDim; ++dir)
  {
    areasBox = geom_box;
    areasBox.surroundingNodes(dir);
    areas[dir].define(areasBox,1);
    m_csh->getAreas(areas[dir], areasBox, dir, m_level);
  }
    
    
  if( m_TimeMethod == TA_RK2 )
  {    
//RK2 step does not contain shock BCs or ugradp for PUIs
    // Create a temp storage for a_U
    FArrayBox Uold( UBox, numStates );
    Uold.copy( a_U );
  
    FArrayBox  FPredictor[CH_SPACEDIM];
    FArrayBox  UPredictor(UBox, numStates);   
    EdgeBox    EPredictor; 
    Box        BPredictor(a_box);    
    
    UPredictor.copy( a_U );
    BPredictor.grow( 2 );    
    if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
      EPredictor.resize(a_E.box());
    
    
    EulerStep( UPredictor, FPredictor, EPredictor, a_divB, scale, areas, a_S,     a_dt, BPredictor);
    a_U.copy( UPredictor );

    EulerStep( a_U,        a_F,        a_E,        a_divB, scale, areas, a_S, 0.5*a_dt, a_box);
    
      
    //#pragma omp parallel for private(dir) //num_threads(CH_SPACEDIM) 
    for( dir = 0; dir < SpaceDim; dir++ )
    {
      a_F[dir]+=FPredictor[dir];
      a_F[dir]*=0.5;
    }
    
    
    if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
    {
      a_E+=EPredictor;
      a_E*=0.5;
    }
    
    int iBGN = 0;
    int iEND = numStates - 1;
    
    FORT_AVERAGESTATE( CHF_FRA(a_U),
                       CHF_CONST_FRA(Uold),
                       CHF_CONST_FRA(UPredictor),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(a_box) );

    
  } else if (m_TimeMethod == TA_Hancock)
  {
    FArrayBox  dummyPhi(UBox, 1);dummyPhi.setVal(0.0);
    FArrayBox  W(UBox, numStates);  
    FArrayBox  UPredictor(UBox, numStates);       
    FArrayBox  WPredictor(UBox, numStates);  
    FluxBox    Bn;  
    UPredictor.copy( a_U );
    
    FArrayBox WMinus[CH_SPACEDIM], WPlus[CH_SPACEDIM];
    
    //scale*=0.5*a_dt;
    predictorHancock(UPredictor, W, WPredictor,      WMinus, WPlus, a_E, Bn, a_divB, scale, areas, a_S, 0.5*a_dt, a_box);
    //scale*=2.0;
    correctorHancock(a_U,        W, WPredictor, a_F, WMinus, WPlus, a_E, Bn, a_divB, scale, areas, a_S,     a_dt, a_box);

  } else if (m_TimeMethod == TA_FirstOrder)
  {
//Euler step does not contain shock BCs or ugradp for PUIs
    //scale*=a_dt;
    EulerStep( a_U, a_F, a_E, a_divB, scale, areas, a_S, a_dt, a_box);
  }

                     // Get next time step for this patch                           
  a_newDt = computeDt( a_U, a_box, a_minDtCell);

  if (m_verbosity >= 3) {
    pout() << "Leave PatchMHDAM::updateState "  << endl;
  }
}

//                                Return true if everything is defined and setup
bool PatchMHDAM::isDefined() const
{
  return m_isDefined        &&
         m_isPPSet          &&
         m_isRSSet          &&
         m_isRecSet         &&
         m_isCurrentTimeSet ;
}

void PatchMHDAM::postTimeStep( void )
{
}

/*
////////////////////////////////////////////////////

//        Update the conserved variables, return the final fluxes used for this,
//                          and return the maximum wave speed on this patch/grid
void PatchMHDAM::updateState(       FArrayBox& a_U,
                                    FArrayBox  a_F[CH_SPACEDIM],
                                    EdgeBox  & a_E,
                                    Real&      a_newDt,
                              IntVect&         a_minDtCell,
                              const FArrayBox& a_S,
                              FArrayBox&       a_divB,
                              const Real&      a_dt,
                              const Box&       a_box             )
{
  if (m_verbosity >= 3) {
    pout() << "PatchMHDAM::updateState ";
    pout() << "(" << D_TERM( a_box.smallEnd()[0], << "," << a_box.smallEnd()[1], << "," << a_box.smallEnd()[2]) << ") ";
    pout() << "(" << D_TERM( a_box.bigEnd()[0],   << "," << a_box.bigEnd()[1],   << "," << a_box.bigEnd()[2])   << ")";
    pout() << endl;
  }

  CH_assert(isDefined());
  CH_assert(a_box == m_currentBox);

#ifndef NDEBUG
  int LCorner[CH_SPACEDIM]={D_DECL(a_box.smallEnd()[0], a_box.smallEnd()[1], a_box.smallEnd()[2])};
  int UCorner[CH_SPACEDIM]={D_DECL(a_box.bigEnd()[0],   a_box.bigEnd()[1],   a_box.bigEnd()[2])};
#endif

                                          // Get the number of various variables
  int numFlux  = m_eqSys->numFluxes();
  int numPrim  = m_eqSys->numPrimitives();
  int numSlope = m_eqSys->numSlopes();
  int numCons  = a_U.nComp();

//////////////////////////////////////////////////////////  All boxes definition

                                                              // The current box
  Box curBox   = m_currentBox;

  Box UBox     = a_U.box();
  Box WBox     = UBox;
  Box WBoxB    = UBox;
      WBoxB   &= m_domain;
      

                      // Boxes for face centered state - used for MUSCL approach
  Box faceBox[SpaceDim];
                                               // Boxes for face centered fluxes
  Box fluxBox[SpaceDim];

  Box valueBox = curBox;
  Box BnBox    = curBox;

                                   // Define the box where slopes will be needed
                                   // (one larger than the final update box)
  Box slopeBox = curBox;
  slopeBox.grow(1);
  slopeBox  &= m_domain;

                                                     // Boxes for electric field
  Box EBox = a_E.box();
  Interval EInt( 0, 0 );

////////////////////////////////////////////////////////  All array declarations

                                                // Create a temp storage for a_U
  FArrayBox Uold( UBox, numCons );
                                                          // Primitive variables
  FArrayBox W   ( WBox, numPrim );
  FArrayBox WOld( WBox, numPrim );
                                             // Intermediate primitive variables                                             
  FArrayBox WMinus[SpaceDim];
  FArrayBox WPlus [SpaceDim];
                                                          // Source terms arrays
  FArrayBox S( UBox, numCons );

  BaseFab<int> Region;

  EdgeBox E;

////////////////////////////////////////////////  Initial filling of some arrays

  Uold.copy( a_U );

               // Calculate the primitive variables from the conserved variables
  m_eqSys->stateToPrim( W, a_U, WBoxB );

  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    m_PhPr->fillGhostCells( W, dir1, m_currentTime );
  }

  Real finalStep  = a_dt;

////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////// Runge-Kutta first stage
  if( timeApproximation() == TA_RK2 )
  {
    finalStep *= 0.5;

    valueBox.grow( 2 );
    valueBox &= m_domain;
    slopeBox = valueBox;
    slopeBox.grow( 1 );

    preprocessing( W, S, Region, slopeBox );

////////////////////////////////////////////////////////////////  Reconstruction

    for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
    {
      WMinus[dir1].clear();
      WPlus [dir1].clear();

                      // Size the intermediate, extrapolated primitive variables
      WMinus[dir1].resize( slopeBox, numPrim  );
      WPlus [dir1].resize( slopeBox, numPrim  );

                                         // Compute primitive variables on faces
      m_Reconstruction->faceValues( W, WMinus[dir1], WPlus[dir1], m_level, dir1, slopeBox, m_eqSys, m_csh );

      m_Reconstruction->checkPositivity( WMinus[dir1], WPlus[dir1], W, slopeBox );
    }

//////////////////////////////////////////////////////  Final fluxes calculation

    for( int dir1 = 0; dir1 < SpaceDim; ++dir1 )
    {
       //  This box is face centered in direction "dir1", is one bigger than the
       // input box in all directions except "dir1" and stays one cell away from
       // the domain boundary in "dir1"
      faceBox[dir1] = valueBox;
      faceBox[dir1].grow(1);
      faceBox[dir1] &= m_domain;
      faceBox[dir1].grow(dir1,-1);
      faceBox[dir1].surroundingNodes(dir1);

                 // This box is face centered in direction "dir1", is one bigger
                 // than the input box in all directions except "dir1"
      fluxBox[dir1] = valueBox;
      fluxBox[dir1].grow(1);
      fluxBox[dir1].grow(dir1,-1);
      fluxBox[dir1] &= m_domain;
      fluxBox[dir1].surroundingNodes(dir1);

      FArrayBox& F = a_F[dir1];
                 F.resize( fluxBox[dir1], numFlux );

      m_Bn[dir1].clear();
      m_Bn[dir1].resize( fluxBox[dir1], 1 );
    }

    for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
    {
      FArrayBox& F = a_F[dir1];

      // Solve the Riemann problem and get fluxes, these final fluxes will be returned
      fluxesRP( F, WPlus[dir1], WMinus[dir1], dir1, faceBox[dir1] );

      computeBn( m_Bn[dir1], W, WMinus[dir1], WPlus[dir1], m_iDivBMethod, dir1, faceBox[dir1]);

                // Use the user supplied PhysBC object to obtain boundary fluxes
      m_PhPr->fluxBC( F, m_Bn[dir1], WMinus[dir1], WPlus [dir1], dir1, Side::Lo, m_currentTime );
      m_PhPr->fluxBC( F, m_Bn[dir1], WMinus[dir1], WPlus [dir1], dir1, Side::Hi, m_currentTime );

      m_csh->multiplyFluxesByArea(F, m_eqSys->consInterval(), F.box(), dir1, m_level);
    }

    FArrayBox scale(valueBox,1);  
    m_csh->scalingFactor(scale, valueBox, m_level);
    scale *= a_dt;
                             // Update conserved variables to the next time step
                                       // using the final flux in each direction
    for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
    {                  
      finalUpdate( a_U, a_F[dir1], scale, dir1, valueBox );
    }

/////////////////////////////////////////////////////// Source terms calculation

    addExplicitSources( a_U, W, a_S, S, a_divB, Region, a_dt, scale, valueBox );

///////////////////////////////////////////////////// Electric field calculation

    if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
    {            
      E.resize( valueBox );

      calculateElectricField( W, WOld, a_F, E, a_dt, valueBox );

      recalculateMagneticField( a_U, Uold, E, a_dt, valueBox );

      a_E.copy( EBox, EInt, E, EInt );
    }

////////////////////////////////////////////////////////////////  Postprocessing

                                              // Problem specific postprocessing
    m_PhPr->postprocessing( a_U, W, a_dt, m_currentTime, valueBox );

    postprocessing( a_U, Uold, a_dt, valueBox );
  }

////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////// Runge-Kutta second stage

  valueBox = curBox;
  slopeBox = valueBox;
  slopeBox.grow( 1 );

  m_eqSys->stateToPrim( W, a_U, WBoxB );

  // Temporally commented to debug spherical !!!!!!!!!
  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    m_PhPr->fillGhostCells( W, dir1, m_currentTime );
  }

  preprocessing( W, S, Region, slopeBox );

  if( timeApproximation() == TA_RK2 )
  {
    int iBGN = 0;
    int iEND = numCons - 1;

    FORT_AVERAGESTATE( CHF_FRA(a_U),
                       CHF_CONST_FRA(Uold),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(curBox) );
  }

////////////////////////////////////////////////////////////////  Reconstruction

  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    WMinus[dir1].clear();
    WPlus [dir1].clear();

                      // Size the intermediate, extrapolated primitive variables
    WMinus[dir1].resize( slopeBox, numPrim  );
    WPlus [dir1].resize( slopeBox, numPrim  );

                                         // Compute primitive variables on faces
    m_Reconstruction->faceValues( W, WMinus[dir1], WPlus[dir1], m_level, dir1, slopeBox, m_eqSys, m_csh );

    m_Reconstruction->checkPositivity( WMinus[dir1], WPlus[dir1], W, slopeBox );
  }

  if( timeApproximation() == TA_Hancock )
  {
////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////// Hancock's method predictor

    for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      updateHancock( a_U, WMinus[dir1], WPlus[dir1], a_dt, dir1, slopeBox );

      BnBox  = slopeBox;
      BnBox.surroundingNodes(dir1);

      m_Bn[dir1].clear();
      m_Bn[dir1].resize( BnBox, 1 );
      m_Bn[dir1].setVal( 0.0 );

      BnBox  = slopeBox;
      BnBox &= m_domain;
      BnBox.surroundingNodes(dir1);

      computeBn( m_Bn[dir1], W, WMinus[dir1], WPlus[dir1], 2, dir1, BnBox );
    }
    
    FArrayBox scale(slopeBox,1);
                                                     // Source terms calculation
    addExplicitSources( a_U, W, a_S, S, a_divB, Region, 0.5*a_dt, scale, slopeBox );

    WOld.copy( W );

               // Calculate the primitive variables from the conserved variables
    m_eqSys->stateToPrim( W, a_U, WBoxB );

    for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
    {
      m_PhPr->fillGhostCells( W, dir1, m_currentTime );
    }

    for (int dir1 = 0; dir1 < SpaceDim; dir1++)
    {
      updateValuesOnFaces( WMinus[dir1], WPlus[dir1], W, WOld, slopeBox );
    }

    a_U.copy( Uold );
  }

//////////////////////////////////////////////////////  Final fluxes calculation

          // faceBox and fluxBox are now a bit smaller for the final corrections
  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    faceBox[dir1]  = valueBox;
    faceBox[dir1].grow(1);
    faceBox[dir1] &= m_domain;
    faceBox[dir1].grow(dir1,-1);
    faceBox[dir1].surroundingNodes(dir1);

    fluxBox[dir1] = valueBox;
    fluxBox[dir1].grow(1);
    fluxBox[dir1].grow(dir1,-1);
    fluxBox[dir1] &= m_domain;
    fluxBox[dir1].surroundingNodes(dir1);

    FArrayBox& F = a_F[dir1];
               F.resize( fluxBox[dir1], numFlux );
               F.setVal( 0.0 );

    m_Bn[dir1].clear();
    m_Bn[dir1].resize( fluxBox[dir1], 1 );
    m_Bn[dir1].setVal( 0.0 );
  }

  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    FArrayBox& F = a_F[dir1];

    // Solve the Riemann problem and get fluxes, these final fluxes will be returned
    fluxesRP( F, WPlus[dir1], WMinus[dir1], dir1, faceBox[dir1] );

    computeBn( m_Bn[dir1], W, WMinus[dir1], WPlus[dir1], m_iDivBMethod, dir1, faceBox[dir1]);

                // Use the user supplied PhysBC object to obtain boundary fluxes
    m_PhPr->fluxBC( F, m_Bn[dir1], WMinus[dir1], WPlus [dir1], dir1, Side::Lo, m_currentTime );
    m_PhPr->fluxBC( F, m_Bn[dir1], WMinus[dir1], WPlus [dir1], dir1, Side::Hi, m_currentTime );       

    m_csh->multiplyFluxesByArea(F, m_eqSys->consInterval(), F.box(), dir1, m_level); 
  }
  
  FArrayBox scale(curBox,1);  
  m_csh->scalingFactor(scale, curBox, m_level);
  scale *= finalStep;

                             // Update conserved variables to the next time step
                                       // using the final flux in each direction
  for( int dir1 = 0; dir1 < SpaceDim; dir1++ )
  {
    finalUpdate( a_U, a_F[dir1], scale, dir1, curBox );
  }

/////////////////////////////////////////////////////// Source terms calculation

  addExplicitSources( a_U, W, a_S, S, a_divB, Region, finalStep, scale,curBox );

///////////////////////////////////////////////////// Electric field calculation

  if ((m_DivergenceCleaning == DC_CT_BS) || (m_DivergenceCleaning == DC_CT_GS0) || (m_DivergenceCleaning == DC_CT_GS1))
  {        
    E.resize( valueBox );

    calculateElectricField( W, WOld, a_F, E, a_dt,  valueBox );

    if( timeApproximation() == TA_RK2 )
    {
      a_E += E;
      a_E *= 0.5;
    }
    else
    {
      a_E.copy( EBox, EInt, E, EInt );
    }
  }

////////////////////////////////////////////////////////////////  Postprocessing

                                              // Problem specific postprocessing
  m_PhPr->postprocessing( a_U, W, finalStep, m_currentTime, curBox );

  postprocessing( a_U, Uold, a_dt, curBox );

                     // Get next time step for this patch
  a_newDt = computeDt( a_U, curBox, a_minDtCell);

  if (m_verbosity >= 3) {
    pout() << "Leave PatchMHDAM::updateState "  << endl;
  }
}*/

