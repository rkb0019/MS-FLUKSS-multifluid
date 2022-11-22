#include <iostream>
#include <iomanip>
#include <stdlib.h>

using std::ifstream;
using std::ios;

#include "DebugOut.H"
#include "DebugF_F.H"

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SWLISMProblemPolar.H"
#include "SWLISMF_F.H"
#include "ChargeExchange2F_F.H"

#include "PatchMHDAMF_F.H"
#include "PatchIdealMHDF_F.H"
#include "RiemannSolver.H"
#include "LGintegrator.H"
#include "MHDAMDefs.H"
#include "TecplotIO.H"
#include "EosCommon.H"
#include "EqSysMHDMF.H"


// Null constructor
SWLISMProblemPolar::SWLISMProblemPolar(ePhysicalModel ePM)
  : MultiFluidProblem()
{
  m_physModel = ePM;
  m_isFortranCommonSet = false;
    
  m_stop_ref = -1.0;
  m_max_levelTS = 10000;
  m_adaptRmin   = 0.0;
}

// Null constructor
SWLISMProblemPolar::SWLISMProblemPolar()
  : MultiFluidProblem()
{
  m_physModel = PP_Undefined;
  m_isFortranCommonSet = false;
    
  m_stop_ref = -1.0;
  m_max_levelTS = 10000;
}


// Input parameters
void SWLISMProblemPolar::input( ParmParse & parser, int verbosity )
{  
  // Read problem specific parameters here
  
  parser.query( "gamma",   m_gamma );
  parser.query( "lismN",   m_lismN );
  parser.query( "lismV",   m_lismV );
  parser.query( "lismT",   m_lismT );
  parser.query( "lismB",   m_lismB );
  parser.query( "sunN",    m_sunN  );
  parser.query( "sunV",    m_sunV  );
  parser.query( "sunT",    m_sunT  );
  parser.query( "initR",   m_initR );
  
  Real rho        = m_lismN*eos_mp;
  Real pref       = rho*m_lismV*m_lismV;
  Real p          = 2.0*eos_k*m_lismN*m_lismT;
  Real lismP      = p/pref;
  m_lismM         = 1.0/sqrt( m_gamma*lismP );
  
  m_Bref = 1e+6*sqrt(m_lismN*eos_mp*m_lismV*m_lismV);
  
  parser.query( "adaptRmin", m_adaptRmin );
  
           
  if (parser.contains("fluids"))
  {
    int iFluids  = 1;
    
    parser.query( "fluids",  iFluids   );

    switch( iFluids ){
    case 2  : m_physModel  = PP_2FluidPM; break;
    case 3  : m_physModel  = PP_3FluidPM; break;
    case 4  : m_physModel  = PP_4FluidPM; break;
    case 5  : m_physModel  = PP_5FluidPM; break;
    default : m_physModel  = PP_MHDPM;}
  }

           
    
  parser.query( "stop_ref",  m_stop_ref);
  parser.query( "max_levelTS", m_max_levelTS);
  
  m_TMLIM = 50000.0; // Temperature separating regions 1 and 2
  parser.query( "TMLIM", m_TMLIM);
  
  int i;
  
  // We need to sligtly modify m_sunXC, m_sunYC, m_sunZC in order 
  // to Sun position lies strictly at grid node
  
  Real domainLength = 1.0;
  parser.query("domain_length",domainLength);
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim); 
  for (i = 0; i < SpaceDim; ++i) numCells[i]=0;
  parser.queryarr("num_cells",numCells,0,SpaceDim);

  Real safeDZ = 0.0;
  parser.query("RegSafeDZ",safeDZ);
  m_RegSafeZ   = domainLength - safeDZ;

  m_netN = 0.0;
  if (nFluids()>1)  
  {
    if (!parser.contains("netN")) MayDay::Error("netN must be defined for nultifluid calculations");    
  }
  parser.query( "netN",    m_netN  );  
    
  
  m_verbosity = verbosity;

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "The SW interaction with the LISM problem input:" << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "lismN     = " << m_lismN     << endl;
    pout() << "lismV     = " << m_lismV     << endl;
    pout() << "sunN      = " << m_sunN      << endl;
    pout() << "sunV      = " << m_sunV      << endl;
    pout() << "sunT      = " << m_sunT      << endl;
    pout() << "initR     = " << m_initR     << endl;
    
    if (nFluids()>1)
      pout() << "netN      = " << m_netN  << endl;
  }

  setFortranCommon();
  
      
}


// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void SWLISMProblemPolar::setFortranCommon( )
{
  CH_assert(m_isFortranCommonSet == false);
  
  Real sunXC = 0.0, sunYC = 0.0, sunZC = 0.0;
  Real R0    = 10.0;
  FORT_SETSWLISM( CHF_CONST_REAL( m_gamma ),
                  CHF_CONST_REAL( m_lismN ),
                  CHF_CONST_REAL( m_lismV ),
                  CHF_CONST_REAL( m_lismT ),
                  CHF_CONST_REAL( m_lismB ),
                  CHF_CONST_REAL( sunXC ),
                  CHF_CONST_REAL( sunYC ),
                  CHF_CONST_REAL( sunZC ),
                  CHF_CONST_REAL( R0    ),
                  CHF_CONST_REAL( m_sunN  ),
                  CHF_CONST_REAL( m_sunV  ),
                  CHF_CONST_REAL( m_sunT  ),
                  CHF_CONST_REAL( m_initR ),
                  CHF_CONST_REAL( m_netN  ),
                  CHF_CONST_REAL( m_TMLIM ),
                  CHF_CONST_REAL( m_RegSafeZ));

  if (nFluids()>=1)
  {                    
    Real scaleLen  = 1.5e+13;

    FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( m_lismV  ),
                           CHF_CONST_REAL( scaleLen ),
                           CHF_CONST_REAL( m_lismN  ) );  
  }
                         
  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void SWLISMProblemPolar::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this object.
 */
void SWLISMProblemPolar::define(const ProblemDomain& a_domain,                                  
                                const int            a_level)
{
  MultiFluidProblem::define(a_domain, a_level);
  
  if (nFluids() > 1)
  {
    CH_assert(m_numGhostCells>=2);
    const Box & dBox = m_domain.domainBox();    
    
    Box R0Box( IntVect(D_DECL(dBox.smallEnd()[0],                 dBox.size(1)/2,   dBox.smallEnd()[2])), 
               IntVect(D_DECL(dBox.smallEnd()[0]+m_numGhostCells, dBox.bigEnd()[1], dBox.bigEnd()[2])));
                
    Vector<Box> vstrips;
    vstrips.push_back(R0Box);
      
    Vector<int> procs;
    procs.push_back(0);  
      
    DisjointBoxLayout dbl(vstrips, procs, m_domain);    
    m_R0Strip.define(dbl, 5);    
    
    if (procID() == 0)
    {
      FArrayBox& FAB = m_R0Strip[m_R0Strip.dataIterator()];
      m_R0FAB.define(R0Box,5,FAB.dataPtr());
    } else
    {  
      m_R0FAB.define(R0Box,5);
    }
  }
    
}

void SWLISMProblemPolar::defineMesh(const ProblemDomain & a_prob_domain,
                                    const Vector<Real>  & a_domainBox)
{
  
  return;
  
  IntVect iv_off (IntVect::Zero);
  iv_off[1] = 1;
        
  Box boxD(a_prob_domain.domainBox().smallEnd()*iv_off, 
           a_prob_domain.domainBox().bigEnd()  *iv_off);
           
  FArrayBox dphi(boxD,1);
  
  //dphi.setVal(d_PI/boxD.size(1));
  
  Real dp = d_PI/boxD.size(1);
  
  
  Box b; IntVect iv1,iv2;
  
  iv1 = a_prob_domain.domainBox().smallEnd()*iv_off;
  iv1.shift(1,4);  
  b.define(a_prob_domain.domainBox().smallEnd()*iv_off,iv1);
  dphi.setVal(dp,b,0);
  
  iv2 = a_prob_domain.domainBox().bigEnd()*iv_off;
  iv2.shift(1,-4);  
  b.define(iv2,a_prob_domain.domainBox().bigEnd()*iv_off);
  dphi.setVal(dp,b,0);
    
  Real phi_min = 4.0*dp;  
  Real phi_max = d_PI-4.0*dp;  
  Real expc = 2.0;
  Real phi1,phi2;
  boxD.define(iv1,iv2);
  
  iv1.shift(1, 1);  
  iv2.shift(1,-1);  
  BoxIterator bit(boxD);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
            
    phi1 = phi_min+(phi_max-phi_min)*(exp(expc*(iv[1] - boxD.smallEnd()[1])/boxD.size(1))-1.0)/(exp(expc)-1.0);
    phi2 = phi_min+(phi_max-phi_min)*(exp(expc*(iv[1]+1-boxD.smallEnd()[1])/boxD.size(1))-1.0)/(exp(expc)-1.0);    
            
    //iv[1] = boxD.bigEnd()[1] - iv[1];    
            
    dphi.set(iv, 0, fabs(phi2 - phi1));         
  }
 
  
  /*Real phi_min = 0.0;  
  Real phi_max = d_PI;  
  Real expc = 2.0;
  Real phi1,phi2;
  
  BoxIterator bit(boxD);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
            
    phi1 = phi_min+(phi_max-phi_min)*(exp(expc*(iv[1] - boxD.smallEnd()[1])/boxD.size(1))-1.0)/(exp(expc)-1.0);
    phi2 = phi_min+(phi_max-phi_min)*(exp(expc*(iv[1]+1-boxD.smallEnd()[1])/boxD.size(1))-1.0)/(exp(expc)-1.0);    
            
    //iv[1] = boxD.bigEnd()[1] - iv[1];    
        
    
    dphi.set(iv, 0, fabs(phi2 - phi1));         
  } */   
  
  //dphi.box().p();
  //pout().flush();
  m_csh->setGridSpacing(dphi,1);
  
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* SWLISMProblemPolar::new_PhysProblem()
{
  SWLISMProblemPolar* retval = new SWLISMProblemPolar(m_physModel);
  
  retval->copy_PhysProblem(this);  
  
  return static_cast<PhysProblem*>(retval);
}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void SWLISMProblemPolar::copy_PhysProblem(const PhysProblem* a_PP)
{
  const SWLISMProblemPolar* PP = dynamic_cast<const SWLISMProblemPolar*>(a_PP);
  if (PP == NULL) MayDay::Error("SWLISMProblemPolar::copy_PhysProblem. Wrong argument");
  
  MultiFluidProblem::copy_PhysProblem(a_PP);
  
  if (PP->m_isFortranCommonSet == true)
  {    

    this->m_lismN      = PP->m_lismN;
    this->m_lismM      = PP->m_lismM;
    this->m_lismV      = PP->m_lismV;
    this->m_lismT      = PP->m_lismT;
    this->m_lismB      = PP->m_lismB;
    this->m_Bref       = PP->m_Bref;
    this->m_sunN       = PP->m_sunN;
    this->m_sunV       = PP->m_sunV;
    this->m_sunT       = PP->m_sunT;
    this->m_initR      = PP->m_initR;
    this->m_stop_ref   = PP->m_stop_ref;
    this->m_max_levelTS = PP->m_max_levelTS;

    this->m_gamma      = PP->m_gamma;
    this->m_netN       = PP->m_netN;
    this->m_verbosity  = PP->m_verbosity;
    this->m_TMLIM      = PP->m_TMLIM;
         
    this->setFortranCommonSet();
  }  
  
  this->m_adaptRmin  = PP->m_adaptRmin;
    
}

void SWLISMProblemPolar::preTimeStep(LevelData<FArrayBox>&  a_U, Real a_time)
{
  CH_assert(m_numGhostCells == a_U.ghostVect()[0]);
  
  int vector[1] = {1};
  
  if (nFluids() > 1)        
  {
    m_R0FAB.setVal(-1.0,0);
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);    
    int urho1  = eqSys->densityIndexCons(1);
    int ueng1  = urho1+UNUM_E-1;
    
    a_U.copyTo(Interval(urho1,ueng1), m_R0Strip, m_R0Strip.interval());
  
#ifdef CH_MPI    
    MPI_Bcast(m_R0FAB.dataPtr(),m_R0FAB.box().numPts()*m_R0FAB.nComp(),MPI_CH_REAL,0,Chombo_MPI::comm);
#endif    

    m_csh->transCartesianVectToCurv(m_R0FAB,vector,1,m_R0FAB.box(),m_level);                    
  }

}


                                                             // Fill ghost cells
void SWLISMProblemPolar::fillGhostCells(     FArrayBox&      a_W,
                                       const FArrayBox&      a_U,
                                       const int&            a_dir,
                                       const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int indW,indD;
  
  Box WBox = a_W.box();

                       // See if this chops off the high side of the input box
  Box WBoxDomain  = WBox;
  WBoxDomain     &= m_domain;  
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);
  int iHP    = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  
  int nGS;
  int nGSMax = 4; // How many two ghost cells must be filled
  
  if (a_dir == 0)
  {
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = MIN(abs(indW-indD),nGSMax);
        
    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
          
      FORT_LISMINITPOLAR(
          CHF_FRA(a_W), 
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),            
          CHF_CONST_INT(iHP),            
          CHF_CONST_INT(m_level),
          CHF_BOX(boundaryBox));            
    }

    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
  
                                                         // Fill the ghost cells      
      FORT_SWINITPOLAR(
          CHF_FRA(a_W),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),         
          CHF_CONST_INT(iHP),              
          CHF_CONST_INT(m_level),
          CHF_BOX(boundaryBox));      
            
      if (nFluids() > 1)
      {
        int j_max = m_domain.domainBox().bigEnd()[1];
        FORT_NEUTRALS_SUNBCPOLAR(
          CHF_FRA(a_W),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_FRA(m_R0FAB),
          CHF_CONST_INT(j_max),
          CHF_BOX(boundaryBox));                  
      }
                
    }    
        
  }
  
  if (a_dir == 1)
  {
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    // Axisymmetric case, top boundary
    if( indW != indD )
    {
      int sign         = 1;
                  
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD + 1, nGS ); // We must fill ghost cell in boundaryBox
      
      FORT_POLARGSAXIS(
         CHF_FRA(a_W),       
         CHF_CONST_INT(iRhoN),
         CHF_CONST_INT(fluids), 
         CHF_CONST_INT(iHP),            
         CHF_CONST_INT(sign),
         CHF_CONST_INT(indD),
         CHF_CONST_INT(nGS),
         CHF_BOX(boundaryBox));
                        
    }
  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    // Axisymmetric case, bottom boundary
    if( indW != indD )
    {
      int sign         = -1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS, nGS );
      
      FORT_POLARGSAXIS(
         CHF_FRA(a_W),       
         CHF_CONST_INT(iRhoN),
         CHF_CONST_INT(fluids), 
         CHF_CONST_INT(iHP),            
         CHF_CONST_INT(sign),
         CHF_CONST_INT(indD),
         CHF_CONST_INT(nGS),
         CHF_BOX(boundaryBox));
              
    }
  }  
  
  #ifndef NDEBUG
    FORT_VIEWBOXDATA(
      CHF_FRA(a_W)
      );
  #endif
}


// Set boundary fluxes
void SWLISMProblemPolar::fluxBC(    FArrayBox&      a_F,
                                  FArrayBox&      a_Bn,
                            const FArrayBox&      a_WMinus,
                            const FArrayBox&      a_WPlus,
                            const int&            a_dir,
                            const Side::LoHiSide& a_side,
                            const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);
  
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
      
      #ifndef NDEBUG
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WPlus)
        );
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WMinus)
        );          
      #endif
      
      // Pogorelov's boundary conditions      
                            
      CH_assert(m_RS!=NULL);        
      if (fluids>1) CH_assert(m_RSGD!=NULL);
                  
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
      
      for (int iFluid = 1; iFluid < fluids; ++iFluid)
      { 
        int iRho = eqSys->densityIndexCons(iFluid);        
        m_RSGD->fluxes( a_F, a_WPlus, a_WMinus, a_dir, iRho, boundaryBox );
      }
                                 
                            
      FArrayBox& shiftWLeft  = (FArrayBox&)(a_WPlus);
      FArrayBox& shiftWRight = (FArrayBox&)(a_WMinus);
                          // Shift the left and right primitive variable boxes
                          // 1/2 cell so they are face centered
      shiftWLeft .shiftHalf(a_dir, 1);
      shiftWRight.shiftHalf(a_dir,-1);
      CH_assert( shiftWLeft.box().contains(boundaryBox) );
      CH_assert( shiftWRight.box().contains(boundaryBox) );        
      FORT_SWLISMBCPOLAR( CHF_FRA(a_F),
            CHF_FRA1(a_Bn,0),              
            CHF_CONST_FRA(shiftWLeft), 
            CHF_CONST_FRA(shiftWRight),
            CHF_CONST_INT(sign),
            CHF_CONST_INT(a_dir), 
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),              
            CHF_BOX(boundaryBox) );        
                          // Shift the left and right primitive variable boxes
                          // back to their original position (no net change is made!)
      shiftWLeft .shiftHalf(a_dir,-1);
      shiftWRight.shiftHalf(a_dir, 1);
      
      // Axisymmetric case
      if (a_dir==1) a_Bn.setVal(0.0,boundaryBox,0); 
                                                             
    }
  }
}

// Set up initial conditions
void SWLISMProblemPolar::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);
  int iHP    = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition on this grid    
    FORT_SWLISMINITPOLAR(CHF_CONST_FRA(U),  
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),                        
                    CHF_CONST_INT(iHP),
                    CHF_CONST_INT(m_level),
                    CHF_BOX(uBox));                        
          
  }
  
}


                                              // Problem specific postprocessing
void SWLISMProblemPolar::postprocessing(       FArrayBox & a_U,
                                    const FArrayBox & a_W,
                                    const Real      & a_dt,
                                    const Real      & a_time,
                                    const Box       & a_box       )
{  
                    
}

// Pogorelov's boundary conditions for neutrals
  
void SWLISMProblemPolar::postTimeStep(LevelData<FArrayBox>&  a_U)
{
  
}


void SWLISMProblemPolar::defineRegions( const FArrayBox    & a_W,
                                            FArrayBox    & a_S,
                                            BaseFab<int> & a_R,
                                      const Box          & a_box)
{
  a_R.resize( a_box, 1 );    
                          
  Real dx = -1.0;
  FORT_DEFINE_REGIONS_2F( CHF_CONST_FRA(a_W),
                          CHF_FIA1(a_R,0),       
                          CHF_CONST_REAL(dx),                   
                          CHF_BOX(a_box) );

}


// Number additional variables for writing to plot file
int SWLISMProblemPolar::numPlotVars()
{
  //return 1; // Plasma temperature
  return 2; // Plasma temperature + region
  if (nFluids() > 1) return 2+4;  // Four source terms
  return 0; 
}
  
// Names of the additional variables for writing to plot file  
Vector<string> SWLISMProblemPolar::plotNames()
{
  Vector<string> retval;  
   
  if (numPlotVars() >= 1)
    retval.push_back("T");
  if (numPlotVars() >= 2)  
    retval.push_back("region");
  if (numPlotVars() >= 2+4)  
  {
    retval.push_back("srho");
    retval.push_back("x-smom");
    retval.push_back("y-smom");
    retval.push_back("seng");    
  }    
        
  return retval;    
}
  
// Calculates variables for plotting using primitive variables
void SWLISMProblemPolar::calcPlotVars(FArrayBox& a_vars,
                           int              a_comp,
                           const FArrayBox& a_W,
                           const Box&       a_box)
{
  if (numPlotVars() == 0) return;  
  
  FArrayBox W(a_box, a_W.nComp());  
  W.copy(a_W);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  
  BaseFab<int> Region;
  FArrayBox dummyFab;      
  defineRegions(W, dummyFab, Region, a_box);
    
  Real coeff = m_gamma*m_lismM*m_lismM*m_lismT;
  
  int iT      = a_comp+0;
  int iRegion = a_comp+1;
  
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    Real T = coeff*W(iv,WPRES)/W(iv,WRHO);
    if (numPlotVars() >= 1)  a_vars.set(iv, iT,      T);    
    if (numPlotVars() >= 2)  a_vars.set(iv, iRegion, (Real)(Region(iv, 0)) );    
  }
  
  int isrho   = a_comp+2;
  if (numPlotVars() >= 2+4)  
  {
    int iFluids = nFluids();
    Real dt = 1.0;
    FArrayBox S( a_box, a_W.nComp() );
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    int iRhoN  = eqSys->densityIndexPrim(1);  
    
    if (iFluids == 2)
      FORT_CHARGE_EXCHANGE_2F( CHF_FRA(S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_FIA1(Region,0),
                               CHF_CONST_REAL(dt),  
                               CHF_CONST_INT(iRhoN),                           
                               CHF_BOX(a_box) );
    if (iFluids == 4)
      FORT_CHARGE_EXCHANGE_4F( CHF_FRA(S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_FIA1(Region,0),
                               CHF_CONST_REAL(dt),
                               CHF_CONST_INT(iRhoN),
                               CHF_BOX(a_box) );
            
    a_vars.copy(S,0,isrho,  3);
    a_vars.copy(S,4,isrho+3,1); // copy s-eng, s-momz is skipped
    
    #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(a_vars)
        );
    #endif
  }
  
}

void SWLISMProblemPolar::primForPlot(      FArrayBox& a_W,
                           const Box&       a_box)
{
  //return;
  //m_csh->transCartesianVectToCurv(a_W, a_box, m_level);
  a_W.mult(m_lismN, WRHO);
  a_W.mult(1e-5*m_lismV, WVELX, 3);  
  a_W.mult(m_Bref , WBX  , 3);
  
  int Fluids = nFluids();
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);
  
  for (int iFluid = 1; iFluid < Fluids; ++iFluid)
  { 
    int iRho = iRhoN+(iFluid-1)*WNUM_E;
    a_W.mult(m_lismN, iRho);
    a_W.mult(1e-5*m_lismV, iRho+1, 3);  
  }
  
}


//                            Return boundary condition flags for all boundaries
void SWLISMProblemPolar::getBCFlags( eBoundaryConditions leftBC,
                                eBoundaryConditions rightBC,
                                eBoundaryConditions bottomBC,
                                eBoundaryConditions topBC,
                                eBoundaryConditions frontBC,
                                eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Continuous;
  bottomBC = BC_Axis;
  topBC    = BC_Continuous;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions SWLISMProblemPolar::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  case 1  : return (a_sd == Side::Lo) ? BC_Axis  : BC_Continuous;
  default : return BC_Undefined;
  }
}

/// Creates tagged cells for dynamic mesh refinement
/**
  Problem specific cells tagging
 */
void SWLISMProblemPolar::tagCells(const FArrayBox&  a_U,
                                    const Box&        a_box,
                                          IntVectSet& a_tags)
{    
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();

  /*if (m_level==0)  
  {
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      a_tags |= iv;      
    }
  }*/
  
  if (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym)     
  {
    /*Real PI = 3.141592654;
    Real dphi = m_csh->dx(1,m_level);
    
    IntVect lo,hi;
    lo[0] = (int)(m_domain.domainBox().size(0)/2);
    lo[1] = (int)( (20.0/180.0)*PI/dphi  );    
        
    hi[0] = (int)(m_domain.domainBox().size(0)/1.5);
    hi[1] = (int)( (60.0/180.0)*PI/dphi  );
    
    
    Box b1, bc;
    
    b1.define(lo,hi);
    bc = b1&a_box;
    
    if (!bc.isEmpty())  
    {
      // Tag all bc
      BoxIterator bit(bc);
      for (bit.begin(); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();                  
        a_tags |= iv;      
      }
    }*/
    
    /*Real r0 = 130.0;
    Real r1 = 140.0;
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      if ((m_dudvdw[0]*iv[0] > r0) && (m_dudvdw[0]*iv[0] < r1)) a_tags |= iv;      
    }*/
  }  
}                                         

/// Check geometrical/problem limitations for grid adaptation
/**
 */   
void SWLISMProblemPolar::lockedCellsRegrid( BaseFab<int> & a_flag,
                                 const FArrayBox&  a_U,
                                 const Box&     a_box)
{
  
  if (m_max_levelTS <= m_level )
  {
    FArrayBox W(a_box,m_eqSys->numPrimitives());
    m_eqSys->stateToPrim(W, a_U, a_box);  
                  
    BaseFab<int> Region;
    FArrayBox dummyFab;      
    m_csh->transCartesianVectToCurv(W, a_box, m_level);    
    defineRegions(W, dummyFab, Region, a_box);  
      
    bool Reg2Present = false;
    bool Reg3Present = false;
    int  reg; RealVect cv_coord;
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      m_csh->getCellCenter(cv_coord,iv,m_level);      
      reg = Region(iv,0);
      if (reg == 2) Reg2Present = true;
      if (reg == 3) Reg3Present = true;      
    }
    if (Reg2Present && Reg3Present) a_flag.setVal(1);
  }
      
  if (m_adaptRmin > 0.0)
  {
    RealVect cv_coord;
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      m_csh->getCellCenter(cv_coord,iv,m_level);      
      if (cv_coord[0] < m_adaptRmin) a_flag(iv,0) = 1;     
    }    
  }
  
 
  
  
  // for Medvedev (regridding near the HP)  
  
/*  FArrayBox W(a_box,m_eqSys->numPrimitives());
  m_eqSys->stateToPrim(W, a_U, a_box);  
                
  BaseFab<int> Region;
  FArrayBox dummyFab;      
  m_csh->transCartesianVectToCurv(W, a_box, m_level);    
  defineRegions(W, dummyFab, Region, a_box);  
    
  bool Reg2Present = false;
  bool Reg3Present = false;
  int  reg; RealVect cv_coord;
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    if ((m_level>=1)&&(cv_coord[0]>220.0)&&(cv_coord[1]>d_PI_2)) a_flag(iv,0)=1;
    
    if ((m_level>=2)&&(cv_coord[1]>(160.0/180.0)*d_PI)) a_flag(iv,0)=1;
    
    reg = Region(iv,0);
    if (reg == 2) Reg2Present = true;
    if (reg == 3) Reg3Present = true;      
  }
  if (Reg2Present && Reg3Present) a_flag.setVal(1);*/
}

// Tags cells that do not limit time step  
void SWLISMProblemPolar::lockedCellsTimeStep( BaseFab<int> & a_flag,
                            const FArrayBox&  a_U,
                            const Box&     a_box)
{
  
}

// Converts dimensionless time to time in seconds 
Real SWLISMProblemPolar::getPhysTime(Real a_time)
{
  return (eos_AU/m_lismV)*a_time;
}
