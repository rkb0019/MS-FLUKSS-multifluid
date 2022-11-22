#include <iomanip>
#include <limits>

#ifdef CH_OMPCPP
  #include <omp.h>
#endif


#include "parstream.H"
#include "ParmParse.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"
#include "SPMD.H"
#include "LoadBalance.H"
#include "LevelFluxRegister.H"
#include "ProblemDomain.H"
#include "AMRIO.H"
#include "AMRLevel.H"
#include "computeSum.H"
#include "LoHiCenter.H"
#include "computeNorm.H"

#include "AMRLevelIdealMHD.H"
#include "PatchMHDAM.H"
#include "LGintegrator.H"
#include "PhysProblem.H"
#include "TecplotIO.H"
#include "PatchIdealMHDF_F.H" // needed for FORT_CORRECT_B
#include "PatchMHDAMF_F.H"    // needed for CYL coord system
#include "SourceCalculator.H"
#include "BCFuncMHDAM.H"
#include "printArrayF_F.H"
#include "DivCleaningF_F.H"
#include "EquationSystem.H"
#include "EqSysMHDMF.H"
#include "CSHSphericalF_F.H"
#include "MHDAMDefs.H"
#include "LevelSetF_F.H"
#include "FineInterpMHDAMF_F.H"

extern FILE* extern_tecplotfile1;

// Constructor
AMRLevelIdealMHD::AMRLevelIdealMHD()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD default constructor" << endl;
  }

  m_cfl                     = 0.8;
  m_domainLength            = RealVect(D_DECL(1.0,1.0,1.0)); 
  m_initial_dt_multiplier   = 0.1;
  m_patchMHDAMFactory       = NULL;
  m_patchMHDAM              = NULL;
  m_output_density_gradient = false;
  m_output_B_gradient       = false;  
  m_output_vecCS            = false;
}

// Destructor
AMRLevelIdealMHD::~AMRLevelIdealMHD()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD destructor" << endl;
  }

  // Get rid of the patch integrators (factory and object)
  if (m_patchMHDAMFactory != NULL)
  {
    delete m_patchMHDAMFactory;
  }

  if (m_patchMHDAM != NULL)
  {
    delete m_patchMHDAM;
  }
}


// This instance should never get called - historical
void AMRLevelIdealMHD::define(       AMRLevel*  a_coarserLevelPtr,
                               const Box&       a_problemDomain,
                                     int        a_level,
                                     int        a_refRatio)
{
  ProblemDomain physdomain(a_problemDomain);

  MayDay::Error("AMRLevelIdealMHD::define -\n\tShould never be called with a Box for a problem domain");
}

// Define new AMR level
void AMRLevelIdealMHD::define(       AMRLevel*      a_coarserLevelPtr,
                               const ProblemDomain& a_problemDomain,
                                     int            a_level,
                                     int            a_refRatio)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::define " << a_level << endl;
  }

  // Call inherited define
  AMRLevel::define(a_coarserLevelPtr,
                   a_problemDomain,
                   a_level,
                   a_refRatio);

  // Get setup information from the next coarser level
  if (a_coarserLevelPtr != NULL)
  {
    AMRLevelIdealMHD* amrMHDPtr = dynamic_cast<AMRLevelIdealMHD*>(a_coarserLevelPtr);

    if (amrMHDPtr != NULL)
    {
      m_cfl = amrMHDPtr->m_cfl;
      m_domainLength = amrMHDPtr->m_domainLength;
      m_tagBufferSize = amrMHDPtr->m_tagBufferSize;
  	  m_output_density_gradient=amrMHDPtr->m_output_density_gradient;
      m_output_B_gradient=amrMHDPtr->m_output_B_gradient;
    }
    else
    {
      MayDay::Error("AMRLevelIdealMHD::define: a_coarserLevelPtr is not castable to AMRLevelIdealMHD*");
    }
  }
    

  // Remove old patch integrator (if any), create a new one, and initialize
  if (m_patchMHDAM != NULL)
  {
    delete m_patchMHDAM;
  }

  m_patchMHDAM = m_patchMHDAMFactory->new_patchMHDAM();
  m_patchMHDAM->define(m_problem_domain,a_level,s_verbosity);
  
  m_csh   = m_patchMHDAM->getPhysProblem()->coordinateSystem();
  m_eqSys = m_patchMHDAM->getPhysProblem()->equationSystem();


  // Get additional information from the patch integrator
  m_numStates  = m_eqSys->numStates();
  m_stateNames = m_eqSys->stateNames();  
  m_numGhost   = m_patchMHDAM->numGhostCells();
  
  switch (m_patchMHDAM->timeApproximation())
  {
    case PatchMHDAM::TA_FirstOrder: m_numSCGhost = 0; break;
    case PatchMHDAM::TA_Hancock:    m_numSCGhost = 1; break;
    case PatchMHDAM::TA_RK2:        m_numSCGhost = 2; break;
    default:                        m_numSCGhost = 0; break;
  }
  
  eDivergenceCleaning dc = m_patchMHDAM->getDivergenceCleaning();  
  m_CTused = ((dc == DC_CT_BS)||(dc == DC_CT_GS0)||(dc == DC_CT_GS1));
  m_CPused = (dc == DC_CP);
       
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
  {
  
    Vector<Box> vstrips, vadjacentStrips;
    int i;
    
    Box domBox = m_problem_domain.domainBox();   
    if (m_problem_domain.isPeriodic(1))
      domBox.grow(1,m_numGhost);
    Box strip = domBox;
    
    // theta = 0 bc    
    strip.setBig(2,domBox.smallEnd()[2]+m_numGhost-1);   
    strip.setBig(0,domBox.smallEnd()[0]);
    // each "strip" runs the length of the domain in the j direction
    for(i=0; i<domBox.size(0); i++)
    {
      vstrips.push_back(strip); // these are inside the domain
      strip.shift(2,-m_numGhost); // shift into m_numGhost-ghost cells in k direction
      vadjacentStrips.push_back(strip); // these are on the ghost cells
      strip.shift(2, m_numGhost); // undo
      strip.shift(0,1); // shift to next j strip.
    }
    // theta = PI bc
    strip = domBox;
    strip.setSmall(2,domBox.bigEnd()[2]-m_numGhost+1);     
    strip.setBig(0,domBox.smallEnd()[0]);
    // each "strip" runs the length of the domain in the j direction
    for(i=0; i<domBox.size(0); i++)
    {
      vstrips.push_back(strip); // these are inside the domain
      strip.shift(2, m_numGhost); // shift into m_numGhost-ghost cells in k direction
      vadjacentStrips.push_back(strip); // these are on the ghost cells
      strip.shift(2,-m_numGhost); // undo
      strip.shift(0,1); // shift to next j strip.
    }
    
    Vector<int> procs;
    LoadBalance(procs, vstrips);
    DisjointBoxLayout strips(vstrips, procs/*, m_problem_domain*/);
    DisjointBoxLayout adjacent(vadjacentStrips, procs/*, m_problem_domain*/);        

    m_lstrips.define(strips,   m_numStates);
    m_ladjstr.define(adjacent, m_numStates);        
            
    if (m_hasCoarser)
    {
      /*m_patcherZ0.define(strips0,
                       a_coarserDisjointBoxLayout,
                       m_numStates,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_level-1,
                       m_numGhost,
                       m_csh);
      m_patcherZPI.define(stripsPI,
                       a_coarserDisjointBoxLayout,
                       m_numStates,
                       coarsen(a_domain,a_refineCoarse),
                       a_refineCoarse,
                       m_level-1,
                       m_numGhost,
                       m_csh);*/
    }
  }    

}



// Initialize grids
#if CHOMBO_VERSION_MAJOR < 4         
void AMRLevelIdealMHD::initialGrid(const Vector<Box>& a_newGrids)
#else
void AMRLevelIdealMHD::initialGrid(const DisjointBoxLayout & a_grid)
#endif
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::initialGrid " << m_level << endl;
  }

#if CHOMBO_VERSION_MAJOR < 4         
  // Save original grids and load balance
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);
#else
  m_grids = a_grid;
#endif

  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;
    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  // Define old and new state data structures
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);
  m_UOld.define(m_grids,m_numStates,ivGhost);  
  
  IntVect STGhost = m_numSCGhost*IntVect::Unit;
  if (m_patchMHDAM->getSourceCalculator()!=NULL)
  {
    if (m_patchMHDAM->getSourceCalculator()->numSourceTerms()>0)
    {
      m_SCData.define(m_grids,m_patchMHDAM->getSourceCalculator()->numSourceTerms(), STGhost);
      for(DataIterator dit = m_SCData.dataIterator();dit.ok(); ++dit) m_SCData[dit()].setVal(0.0);
    }
  }

  eDivergenceCleaning dc = m_patchMHDAM->getDivergenceCleaning();
  if (dc == DC_ProjHancock)
  {
    m_divB.define(m_grids,1,IntVect::Zero);
    m_phi.define(m_grids,1,IntVect::Unit);
  }
  if ((dc == DC_CT_BS)||(dc == DC_CT_GS0)||(dc == DC_CT_GS1))
  {
    m_E   .define( m_grids, 1 );
    m_EAcc.define( m_grids, 1 );
  }
  
  if (getFinerLevel()!=NULL)
  {
    m_invVolumes.define(m_grids,1,IntVect::Zero);
  }

  // Set up data structures
  levelSetup();
}

// Initialize data
void AMRLevelIdealMHD::initialData()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::initialData " << m_level << endl;
  }

  PhysProblem* PhysProblemPtr = m_patchMHDAM->getPhysProblem();
  PhysProblemPtr->initialize(m_UNew);

  if (m_patchMHDAM->getSourceCalculator()!=NULL)
    m_patchMHDAM->getSourceCalculator()->initialize(m_SCData);

  //loadfort55();
}

// Things to do after initialization
void AMRLevelIdealMHD::postInitialize()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::postInitialize " << m_level << endl;
  }
}

// Setup menagerie of data structures
void AMRLevelIdealMHD::levelSetup()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::levelSetup " << m_level << endl;
  }

  AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();
  AMRLevelIdealMHD* amrMHDFinerPtr   = getFinerLevel();

  m_hasCoarser = (amrMHDCoarserPtr != NULL);
  m_hasFiner   = (amrMHDFinerPtr   != NULL);

  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();

    m_coarseAverage.define(m_grids,
                           m_numStates,
                           nRefCrse,
                           m_level-1,
                           m_csh);
// Original Chombo class
//    m_coarseAverage.define(m_grids,
//                           m_numStates,
//                           nRefCrse);

    m_coarseAverageEdge.define( m_grids, 1, nRefCrse);

    m_fineInterp.define(m_grids,
                        m_numStates,
                        nRefCrse,                        
                        m_problem_domain,
                        m_level-1,
                        m_csh);
                        
// Original Chombo class                        
//    m_fineInterp.define(m_grids,
//                        m_numStates,
//                        nRefCrse,                        
//                        m_problem_domain);

    if (s_verbosity >= 5)
    {
      pout() << "         m_patcher.define,  m_problem_domain " << m_problem_domain  << " nRefCrse " << nRefCrse << endl;
      //m_grids.physDomain().dumpOn(pout());
      m_grids.physDomain().domainBox().p();
    }
    
    m_patcher.define(m_grids,
                       amrMHDCoarserPtr->m_grids,
                       m_numStates,
                       coarsen(m_problem_domain,nRefCrse),
                       nRefCrse,
                       m_level-1,
                       m_numGhost,
                       m_csh);

    if (( m_patchMHDAM->getPhysProblem()->physicalModel() != PhysProblem::PP_EulerPM ) &&
        ( m_patchMHDAM->getDivergenceCleaning() != 0))
    {
      int iBx  = UBX;
      m_fineInterp.setBXIndex( iBx );
      m_patcher.setBXIndex( iBx );
    }        
                          
    //const DisjointBoxLayout& coarserLevelDomain = amrMHDCoarserPtr->m_grids;

    // Maintain levelMHDAM
    //m_levelMHDAM.define( m_grids,
    //                     coarserLevelDomain,
    //                     m_problem_domain,
    //                     nRefCrse,                         
    //                     m_level,
    //                     m_patchMHDAMFactory,
    //                     m_hasCoarser,
    //                     m_hasFiner,
		//				             s_verbosity);

    // This may look twisted but you have to do this this way because the
    // coarser levels get setup before the finer levels so, since a flux
    // register lives between this level and the next FINER level, the finer
    // level has to do the setup because it is the only one with the
    // information at the time of construction.

    // Maintain flux registers
    amrMHDCoarserPtr->m_fluxRegister.define(m_grids,
                                            amrMHDCoarserPtr->m_grids,
                                            m_problem_domain,
                                            amrMHDCoarserPtr->m_ref_ratio,
                                            m_numStates, m_csh->scaleFineFluxes() );
    amrMHDCoarserPtr->m_fluxRegister.setToZero();
  }
  else
  {
    //m_levelMHDAM.define( m_grids,
    //                     DisjointBoxLayout(),
    //                     m_problem_domain,
    //                     m_ref_ratio,                         
    //                     m_level,
    //                     m_patchMHDAMFactory,
    //                     m_hasCoarser,
    //                     m_hasFiner,
		//				             s_verbosity);
  }
  
  if (m_hasFiner)
  {
    // Compute volumes
    DataIterator dit = m_invVolumes.dataIterator();
    for(;dit.ok(); ++dit)
    {
      FArrayBox& volBox = m_invVolumes[dit()];      
      m_csh->scalingFactor(volBox,volBox.box(),m_level);            
    }
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
  {
  
    Vector<Box> vStrips_0, vStripsPI;    
    Vector<int> procs_0,   procsPI;    
    
    Box domBox = m_problem_domain.domainBox();   
    
    int iLoDomain = domBox.smallEnd()[2];
    int iHiDomain = domBox.bigEnd()[2];        
            
    Box box_0, boxPI; IntVect ivLo, ivHi;
        
    LayoutIterator lit = m_grids.layoutIterator();    
    for (lit.begin();lit.ok();++lit)
    {
      const Box & b_orig = m_grids[lit()];      
      Box b = b_orig;
      bool bLo = (b.smallEnd()[2] == iLoDomain);
      bool bHi = (b.bigEnd()[2]   == iHiDomain);
      
      if ((bLo == false) && (bHi == false)) continue;
      
      int pID = m_grids.procID(lit.i());
      
      if (m_problem_domain.isPeriodic(1))
      {
        if (b.smallEnd()[1] == domBox.smallEnd()[1])
          b.growDir(1,Side::Lo,m_numGhost);
        if (b.bigEnd()[1]   == domBox.bigEnd()[1])      
          b.growDir(1,Side::Hi,m_numGhost);      
      }
      
      //b.grow(0,m_numGhost);
      //b.grow(1,m_numGhost);
      
      if (bLo)
      {
        ivLo = b.smallEnd();ivLo[2] = iLoDomain;
        ivHi = b.bigEnd();  ivHi[2] = iLoDomain + m_numGhost - 1;
        box_0.define(ivLo,ivHi);        
      }
      if (bHi)
      {
        ivLo = b.smallEnd();ivLo[2] = iHiDomain - m_numGhost + 1;
        ivHi = b.bigEnd();  ivHi[2] = iHiDomain;
        boxPI.define(ivLo,ivHi);        
      }
      
      // We do need this if we have only one level
      /*LayoutIterator lit_in = m_grids.layoutIterator();    
      for (lit_in.begin();lit_in.ok();++lit_in)
      {
        const Box& b_in = m_grids[lit_in()];
        
        bool bLo_in = (b_in.smallEnd()[2] == iLoDomain);
        bool bHi_in = (b_in.bigEnd()[2]   == iHiDomain);      
        if ((bLo_in == false) && (bHi_in == false)) continue;        
        if (b_in == b_orig) continue; 
        
        Box b_inter = b & b_in;
        if (b_inter.isEmpty()) continue;
        
        if (bLo && bLo_in)
        {
          for (int dir = 0; dir<2; ++dir)
          if (b_inter.size(dir) == m_numGhost)
          {
            if (b_inter.smallEnd()[dir] == box_0.smallEnd()[dir])
              box_0.setSmall(dir, b_inter.bigEnd()[dir]+1);
            
            if (b_inter.bigEnd()[dir] == box_0.bigEnd()[dir])
              box_0.setBig(dir, b_inter.smallEnd()[dir]-1);            
          }
          CH_assert(!box_0.contains(b_inter));
        }
        if (bHi && bHi_in)
        {
          for (int dir = 0; dir<2; ++dir)
          if (b_inter.size(dir) == m_numGhost)
          {
            if (b_inter.smallEnd()[dir] == boxPI.smallEnd()[dir])
              boxPI.setSmall(dir, b_inter.bigEnd()[dir]+1);
            
            if (b_inter.bigEnd()[dir] == boxPI.bigEnd()[dir])
              boxPI.setBig(dir, b_inter.smallEnd()[dir]-1);            
          }                   
          CH_assert(!boxPI.contains(b_inter));
        }
                
      }*/
      if (bLo) 
      {        
        //box_0 &= m_problem_domain.domainBox();        
        CH_assert(!box_0.isEmpty());
        vStrips_0.push_back(box_0);        
        procs_0.push_back(pID);        
      }
      if (bHi) 
      {       
        //boxPI &= m_problem_domain.domainBox();
        CH_assert(!boxPI.isEmpty());
        vStripsPI.push_back(boxPI);        
        procsPI.push_back(pID);        
      }          
    }
      
        
    DisjointBoxLayout  dbl_0 ( vStrips_0,    procs_0/*, m_problem_domain*/);  
    m_lstrips_0.define(dbl_0, m_numStates);
        
    DisjointBoxLayout  dblPI ( vStripsPI,    procsPI/*, m_problem_domain*/);  
    m_lstripsPI.define(dblPI, m_numStates);    
  }
  
  // Define data for projection scheme
  if (m_patchMHDAM->getDivergenceCleaning() == DC_ProjHancock)        
  {
    ProblemDomain baseDomain(m_problem_domain); // on this level
    int numSolverLevels;
    AMRPoissonOpFactory  PoissonOpFactory;  
    
    PhysProblem* PhysProblemPtr = m_patchMHDAM->getPhysProblem();
    BCHolder bc( RefCountedPtr<BCFuncMHDAM>(new BCFuncMHDAM(m_patchMHDAM->getPhysProblem())));
    
    if (m_hasCoarser)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      numSolverLevels = 2;
      baseDomain.coarsen(nRefCrse);
      Vector<DisjointBoxLayout> allGrids(2);
      allGrids[0] = amrMHDCoarserPtr->m_grids;
      allGrids[1] = m_grids;
      Vector<int> refRatios(1, nRefCrse);
      
      Real dxCrse = m_csh->dx(0,m_level-1);
#if CHOMBO_VERSION_MAJOR < 4                   
      PoissonOpFactory.define(baseDomain,
                                   allGrids,
                                   refRatios,
                                   dxCrse,
                                   bc);
#else
      PoissonOpFactory.define(baseDomain,
                                   allGrids,
                                   refRatios,
                                   dxCrse,
                                   bc,
                                   false // useBoxAgglomeration, I do now know what is this
                                   );
            
#endif
                                   
    } else
    {
#if CHOMBO_VERSION_MAJOR < 4             
        // no coarser level:  define solver on only one level
      numSolverLevels = 1;
      PoissonOpFactory.define(m_problem_domain,
                                   m_grids,
                                   m_csh->dx(0,m_level),
                                   bc);
#else
      // Add code here
#endif
                                   
    }
    
    if ((numSolverLevels == 1) || (m_grids.size()>0)) // Prevent infinite loop when m_grids is empty at the beginning
      m_solver.define(baseDomain, PoissonOpFactory,  &m_bottomSolver, numSolverLevels);
      
    //m_solver.m_eps = 1e-10;    
    //m_solver.m_verbosity = 0;
    //m_bottomSolver.m_verbosity = 0;
    
    
    /*Vector< DisjointBoxLayout > grids;
    Vector< int >               refRatio;
    ProblemDomain               coarsestDomain;
    
    grids.   resize(m_level+1);
    refRatio.resize(m_level+1);
    
    AMRLevelIdealMHD * tmpLevel = this;
    for (int level = m_level;level>=0;--level)
    {
      grids[level] = tmpLevel->m_grids;
      refRatio[level] = tmpLevel->refRatio();
      
      if (level == 0) coarsestDomain = tmpLevel->problemDomain();
      
      tmpLevel = tmpLevel->getCoarserLevel();
    }
    Real coarsestDx = m_csh->dx(0,0);    
    m_opFactory.define(coarsestDomain,
                       grids,
                       refRatio,
                       coarsestDx,
                       bc);                     
    m_solver.define(coarsestDomain, m_opFactory,  &m_bottomSolver, m_level+1);*/
  
  }
  
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::levelSetup " << m_level << "completed" << endl;
  }
  
}

// Return number of points in grid      
int	AMRLevelIdealMHD::numPts()
{
  return (int)(m_grids.numCells());
}
  
// Return the grids at this level 
const DisjointBoxLayout & AMRLevelIdealMHD::boxes()
{
  return m_grids;
}

void AMRLevelIdealMHD::preTimeStep()
{
  if (s_verbosity >= 2)  pout() << "Enter  AMRLevelIdealMHD::preTimeStep"  << endl;
  
  if ((m_finer_level_ptr != NULL) || ( m_CTused == true ))
  {
    DataIterator dit = m_UNew.dataIterator();
    for( ; dit.ok(); ++dit )
    {
      m_UOld[dit()].copy(m_UNew[dit()]);
    }
  }
  
  // A finer level exists
  if( m_hasFiner == true )
  {    
    if( m_CTused == true )
    {
      AMRLevelIdealMHD* pFiner = getFinerLevel();

      for( DataIterator ditF = pFiner->m_EAcc.dataIterator(); ditF.ok(); ++ditF )
      {
        pFiner->m_EAcc[ditF()].setVal( 0.0 );
      }
    }  
    // Recall that my flux register goes between my level and the next
    // finer level      
    m_fluxRegister.setToZero();
  }
  
  // Setup an interval corresponding to the conserved variables
  Interval UInterval(0,m_numStates-1);  
        
  // Fill a_U's ghost cells using fillInterp
  if (m_hasCoarser)
  {
    if (s_verbosity >= 4)  pout() << "Fill a_U's ghost cells using fillInterp "  << endl;
            
    AMRLevelIdealMHD* coarserPtr = getCoarserLevel();
    
    const LevelData<FArrayBox> & UCoarseOld = coarserPtr->m_UOld;
    const LevelData<FArrayBox> & UCoarseNew = coarserPtr->m_UNew;    
    
    Real tCoarserNew = coarserPtr->m_time;
    Real tCoarserOld = tCoarserNew - coarserPtr->m_dt;
    
    int nRefCrse = m_coarser_level_ptr->refRatio();
            
    // Fraction "a_time" falls between the old and the new coarse times
    Real alpha = (m_time - tCoarserOld) / (tCoarserNew - tCoarserOld);

    // Truncate the fraction to the range [0,1] to remove floating-point
    // subtraction roundoff effects
    Real eps = 0.04 * m_dt / nRefCrse;

    if (Abs(alpha) < eps)
      {
        alpha = 0.0;
      }

    if (Abs(1.0-alpha) < eps)
      {
        alpha = 1.0;
      }

    // Current time before old coarse time
    if (alpha < 0.0)
      {
        MayDay::Error( "AMRLevelIdealMHD::step: alpha < 0.0");
      }

    // Current time after new coarse time
    if (alpha > 1.0)
      {
        MayDay::Error( "AMRLevelIdealMHD::step: alpha > 1.0");
      }

    // Interpolate ghost cells from next coarser level using both space
    // and time interpolation
    m_patcher.fillInterp(m_UNew,
                         UCoarseOld,
                         UCoarseNew,
                         alpha,
                         0,0,m_numStates);
                         
    if (m_csh->coordinateSystem() != CoordinateSystemHandler::CS_Spherical) 
    {
      /*m_patcherZ0.fillInterp(m_lstrips0,
                         a_UCoarseOld,
                         a_UCoarseNew,
                         alpha,
                         0,0,m_numStates);
      m_patcherZPI.fillInterp(m_lstripsPI,
                         a_UCoarseOld,
                         a_UCoarseNew,
                         alpha,
                         0,0,m_numStates);*/
    }
    
    if (s_verbosity >= 4)  pout() << "Filling complete"  << endl;
  }
        
  // Exchange all the data between grids at this level  
  m_UNew.exchange(UInterval);  
  
  // z-axis bc
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
  {  
    int sign;
          
    //FArrayBox& curU = m_UNew[m_grids.dataIterator()];
    //FORT_VIEWBOXDATACONST(
    //      CHF_CONST_FRA(curU)
    //      );
    
    int iLoDomain = m_problem_domain.domainBox().smallEnd()[2];
    int iHiDomain = m_problem_domain.domainBox().bigEnd()[2];

    DataIterator ditU   = m_UNew.dataIterator();
    DataIterator dit_0  = m_lstrips_0.dataIterator();  
    DataIterator ditPI  = m_lstripsPI.dataIterator();  
    for(; ditU.ok(); ++ditU)
    {
      const FArrayBox& U = m_UNew[ditU];
      const Box&    boxU = U.box();
      
      bool bLo = (boxU.smallEnd()[2] <= iLoDomain);
      bool bHi = (boxU.bigEnd()[2]   >= iHiDomain);
        
      if ((bLo == false) && (bHi == false)) continue;
      
      if (bLo)
      {
        FArrayBox& from = m_lstrips_0[dit_0];      
        
        CH_assert(boxU.contains(from.box()));
        
        from.copy(U);      
        m_csh->transCartesianVectToCurv(from,from.box(),m_level);    
                     
        ++dit_0;      
      }
      
      if (bHi)
      {
        FArrayBox& from = m_lstripsPI[ditPI];      
        
        CH_assert(boxU.contains(from.box()));
        
        from.copy(U);      
        m_csh->transCartesianVectToCurv(from,from.box(),m_level);    
                     
        ++ditPI;      
      }
      
    }
    
    m_lstrips_0.copyTo(m_lstrips.interval(), m_lstrips, m_lstrips.interval());
    m_lstripsPI.copyTo(m_lstrips.interval(), m_lstrips, m_lstrips.interval());
    
    int kmin = m_problem_domain.domainBox().smallEnd()[2];
    DataIterator dit = m_lstrips.dataIterator();
    DataIterator dit2= m_ladjstr.dataIterator();  
    for(; dit.ok(); ++dit,++dit2)
    {
      FArrayBox& from = m_lstrips[dit];
      FArrayBox& to   = m_ladjstr[dit2];
              
      sign = (from.box().smallEnd()[2] == kmin ? -1 : 1);
      
      // copy with a modular pi shift.
      FORT_ZAXISCOPY(
         CHF_CONST_FRA(from),
         CHF_FRA(to),
         CHF_CONST_INT(m_numGhost),
         CHF_CONST_INT(sign) ); 
    }

    m_ladjstr.copyTo(m_UNew.interval(), m_UNew, m_UNew.interval());
  }
  
  PhysProblem* pPhysPr = m_patchMHDAM->getPhysProblem();    
  pPhysPr->preTimeStep( m_UNew, m_time );  
  
  // Potentially used in boundary conditions
  m_patchMHDAM->setCurrentTime(m_time);
  
  if (s_verbosity >= 2)  pout() << "Leave AMRLevelIdealMHD::preTimeStep"  << endl;
}


Real AMRLevelIdealMHD::stepHancockProjection()
{

// Dummy source used if source term passed in is empty
  FArrayBox dummyFAB;
  
  Real local_dtNew = numeric_limits<Real>::max();
  IntVect minDtCell;

  if (s_verbosity >= 4)  pout() << "Beginning of loop through patches/grids. "  << endl;
  
  timeval time_before,time_after;    
  int    iBox      = 0;
  long   num_cells = 0;
  double stepTime  = 0.0;
      
  LevelData<FluxBox> Bn(m_grids,1,2*IntVect::Unit);
      
  // Predictor
  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit,++iBox)
  {
    // The current box
    Box curBox = m_grids.get(dit());
    
    FluxBox & curBn = Bn[dit()];
    
    pout() << "Box: " << iBox << endl;

#ifndef NDEBUG          
    num_cells+=curBox.volume();
#endif      

  // The current grid of conserved variables
    FArrayBox& curU = m_UNew[dit()];

    // The current source terms if they exist
    const FArrayBox* source = &dummyFAB;
    if (m_SCData.isDefined())
    {
      source = &m_SCData[dit()];
    }
    FArrayBox* divB = &dummyFAB;
    if (m_divB.isDefined())
    {
      divB = &m_divB[dit()];
    }

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FArrayBox flux[SpaceDim];

    EdgeBox E;
                
    if (s_verbosity >= 4)  pout() << "Box: " << iBox << endl;

    // Update the current grid's conserved variables, return the final
    // fluxes used for this, and the maximum wave speed for this grid
    
    // Calculate geometrical info
    FArrayBox scale;                        // Inversion of cell volumes
    FArrayBox areas[CH_SPACEDIM];           // Cell areas
    Box geom_box = curBox;
    geom_box.grow(2); 
    scale.define(geom_box,1);
    m_csh->scalingFactor(scale, geom_box, m_level);
    Box areasBox;
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      areasBox = geom_box;
      areasBox.surroundingNodes(dir);
      areas[dir].define(areasBox,1);
      m_csh->getAreas(areas[dir], areasBox, dir, m_level);
    }
    
    Box UBox = curU.box();
    FArrayBox  W(UBox, m_numStates);  
    FArrayBox  UPredictor(UBox, m_numStates);       
    FArrayBox  WPredictor(UBox, m_numStates);       
    UPredictor.copy( curU );
    
    FArrayBox WMinus[CH_SPACEDIM], WPlus[CH_SPACEDIM];
    
    //scale*=0.5*m_dt;    
    
    m_patchMHDAM->predictorHancock(UPredictor, W, WPredictor, WMinus, WPlus, E, curBn, *divB,  scale, areas, *source, 0.5*m_dt, curBox);      
            
    /*if (iBox == 5)
    {
      FILE* tfile = OpenTecplotFile("divB_box5.dat","TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"divB\"  ");
      WriteFArrayBoxToTecplotFile(tfile, *divB, curBox, Interval(0, 0), m_csh->dx(0,m_level));
      CloseTecplotFile(tfile);
    }*/
              
                                  
  }
  // Projection
  {  
    
#ifndef NDEBUG
    Real sumRHS = computeNorm(m_divB, NULL, 2, m_csh->dx(0,m_level));
    Real maxRHS = computeMax(m_divB, NULL, 2);
    pout() << "Before projection, sum of divB = " << sumRHS << ", max of divB =" << maxRHS << endl;
#endif

    DataIterator dit;

    //solve lapl phi = div(u).
    
    int numSolverLevels = (m_level > 0 ? 2 : 1);
    
    Vector<LevelData<FArrayBox> *> rhsHier(numSolverLevels, NULL), phiHier(numSolverLevels, NULL);    
    if(m_level > 0)
    {        
      rhsHier[1] = &m_divB;
      phiHier[1] = &m_phi;
      phiHier[0] = &(getCoarserLevel()->m_phi);
    } else
    {
      rhsHier[0] = &m_divB;
      phiHier[0] = &m_phi;
    }
    
    iBox = 0;
    for(dit = m_phi.dataIterator(); dit.ok(); ++dit,iBox++)
    {
      m_phi[dit()].setVal(0.0);
      
      FArrayBox & divB = m_divB[dit()];
      Real maxdivB = MAX(fabs(divB.min()), fabs(divB.max()) );
      Real sumdivB = sqrt(divB.sumPow(divB.box(),2));
//      pout() << "Box " << iBox << ", max divB = " << maxdivB << " sqrt(sum(divB**2)) = " << sumdivB << endl;      
    }
          
       
    //PatchMHDMF * pMHDMF = static_cast<PatchMHDMF*>(m_patchMHDAM);
    
    Box fluxBox[SpaceDim];
    
    for (int iProj = 0;iProj<1;++iProj)
    {
      m_solver.solve(phiHier, rhsHier, numSolverLevels-1, numSolverLevels-1);
      // apply appropriate BC's
      /*for (dit.reset(); dit.ok(); ++dit)
      {
        pPhysPr->projectBbc(m_phi[dit],
                            m_grids[dit],
                            m_problem_domain,
                            m_csh->dx(0,m_level),
                            false); // not homogeneous
      }*/
    
      Interval phiComps(0,0);
      m_phi.exchange(phiComps);
      
      pout() << "Projection iteration " << iProj << endl;
    
      iBox = 0;
      for(dit = m_phi.dataIterator(); dit.ok(); ++dit,iBox++)
      {
        Box b = m_grids[dit()];
        
        FluxBox & curBn = Bn[dit()];
        
        FArrayBox & divB = m_divB[dit()];
        
        int dir;
        for(dir = 0; dir < SpaceDim; ++dir )
        {
          fluxBox[dir] = b;
          fluxBox[dir].grow(1);
          fluxBox[dir].grow(dir,-1);
          fluxBox[dir] &= m_problem_domain;
          fluxBox[dir].surroundingNodes(dir);    
        }
        
        for (dir = 0; dir<SpaceDim;++dir)
        m_patchMHDAM->correctBn( curBn[dir], m_phi[dit()], dir, fluxBox[dir]);
        
        m_patchMHDAM->computeDivB( divB, curBn, dummyFAB, b);
        
        Real maxdivB = MAX(fabs(divB.min()), fabs(divB.max()) );
        Real sumdivB = sqrt(divB.sumPow(divB.box(),2));
//        pout() << "Box " << iBox << ", max divB = " << maxdivB << " sqrt(sum(divB**2)) = " << sumdivB << endl;
        
        //m_phi[dit()].setVal(0.0);
      }
    }                    
  }
  iBox=0;
  
  //Corrector  
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit,++iBox)
  {
    // The current box
    Box curBox = m_grids.get(dit());
    
    FluxBox & curBn = Bn[dit()];
    
    pout() << "Box: " << iBox << endl;

#ifndef NDEBUG          
    num_cells+=curBox.volume();
#endif      

  // The current grid of conserved variables
    FArrayBox& curU = m_UNew[dit()];

    // The current source terms if they exist
    const FArrayBox* source = &dummyFAB;
    if (m_SCData.isDefined())
    {
      source = &m_SCData[dit()];
    }
    FArrayBox* divB = &dummyFAB;
    if (m_divB.isDefined())
    {
      divB = &m_divB[dit()];
    }

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FArrayBox flux[SpaceDim];

    EdgeBox E;
    
    Real    local_dtNewGrid;
    IntVect minDtCellGrid;      
    
    if (s_verbosity >= 4)  pout() << "Box: " << iBox << endl;
     
    // Calculate geometrical info
    FArrayBox scale;                        // Inversion of cell volumes
    FArrayBox areas[CH_SPACEDIM];           // Cell areas
    Box geom_box = curBox;
    geom_box.grow(2); 
    scale.define(geom_box,1);
    m_csh->scalingFactor(scale, geom_box, m_level);
    Box areasBox;
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      areasBox = geom_box;
      areasBox.surroundingNodes(dir);
      areas[dir].define(areasBox,1);
      m_csh->getAreas(areas[dir], areasBox, dir, m_level);
    }
    
    Box UBox = curU.box();
    FArrayBox  W(UBox, m_numStates);  
    FArrayBox  UPredictor(UBox, m_numStates);       
    FArrayBox  WPredictor(UBox, m_numStates);       
    FluxBox    dummyBn;
    UPredictor.copy( curU );
    
    FArrayBox WMinus[CH_SPACEDIM], WPlus[CH_SPACEDIM];
    
    /*if (iBox == 5)
    {
      FILE* tfile = OpenTecplotFile("phi_box5.dat","TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"phi\"  ");
      WriteFArrayBoxToTecplotFile(tfile, m_phi[dit()], curBox, Interval(0, 0), m_csh->dx(0,m_level));
      CloseTecplotFile(tfile);
    }*/
    
    //scale*=0.5*m_dt;
    m_patchMHDAM->predictorHancock(UPredictor, W, WPredictor,       WMinus, WPlus, E, dummyBn, *divB, scale, areas, *source, 0.5*m_dt, curBox);
    //scale*=2.0;    
    m_patchMHDAM->correctorHancock(curU,       W, WPredictor, flux, WMinus, WPlus, E, curBn,   *divB, scale, areas, *source,     m_dt, curBox);
        
    // Get next time step for this patch
    local_dtNewGrid = m_patchMHDAM->computeDt( curU, curBox, minDtCellGrid);     
        
    if (local_dtNewGrid < local_dtNew)  
    {
      local_dtNew = local_dtNewGrid;
      minDtCell   = minDtCellGrid;        
    }
      
    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++)
    {                          
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner)
        {
          FArrayBox fluxC(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementCoarse(fluxC, flux[idir], idir, m_level);
          m_fluxRegister.incrementCoarse(fluxC,m_dt,dit(),
                                              m_eqSys->consInterval(),
                                              m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_finerFluxRegister.poutCoarseRegisters();
        }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser)
        {
          AMRLevelIdealMHD* coarserPtr = getCoarserLevel();    
          LevelFluxRegister & coarserFR = coarserPtr->m_fluxRegister;
          FArrayBox fluxM(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementFine(fluxM, flux[idir], idir, m_level);
          coarserFR.incrementFine(fluxM,m_dt,dit(),
                                             m_eqSys->consInterval(),
                                             m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_coarserFluxRegister.poutFineRegisters();      
        }
    }
  }
  
#ifndef NDEBUG
    Real sumRHS = computeNorm(m_divB, NULL, 2, m_csh->dx(0,m_level));
    Real maxRHS = computeMax(m_divB, NULL, 2);
    pout() << "After updatestate, sum of divB = " << sumRHS << ", max of divB =" << maxRHS << endl;
#endif
                
  //pout() << "boxes: " << iBox <<", cells: " << num_cells << endl;  
  
  // Return the maximum stable time step
  return local_dtNew;
      
}

Real AMRLevelIdealMHD::step(IntVect &  a_minDtCell)
{
  // Dummy source used if source term passed in is empty
  FArrayBox dummyFAB;
  
  Real local_dtNew = numeric_limits<Real>::max();  
    

  if (s_verbosity >= 2)  pout() << "Beginning of loop through patches/grids. "  << endl;
  
  timeval time_before,time_after;    
  int    iBox      = 0;
  long   num_cells = 0;
  double stepTime  = 0.0;
  
  PhysProblem* pPhysPr = m_patchMHDAM->getPhysProblem();
      
  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit,++iBox)
  {
    // The current box
    Box curBox = m_grids.get(dit());

#ifndef NDEBUG          
    num_cells+=curBox.volume();
#endif      

  // The current grid of conserved variables
    FArrayBox& curU = m_UNew[dit()];

    // External source terms if they exist
    const FArrayBox* source = &dummyFAB;
    if (m_SCData.isDefined())
    {
      source = &m_SCData[dit()];
    }
    FArrayBox* divB = &dummyFAB;
    if (m_divB.isDefined())
    {
      divB = &m_divB[dit()];
    }

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FArrayBox flux[SpaceDim];

    EdgeBox E;

    if( m_CTused )
    {
      E.define( curU.box() );
    }

    Real    local_dtNewGrid;
    IntVect minDtCellGrid;      
    
    if (s_verbosity >= 4)  pout() << "Box: " << iBox << endl;

    // Update the current grid's conserved variables, return the final
    // fluxes used for this, and the minimum time step for this patch
        
    gettimeofday(&time_before,NULL);    
  
    m_patchMHDAM->updateState(curU,
                              flux,
                              E,
                              local_dtNewGrid,
                              minDtCellGrid,
                              *source,
                              *divB,
                              m_dt,
                              curBox);
                              
    gettimeofday(&time_after, NULL);
    stepTime += (double)(time_after.tv_sec-time_before.tv_sec) + 1e-6*((double)(time_after.tv_usec-time_before.tv_usec));      
                                        
    if( m_CTused )
    {
      EdgeBox & curE = m_E[dit()];
      Interval EInt( 0, 0 );
      curE.copy( curBox, EInt, E, EInt ); 
    }

    if (local_dtNewGrid < local_dtNew)  
    {
      local_dtNew = local_dtNewGrid;
      a_minDtCell   = minDtCellGrid;        
    }
          

    // Do flux register updates
    if (m_patchMHDAM->getLSonlyFlag() == false)
    for (int idir = 0; idir < SpaceDim; idir++)
    {                          
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner)
        {
          FArrayBox fluxC(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementCoarse(fluxC, flux[idir], idir, m_level);
          m_fluxRegister.incrementCoarse(fluxC,m_dt,dit(),
                                              m_eqSys->consInterval(),
                                              m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_finerFluxRegister.poutCoarseRegisters();
        }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser)
        {
          AMRLevelIdealMHD* coarserPtr = getCoarserLevel();    
          LevelFluxRegister & coarserFR = coarserPtr->m_fluxRegister;
          FArrayBox fluxM(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementFine(fluxM, flux[idir], idir, m_level);
          coarserFR.incrementFine(fluxM,m_dt,dit(),
                                             m_eqSys->consInterval(),
                                             m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_coarserFluxRegister.poutFineRegisters();      
        }
    }
  }
          
    
  
  //pout() << "boxes: " << iBox <<", cells: " << num_cells << endl;  
      
  
//#ifndef NDEBUG
//    Real sumRHS = computeNorm(m_divB, NULL, 2, m_csh->dx(0,m_level));
//    Real maxRHS = computeMax(m_divB, NULL, 2);
//    pout() << "sum of divB = " << sumRHS << ", max of divB =" << maxRHS << endl;
//#endif
  

  // Return the maximum stable time step
  return local_dtNew;
}


#ifdef CH_OMPCPP
Real AMRLevelIdealMHD::stepOMP(IntVect &  a_minDtCell)
{

  Real local_dtNew = numeric_limits<Real>::max();  
#pragma omp parallel
{
  // Dummy source used if source term passed in is empty
  
  int num_threads = omp_get_num_threads();
  int thread_num =  omp_get_thread_num();
  
  FArrayBox dummyFAB;
  
  Real local_dtNewThread = numeric_limits<Real>::max();  
    

  if (s_verbosity >= 4)  pout() << "Beginning of loop through patches/grids. "  << endl;
  
  timeval time_before,time_after;    
  int    iBox      = 0;
  long   num_cells = 0;
  double stepTime  = 0.0;
  
  PhysProblem* pPhysPr = m_patchMHDAM->getPhysProblem();
      
  // Beginning of loop through patches/grids.
  for (DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit,++iBox)
  if (iBox%num_threads == thread_num)
  {
    // The current box
    Box curBox = m_grids.get(dit());

#ifndef NDEBUG          
    num_cells+=curBox.volume();
#endif      

  // The current grid of conserved variables
    FArrayBox& curU = m_UNew[dit()];

    // The current source terms if they exist
    const FArrayBox* source = &dummyFAB;
    if (m_SCData.isDefined())
    {
      source = &m_SCData[dit()];
    }
    FArrayBox* divB = &dummyFAB;
    if (m_divB.isDefined())
    {
      divB = &m_divB[dit()];
    }

    // The fluxes computed for this grid - used for refluxing and returning
    // other face centered quantities
    FArrayBox flux[SpaceDim];

    EdgeBox E;

    if( m_CTused )
    {
      E.define( curU.box() );
    }

    Real    local_dtNewGrid;
    IntVect minDtCellGrid;      
    
    if (s_verbosity >= 4)  pout() << "Box: " << iBox << endl;

    // Update the current grid's conserved variables, return the final
    // fluxes used for this, and the maximum wave speed for this grid  
        
    gettimeofday(&time_before,NULL);    
  
    m_patchMHDAM->updateState(curU,
                              flux,
                              E,
                              local_dtNewGrid,
                              minDtCellGrid,
                              *source,
                              *divB,
                              m_dt,
                              curBox);
                              
    gettimeofday(&time_after, NULL);
    stepTime += (double)(time_after.tv_sec-time_before.tv_sec) + 1e-6*((double)(time_after.tv_usec-time_before.tv_usec));      
                                        
    if( m_CTused )
    {
      EdgeBox & curE = m_E[dit()];
      Interval EInt( 0, 0 );
      curE.copy( curBox, EInt, E, EInt ); 
    }

    if (local_dtNewGrid < local_dtNewThread)  
    {
      local_dtNewThread = local_dtNewGrid;
      a_minDtCell   = minDtCellGrid;        
    }
          

    // Do flux register updates
    for (int idir = 0; idir < SpaceDim; idir++)
    {                          
      // Increment coarse flux register between this level and the next
      // finer level - this level is the next coarser level with respect
      // to the next finer level
      if (m_hasFiner)
        {
          FArrayBox fluxC(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementCoarse(fluxC, flux[idir], idir, m_level);
          m_fluxRegister.incrementCoarse(fluxC,m_dt,dit(),
                                              m_eqSys->consInterval(),
                                              m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_finerFluxRegister.poutCoarseRegisters();
        }

      // Increment fine flux registers between this level and the next
      // coarser level - this level is the next finer level with respect
      // to the next coarser level
      if (m_hasCoarser)
        {
          AMRLevelIdealMHD* coarserPtr = getCoarserLevel();    
          LevelFluxRegister & coarserFR = coarserPtr->m_fluxRegister;
          FArrayBox fluxM(flux[idir].box(),flux[idir].nComp());
          m_csh->modifyFluxesFRincrementFine(fluxM, flux[idir], idir, m_level);
          coarserFR.incrementFine(fluxM,m_dt,dit(),
                                             m_eqSys->consInterval(),
                                             m_eqSys->consInterval(),idir);
          //pout() << "LevelMHDAM::step, Level " << m_level << endl;
          //a_coarserFluxRegister.poutFineRegisters();      
        }
    }
  }
          
  #pragma omp critical
  {
    if (local_dtNewThread < local_dtNew)  
      {
        local_dtNew = local_dtNewThread;
      }
  } 
} // #pragma omp parallel 
  
  //pout() << "boxes: " << iBox <<", cells: " << num_cells << endl;  
      
  
//#ifndef NDEBUG
//    Real sumRHS = computeNorm(m_divB, NULL, 2, m_csh->dx(0,m_level));
//    Real maxRHS = computeMax(m_divB, NULL, 2);
//    pout() << "sum of divB = " << sumRHS << ", max of divB =" << maxRHS << endl;
//#endif
  

  // Return the maximum stable time step
  return local_dtNew;
}
#endif

// Advance by one timestep
Real AMRLevelIdealMHD::advance()
{
  if ((s_verbosity >= 2) || ((s_verbosity <= -2)&&(procID() == 0)))
  {
    pout() << "  AMRLevelIdealMHD::advance level " << m_level << " to time " << m_time + m_dt << " (dt=" << m_dt << ")" << endl;
  }
  
  preTimeStep();
  
  Real local_dtNew = 0.0;
  IntVect minDtCell;
          
  if (m_patchMHDAM->getDivergenceCleaning() == DC_ProjHancock)        
    local_dtNew = stepHancockProjection();
   else      
   {
#ifndef CH_OMPCPP
//    local_dtNew = step(minDtCell); 
#else    
//    local_dtNew = stepOMP(minDtCell); 
#endif    
    local_dtNew = step(minDtCell); 
   }
    
  m_patchMHDAM->postTimeStep();
  PhysProblem* pPhysPr = m_patchMHDAM->getPhysProblem();    
  pPhysPr->postTimeStep( m_UNew );
    
  
// Find the minimum of dt's over this level
  
  //local_dtNew = 0.04/0.2; When uncomment this, modify the same line in LevelMHDAM::computeDt
  Real dtNew;

#ifdef CH_MPI
  int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL, MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){ //bark!!!
	MayDay::Error("sorry, but I had a communcation error on new dt");
  }
#else
  dtNew = local_dtNew;
#endif

  if ((m_level == 0) && (m_CPused == true))
  {
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
    CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();        
    Real Ch,Cp,factorCh,factorCp;    
    eqSys->getCorrectionPotentialParams(factorCh,factorCp);
    if ((CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
        (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric) ||
        (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) )
    
    {
      Real dx = m_csh->dx(0,m_level);
      Ch = factorCh*dx/dtNew;        
      Cp = sqrt(Ch*dx/factorCp);
    }
    if (procID() == 0) pout() << "  Ch = " << Ch << "  Cp = " << Cp << endl;            
    eqSys->setCorrectionPotential(Ch,Cp);    
  }
  
  if (abs(s_verbosity) >= 2)  
  {        
#ifdef CH_MPI    
    int cpuID = numProc()+1;
    if (local_dtNew == dtNew) cpuID = procID();
    
    int mindt_cpu = 0;
    int result = MPI_Allreduce(&cpuID, &mindt_cpu, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
    if (mindt_cpu!=0)
    {
      if (procID() == mindt_cpu)             
        MPI_Send(minDtCell.dataPtr(), CH_SPACEDIM, MPI_INT, 0, 0, Chombo_MPI::comm);
      if (procID() == 0)             
        MPI_Recv(minDtCell.dataPtr(), CH_SPACEDIM, MPI_INT, mindt_cpu, 0, Chombo_MPI::comm, MPI_STATUS_IGNORE);
    }
#endif   

    if (procID() == 0)
    {      
      //long pos = pout().tellp();
      //pos--;
      //pout().seekp(pos);
            
      pout() << " new dt = "  << m_cfl*dtNew << " at "; minDtCell.p(); 
      //pout() << endl; 
                  
      //pos = pout().tellp();
      //pos--;
      //pout().seekp(pos);      
      //char buf[20];      
      //sprintf(buf,"%.4g",stepTime);
      //pout()  << " (" << buf << " sec)" << endl;
    }
  }  
      
    
  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * dtNew;

  m_dtNew = returnDt;
  
  if (s_verbosity >= 2) pout() << "Leave AMRLevelIdealMHD::advance "  << endl;     

  return returnDt;
}


// Things to do after a timestep
void AMRLevelIdealMHD::postTimeStep()
{
  // Used for conservation tests
  static Real orig_integral = 0.0;
  static Real last_integral = 0.0;
  static bool first = true;

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::postTimeStep " << m_level << endl;
  }
  
  /*if (m_level == 0)
  { 
    pout() << "Level " << m_level << endl;
    m_fluxRegister.poutCoarseRegisters();
    m_fluxRegister.poutFineRegisters();      
    //pout() << "Level " << m_level << " dir=1 " << endl;
    //m_fluxRegister[1].poutCoarseRegisters();
    //m_fluxRegister[1].poutFineRegisters();      
  }
  if (m_level == 1)
  {
    pout() << "Level " << m_level << " dir=0 " <<  endl;
    getCoarserLevel()->m_fluxRegister[0].poutCoarseRegisters();
    getCoarserLevel()->m_fluxRegister[0].poutFineRegisters();      
    //pout() << "Level " << m_level << " dir=1 " <<  endl;
    //getCoarserLevel()->m_fluxRegister[1].poutCoarseRegisters();
    //getCoarserLevel()->m_fluxRegister[1].poutFineRegisters();      
  }*/



  if (m_hasFiner)
  {
    DataIterator dit;
    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
    
    // Reflux    
    if (m_patchMHDAM->getLSonlyFlag() == false)
    {
      Interval fluxInt = m_eqSys->consInterval();    
      m_fluxRegister.reflux(m_UNew, -1.0, fluxInt, fluxInt, m_invVolumes);              
    }
        
    // Average from finer level data
    AMRLevelIdealMHD* amrMHDFinerPtr = getFinerLevel();
    amrMHDFinerPtr->m_coarseAverage.averageToCoarse( m_UNew, amrMHDFinerPtr->m_UNew );
  }

  if( m_CTused == true )
  {
    if( m_hasFiner == true )
    {
      AMRLevelIdealMHD* amrMHDFinerPtr = getFinerLevel();

      double inv_dt  = 1.0/m_dt;
      amrMHDFinerPtr->m_coarseAverageEdge.averageToCoarse( m_E, amrMHDFinerPtr->m_EAcc, inv_dt );
    }

    int iBX = UBX;

    for( DataIterator dit = m_grids.dataIterator(); dit.ok(); ++dit )
    {
      Box curBox = m_grids.get( dit() );

      FArrayBox& U    = m_UNew[dit()];
      FArrayBox& UOld = m_UOld[dit()];

      EdgeBox& E    = m_E[dit()];

///   FORT_CHECKDIVB( CHF_CONST_FRA(U), CHF_CONST_INT(iBX), CHF_BOX(curBox) );

      m_patchMHDAM->recalculateMagneticField( U, UOld, E, m_dt, curBox );

///   FORT_CHECKDIVB( CHF_CONST_FRA(U), CHF_CONST_INT(iBX), CHF_BOX(curBox) );

    }

    if( m_hasCoarser == true )
    {
      for( DataIterator dit = m_E.dataIterator(); dit.ok(); ++dit )
      {
        EdgeBox& E    = m_E   [dit()];
        EdgeBox& EAcc = m_EAcc[dit()];

        E    *= m_dt;
        EAcc += E;
      }
    }
  }

  if (s_verbosity >= 4 && m_level == 0)
  {
    int nRefFine = 1;

    pout() << "AMRLevelIdealMHD::postTimeStep:" << endl;
    pout() << "  Sums:" << endl;

    if (0)
    for (int comp = 0; comp < m_numStates; comp++)
    {
      Real dx = 1.0;
      Interval curComp(comp,comp);
      Real integral = computeSum(m_UNew,NULL,nRefFine,dx,curComp);

      pout() << "    " << std::setw(23)
                       << std::setprecision(16)
                       << std::setiosflags(std::ios::showpoint)
                       << std::setiosflags(std::ios::scientific)
                       << integral
             << " --- " << m_stateNames[comp];

      if (comp == 0 && !first) {
        pout() << " (" << std::setw(23)
                       << std::setprecision(16)
                       << std::setiosflags(std::ios::showpoint)
                       << std::setiosflags(std::ios::scientific)
                       << (integral-last_integral)
               << " " << std::setw(23)
                      << std::setprecision(16)
                      << std::setiosflags(std::ios::showpoint)
                      << std::setiosflags(std::ios::scientific)
                      << (integral-orig_integral)
               << ")";
      }

      pout() << endl;

      if (comp == 0)
      {
        if (first)
        {
          orig_integral = integral;
          first = false;
        }

        last_integral = integral;
      }
    }
  }

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::postTimeStep " << m_level << " finished" << endl;
  }
}

void AMRLevelIdealMHD::fillGhostCellsCons(Interval & a_interval)
{  
  Interval ival = a_interval;
  if (a_interval.size() == 0)
    ival = m_eqSys->consInterval();
    
  if (ival.size() > 0)
  {
    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
     
    // If there is a coarser level interpolate undefined ghost cells
    if (m_hasCoarser)
    { 
      AMRLevelIdealMHD* coarserPtr = getCoarserLevel();
      const LevelData<FArrayBox> & UCoarseNew = coarserPtr->m_UNew;    
          
      PiecewiseLinearFillPatchMHDAM pwl_main(levelDomain,
                               coarserPtr->m_UNew.disjointBoxLayout(),
                               m_numStates,
                               coarserPtr->m_problem_domain,
                               coarserPtr->m_ref_ratio,
                               m_level-1,
                               2,m_csh);

      pwl_main.fillInterp(m_UNew, UCoarseNew, UCoarseNew, 1.0,
          ival.begin(), ival.begin(), ival.size());
    }

    m_UNew.exchange(ival);
  }
      

}

// Create tags for regridding
void AMRLevelIdealMHD::tagCells(IntVectSet& a_tags)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::tagCells " << m_level << endl;
  }

  // Since tags are calculated using only current time step data, use
  // the same tagging function for initialization and for regridding.
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::tagCellsInit " << m_level << endl;
  }

  Interval InterpInterval;
  m_patchMHDAM->tagCellVarsInterval(InterpInterval);

  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();

  if (InterpInterval.size() > 0)
  {
    // If there is a coarser level interpolate undefined ghost cells
    if (m_hasCoarser)
    {
      const AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();

      PiecewiseLinearFillPatchMHDAM pwl_main(levelDomain,
                               amrMHDCoarserPtr->m_UNew.disjointBoxLayout(),
                               m_numStates,
                               amrMHDCoarserPtr->m_problem_domain,
                               amrMHDCoarserPtr->m_ref_ratio,
                               m_level-1,
                               1,m_csh);


      pwl_main.fillInterp(m_UNew, amrMHDCoarserPtr->m_UNew, amrMHDCoarserPtr->m_UNew, 1.0,
          InterpInterval.begin(), InterpInterval.begin(), InterpInterval.size());
    }

    m_UNew.exchange(InterpInterval);
  }

  IntVectSet localTags;

  if (s_verbosity >= 3)
  pout() << "Start: Compute relative gradient" << endl;


  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    const FArrayBox& UFab = m_UNew[dit()];

    m_patchMHDAM->tagCells(UFab, b, localTags);

	
  }
  if (s_verbosity >= 3)
  pout() << "End: Compute relative gradient"<< endl;


  localTags.grow(m_tagBufferSize);

  // Need to do this in two steps unless a IntVectSet::operator &=
  // (ProblemDomain) operator is defined
  Box localTagsBox = localTags.minBox();
  localTagsBox &= m_problem_domain;
  localTags &= localTagsBox;

  a_tags = localTags;

  /*
  Vector<IntVectSet> allTags;

  const int destProc = uniqueProc(SerialTask::compute);

  gather(allTags,localTags,destProc);

  if (procID() == uniqueProc(SerialTask::compute))
  {
    for (int i = 0; i < allTags.size(); ++i)
    {
      a_tags |= allTags[i];
    }
  }
  */
}

// Create tags at initialization
void AMRLevelIdealMHD::tagCellsInit(IntVectSet& a_tags)
{
  tagCells(a_tags);
}


void AMRLevelIdealMHD::reinitializeLevelSet()
{

  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::reinitializeLevelSet " << m_level << endl;
  }
  

  Interval InterpInterval = m_eqSys->lvlsStateInterval();  

  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();

  if (InterpInterval.size() > 0)
  {
    // If there is a coarser level interpolate undefined ghost cells
    if (m_hasCoarser)
    {
      const AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();

      PiecewiseLinearFillPatchMHDAM pwl_main(levelDomain,
                               amrMHDCoarserPtr->m_UNew.disjointBoxLayout(),
                               m_numStates,
                               amrMHDCoarserPtr->m_problem_domain,
                               amrMHDCoarserPtr->m_ref_ratio,
                               m_level-1,
                               1,m_csh);


      pwl_main.fillInterp(m_UNew, amrMHDCoarserPtr->m_UNew, amrMHDCoarserPtr->m_UNew, 1.0,
          InterpInterval.begin(), InterpInterval.begin(), InterpInterval.size());
    }

    m_UNew.exchange(InterpInterval);
  }
  

  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b = levelDomain[dit()];
    FArrayBox& UFab = m_UNew[dit()];
    
    for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
    {
      int ils  = m_eqSys->lsIndexCons(i);
    
      FORT_LSREINITIALIZE(CHF_FRA(UFab),               
        CHF_CONST_INT(ils),     
        CHF_BOX(b));        	
    }
  }
  
}

// Set up data on this level after regridding
#if CHOMBO_VERSION_MAJOR < 4         
void AMRLevelIdealMHD::regrid(const Vector<Box>& a_newGrids)
#else
void AMRLevelIdealMHD::regrid(const DisjointBoxLayout & a_grid)
#endif
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::regrid " << m_level << endl;
  }

  //pout() << "AMRLevelIdealMHD::regrid " << m_level << endl;

  // Save original grids and load balance
#if CHOMBO_VERSION_MAJOR < 4           
  m_level_grids = a_newGrids;
  m_grids = loadBalance(a_newGrids);
#else
  m_grids = a_grid;
#endif  


  if (s_verbosity >= 4)
  {
    // Indicate/guarantee that the indexing below is only for reading
    // otherwise an error/assertion failure occurs
    const DisjointBoxLayout& constGrids = m_grids;

    pout() << "new grids: " << endl;

    for (LayoutIterator lit = constGrids.layoutIterator(); lit.ok(); ++lit)
    {
      pout() << constGrids[lit()] << endl;
    }
  }

  bool bSTExist = false;  // if external source terms exist in m_SCData
  if (m_patchMHDAM->getSourceCalculator()!=NULL)
  if (m_patchMHDAM->getSourceCalculator()->numSourceTerms()>0) bSTExist = true;

  // Save data for later

  DataIterator dit = m_UNew.dataIterator();
  for(;dit.ok(); ++dit)
  {
    m_UOld[dit()].copy(m_UNew[dit()]);
  }

  LevelData<FArrayBox> SourceTermsOld;
  if (bSTExist) SourceTermsOld.define(m_SCData);

  // Reshape state with new grids
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  m_UNew.define(m_grids,m_numStates,ivGhost);

  if (bSTExist)
  {
    IntVect STGhost = m_numSCGhost*IntVect::Unit;
    m_SCData.define(m_grids,m_patchMHDAM->getSourceCalculator()->numSourceTerms(), STGhost);
  }
  
  eDivergenceCleaning dc = m_patchMHDAM->getDivergenceCleaning();
  if (dc == DC_ProjHancock)
  {
    m_divB.define(m_grids,1,IntVect::Zero);
    m_phi.define(m_grids,1,IntVect::Unit);
  }
  if ((dc == DC_CT_BS)||(dc == DC_CT_GS0)||(dc == DC_CT_GS1))
  {
    m_E   .define( m_grids, 1 );
    m_EAcc.define( m_grids, 1 );
  }
  
  if (getFinerLevel()!=NULL)
  {
    m_invVolumes.define(m_grids,1,IntVect::Zero);
  }
    


  // Set up data structures
  levelSetup();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();

    m_fineInterp.interpToFine(m_UNew,amrMHDCoarserPtr->m_UNew);
    if (bSTExist)
    {
      int nRefCrse = m_coarser_level_ptr->refRatio();
      FineInterp fineInterpST;
      fineInterpST.define(m_grids,
                        m_SCData.nComp(),
                        nRefCrse,
                        m_problem_domain);
      fineInterpST.interpToFine(m_SCData,amrMHDCoarserPtr->m_SCData);
    }
  }

  // Copy from old state
  m_UOld.copyTo( m_UOld.interval(), m_UNew, m_UNew.interval() );

  m_UOld.define(m_grids,m_numStates,ivGhost);
  

  if (bSTExist)
    SourceTermsOld.copyTo(SourceTermsOld.interval(),
				m_SCData,
				m_SCData.interval());

}

void AMRLevelIdealMHD::postRegrid(int a_base_level)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::postRegrid " << m_level << endl;
  }
  
  FArrayBox dummyW;
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& CurBox = levelDomain[dit()];
    FArrayBox& CurU = m_UNew[dit()];

    m_patchMHDAM->getPhysProblem()->postprocessing( CurU, dummyW, dt(), time(), CurBox );
   // m_patchMHDAM->postprocessing(CurU,CurBox);
  }
  
  Real dtNew = computeNewDt(m_cfl);
    
  if (s_verbosity >= 3)
  {
    pout() << "postRegrid dt=" << dtNew << endl;
  }
  //m_dtNew = dtNew;
  //dt(dtNew);

  // Filling ghost cells of m_SCData
  // We can not do it in regrid because coarser levels are regrided later then fine levels.  
  if (m_SCData.isDefined() )
  {
    AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();
    int numGhost = m_SCData.ghostVect()[0];
    if (m_hasCoarser)
    {
      const DisjointBoxLayout& levelDomain = m_SCData.disjointBoxLayout();

      PiecewiseLinearFillPatchMHDAM pwl(levelDomain,
                               amrMHDCoarserPtr->m_SCData.disjointBoxLayout(),
                               m_SCData.nComp(),
                               amrMHDCoarserPtr->m_problem_domain,
                               amrMHDCoarserPtr->m_ref_ratio,
                               m_level-1,
                               numGhost,
                               m_csh);

      if ( (m_patchMHDAM->getPhysProblem()->physicalModel() != PhysProblem::PP_EulerPM ) &&
           (m_patchMHDAM->getDivergenceCleaning() != 0))
      {
        int iBx  = UBX;
        pwl.setBXIndex( iBx );
      }

      pwl.fillInterp( m_SCData, amrMHDCoarserPtr->m_SCData, amrMHDCoarserPtr->m_SCData,
                      1.0, 0, 0, m_SCData.nComp() );
    }

    m_SCData.exchange(Interval(0,m_SCData.nComp()-1));
  }
}


PatchMHDAM* AMRLevelIdealMHD::getpatchMHDAM()
{
   return m_patchMHDAM;
}


// Create a local PatchMHDAM factory using the argument as a factory
void AMRLevelIdealMHD::patchMHDAM(const PatchMHDAM* const a_patchMHDAM)
{
  if (m_patchMHDAMFactory != NULL)
  {
    delete m_patchMHDAMFactory;
  }

  m_patchMHDAMFactory = a_patchMHDAM->new_patchMHDAM();
}

Real AMRLevelIdealMHD::computeNewDt(Real a_cfl)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::computeNewDt " << m_level << endl;
  }
  
  const DisjointBoxLayout& disjointBoxLayout = m_UNew.disjointBoxLayout();
  DataIterator dit = disjointBoxLayout.dataIterator();

  Real local_dtNew = numeric_limits<Real>::max();
  IntVect minDtCell;

  // This computation doesn't need involve a time but the time being set
  // is checked by PatchMHDAM::getMaxWaveSpeed so we have to set it
  // to something...
  m_patchMHDAM->setCurrentTime(0.0);
  
  int iBox = 0;
  
  
  // Loop over all grids to get the maximum wave speed
  for (dit.begin(); dit.ok(); ++dit,++iBox)
  {
    Real    local_dtNewGrid;
    IntVect minDtCellGrid;      
      
    const Box& currentBox = disjointBoxLayout.get(dit());

    // Set the current box and get the maximum wave speed on the current grid    
    local_dtNewGrid = m_patchMHDAM->computeDt(m_UNew[dit()], currentBox, minDtCellGrid);
    
    if (local_dtNewGrid < local_dtNew)  
    {
      local_dtNew = local_dtNewGrid;
      minDtCell   = minDtCellGrid;        
    }        
  }

  //local_dtNew = 0.04/0.2;
  Real dtNew;

#ifdef CH_MPI
  int result = MPI_Allreduce(&local_dtNew, &dtNew, 1, MPI_CH_REAL,
                                 MPI_MIN, Chombo_MPI::comm);
  if(result != MPI_SUCCESS){ //bark!!!
    MayDay::Error("sorry, but I had a communcation error on new dt");
  }
#else
  dtNew = local_dtNew;
#endif

  if ((m_level == 0) && (m_CPused == true))
  {
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);        
    CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
    Real Ch,Cp,factorCh,factorCp;    
    eqSys->getCorrectionPotentialParams(factorCh,factorCp);
    if ((CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
        (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric) ||
        (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) )
    {
      Real dx = m_csh->dx(0,m_level);
      Ch = factorCh*dx/dtNew;        
      Cp = sqrt(Ch*dx/factorCp);
    }
    if (procID() == 0) pout() << "  Ch = " << Ch << "  Cp = " << Cp << endl;            
    eqSys->setCorrectionPotential(Ch,Cp);    
  }

  
  if (abs(s_verbosity) >= 2)  
  {        
#ifdef CH_MPI    
    int cpuID = numProc()+1;
    if (local_dtNew == dtNew) cpuID = procID();
    
    int mindt_cpu = 0;
    int result = MPI_Allreduce(&cpuID, &mindt_cpu, 1, MPI_INT, MPI_MIN, Chombo_MPI::comm);
    if (mindt_cpu!=0)
    {
      if (procID() == mindt_cpu)             
        MPI_Send(minDtCell.dataPtr(), CH_SPACEDIM, MPI_INT, 0, 0, Chombo_MPI::comm);
      if (procID() == 0)             
        MPI_Recv(minDtCell.dataPtr(), CH_SPACEDIM, MPI_INT, mindt_cpu, 0, Chombo_MPI::comm, MPI_STATUS_IGNORE);
    }
#endif   

    if (procID() == 0)
    {      
      //long pos = pout().tellp();
      //pos--;
      //pout().seekp(pos);
            
      pout() << "Level " << m_level << "  dt = "  << a_cfl*dtNew << " at "; minDtCell.p(); 
      //pout() << endl; 
                  
      //pos = pout().tellp();
      //pos--;
      //pout().seekp(pos);      
      //char buf[20];      
      //sprintf(buf,"%.4g",stepTime);
      //pout()  << " (" << buf << " sec)" << endl;
    }
  }  
  
  

  dtNew *= a_cfl;
  
  return dtNew;  
}

// Compute dt using initial data
Real AMRLevelIdealMHD::computeInitialDt()
{
  return computeNewDt(m_initial_dt_multiplier);
}

// Returns the dt computed earlier for this level
Real AMRLevelIdealMHD::computeDt()
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::computeDt, level = " << m_level << endl;
  }

  Real newDt;
  newDt = m_dtNew;

  return newDt;
}


// Set the CFL number
void AMRLevelIdealMHD::CFL(Real a_cfl)
{
  m_cfl = a_cfl;
}

// Set the physical dimension of the longest side of the domain
void AMRLevelIdealMHD::domainLength(RealVect a_domainLength)
{
  m_domainLength = a_domainLength;
}



// Set the tag buffer size
void AMRLevelIdealMHD::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize = a_tagBufferSize;
}

// Set the density gradient output flag
void AMRLevelIdealMHD::output_density_gradient( bool a_output_density_gradient )
{	m_output_density_gradient=a_output_density_gradient;}

// Set the B gradient output flag
void AMRLevelIdealMHD::output_B_gradient( bool a_output_B_gradient )
{	m_output_B_gradient=a_output_B_gradient;}

bool AMRLevelIdealMHD::output_density_gradient() { return  m_output_density_gradient;}
bool AMRLevelIdealMHD::output_B_gradient(){ return m_output_B_gradient;}

// Create a load-balanced DisjointBoxLayout from a collection of Boxes
DisjointBoxLayout AMRLevelIdealMHD::loadBalance(const Vector<Box>& a_grids)
{
  // Load balance and create boxlayout
  Vector<int> procMap;
  //if (procID() == uniqueProc(SerialTask::compute))
  //{
  //  LoadBalance(procMap,a_grids);
  //}
  //broadcast(procMap,uniqueProc(SerialTask::compute));

  // appears to be faster for all procs to do the loadbalance (ndk)
  LoadBalance(procMap,a_grids);

  if (s_verbosity >= 4)
  {
    pout() << "AMRLevelIdealMHD::loadBalance: procesor map: " << endl;
    for (int igrid = 0; igrid < (int)a_grids.size(); ++igrid)
    {
      pout() << igrid << ": " << procMap[igrid] << "  " << endl;
    }
    pout() << endl;
  }

  DisjointBoxLayout dbl(a_grids,procMap,m_problem_domain);
  dbl.close();

  return dbl;
}


// Get the next coarser level
AMRLevelIdealMHD* AMRLevelIdealMHD::getCoarserLevel() const
{
  AMRLevelIdealMHD* amrMHDCoarserPtr = NULL;

  if (m_coarser_level_ptr != NULL)
  {
    amrMHDCoarserPtr = dynamic_cast<AMRLevelIdealMHD*>(m_coarser_level_ptr);

    if (amrMHDCoarserPtr == NULL)
    {
      MayDay::Error("AMRLevelIdealMHD::getCoarserLevel: dynamic cast failed");
    }
  }

  return amrMHDCoarserPtr;
}

// Get the next finer level
AMRLevelIdealMHD* AMRLevelIdealMHD::getFinerLevel() const
{
  AMRLevelIdealMHD* amrMHDFinerPtr = NULL;

  if (m_finer_level_ptr != NULL)
  {
    amrMHDFinerPtr = dynamic_cast<AMRLevelIdealMHD*>(m_finer_level_ptr);

    if (amrMHDFinerPtr == NULL)
    {
      MayDay::Error("AMRLevelIdealMHD::getFinerLevel: dynamic cast failed");
    }
  }

  return amrMHDFinerPtr;
}



// Interpolate source terms after they have been calculated at "0" level.
// Also fills the ghost cells if needed.
/*void AMRLevelIdealMHD::InterpolateSourceTerms()
{

  if (m_patchMHDAM->getSourceCalculator()==NULL) return;
  if (m_patchMHDAM->getSourceCalculator()->numSourceTerms()==0) return;

  AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();

  // Interpolate from coarser level
  if (m_hasCoarser)
  {
    int nRefCrse = m_coarser_level_ptr->refRatio();

    FineInterp fineInterpST;
    fineInterpST.define(m_grids,
                        m_SCData.nComp(),
                        nRefCrse,
                        m_problem_domain);

    fineInterpST.interpToFine(m_SCData,amrMHDCoarserPtr->m_SCData);
  }

  // Filling ghost cells
  if (m_SCData.ghostVect() != IntVect::Zero)
  {
    int numGhost = m_SCData.ghostVect()[0];
    if (m_hasCoarser)
    {
      const DisjointBoxLayout& levelDomain = m_SCData.disjointBoxLayout();

      PiecewiseLinearFillPatch pwl(levelDomain,
                               amrMHDCoarserPtr->m_SCData.disjointBoxLayout(),
                               m_SCData.nComp(),
                               amrMHDCoarserPtr->m_problem_domain,
                               amrMHDCoarserPtr->m_ref_ratio,
                               numGhost);

      pwl.fillInterp(m_SCData,
                 amrMHDCoarserPtr->m_SCData,
                 amrMHDCoarserPtr->m_SCData,
                 1.0,
                 0,
                 0,
                 m_SCData.nComp());
    }

    m_SCData.exchange(Interval(0,m_SCData.nComp()-1));
  }

//  char filename[60];filename[0]=0;
//  sprintf(filename,"mhdk_st%i.dat",m_level);
//  FILE* tfile = OpenTecplotFile(filename,"TITLE = \"source terms\"\n VARIABLES=\"X\" \"Y\" \"SMOM\" \"SENG\"");
//
//  const DisjointBoxLayout& levelDomain = m_SCData.disjointBoxLayout();
//  DataIterator dit = m_SCData.dataIterator();
//  for (dit.begin(); dit.ok(); ++dit)
//  {
//    const Box& b = levelDomain[dit()];
//    FArrayBox& SFab = m_SCData[dit()];
//    WriteFArrayBoxToTecplotFile(tfile, SFab, SFab.box(), Interval(0,1), m_dx);
//  }
//  CloseTecplotFile(tfile);

}*/


///
const LevelData<FArrayBox>& AMRLevelIdealMHD::getStateNew() const
{
  return m_UNew;
}

///
const LevelData<FArrayBox>& AMRLevelIdealMHD::getStateOld() const
{
  return m_UOld;
}

LevelData<FArrayBox>& AMRLevelIdealMHD::getSCData()
{
  return m_SCData;
}

const LevelData<FArrayBox>& AMRLevelIdealMHD::getdivB() const
{
  if (!m_divB.isDefined())
    MayDay::Error("AMRLevelIdealMHD::getdivB: m_divB is not defined");

  return m_divB;
}

/// Clear divB by applying projection scheme (solving Poisson equation)
  /**
   */
void AMRLevelIdealMHD::ApplyProjectionScheme(const LevelData<FArrayBox>& a_phi)
{
  Real dx = m_csh->dx(0,m_level);
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    const Box& b         = levelDomain[dit()];
    FArrayBox& UFab      = m_UNew[dit()];
    const FArrayBox& phi = a_phi[dit()];

    FORT_CORRECT_B(
         CHF_FRA(UFab),
         CHF_CONST_FRA1(phi,0),
         CHF_CONST_REAL(dx),
         CHF_BOX(b));                                	
  }
}

void AMRLevelIdealMHD::interpolateData(
  const RealVect & a_rv,
  const IntVect & a_iv,
  const FArrayBox & a_data,  
  Real* a_value)
{
  FArrayBox slopes[3],dx[3]; int dir,ncomp;
  
  Box b(a_iv,a_iv);
  
  ncomp = a_data.nComp();
  
  for (dir=0; dir<3; dir++)
  {
    slopes[dir].define(b, ncomp);        
    dx[dir].define(Box(IntVect::Zero,IntVect::Zero), 1);
    slopes[dir].setVal(-666.666);        
  }
             
  for (dir = 0; dir < SpaceDim; ++dir)
  {
    IntVect iv_off(IntVect::Zero);
    iv_off[dir]=1;
    
    Box dxBox(m_problem_domain.domainBox().smallEnd()*iv_off, 
              m_problem_domain.domainBox().bigEnd()  *iv_off);
    dxBox.grow(dir,1);
    dxBox = dxBox & m_problem_domain;
    
    dx[dir].define(dxBox,1);                     
    m_csh->dx(dx[dir],dx[dir].box(),dir,m_level);          
  }

  
  for (dir = 0; dir < SpaceDim; ++dir)
  {    
    FArrayBox& dir_slope = slopes[dir];        

    const Box bcenter = grow(m_problem_domain,-BASISV(dir)) & b;
    if (!bcenter.isEmpty())
    {              
      if (m_csh->constStep(dir) == true)
      {
        Real ldx = dx[dir].get(dx[dir].smallEnd(),0);
        FORT_INTERPCENTRALDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( bcenter ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( ldx )
                        );
      } else
      {          
        FORT_INTERPCENTRALDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( bcenter ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( dx[dir], 0 )
                        );                   
      }        
    }
    const Box blo = b & adjCellLo(grow(m_problem_domain,-BASISV(dir)),dir);
    if (!blo.isEmpty())
    {
      if (m_csh->constStep(dir) == true)
      {                    
        Real ldx = dx[dir].get(dx[dir].smallEnd(),0);
        FORT_INTERPHISIDEDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( blo ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( ldx )
                        );
      } else
      {          
        FORT_INTERPHISIDEDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( blo ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( dx[dir], 0 )
                        );                   
      }                
    }
    const Box bhi = b & adjCellHi(grow(m_problem_domain,-BASISV(dir)),dir);
    if (!bhi.isEmpty())
    {
      if (m_csh->constStep(dir) == true)
      {     
        Real ldx = dx[dir].get(dx[dir].smallEnd(),0);
        FORT_INTERPLOSIDEDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( bhi ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( ldx )
                        );
      } else
      {          
        FORT_INTERPLOSIDEDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_data ),
                        CHF_BOX ( bhi ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( dx[dir], 0 )
                        );                   
      }                          
    }    
  } 
  
  // to do limits, we need to have a box which includes
  // the neighbors of a given point (to check for the
  // local maximum...
  Box neighborBox(-1*IntVect::Unit,
                  IntVect::Unit);

  // GHM 7/12/01
  // interplimit iterates over box b_mod (was b), but cells within
  // 1 of the physical boundary never enter result (and this
  // wasted calculation may call upon uninitialized memory).
  // DFM 10/8/01
  // note that this turns off slope limiting for cells adjacent to the
  // boundary -- may want to revisit this in the future
  Box b_mod(b);
  b_mod.grow(1);
  b_mod = m_problem_domain & b_mod;
  b_mod.grow(-1);

  // create a box grown big enough to remove periodic BCs from domain
  Box domBox = grow(b, 2);
  domBox = m_problem_domain & domBox;
    
  FORT_INTERPLIMIT_V ( CHF_FRA ( slopes[0] ),
                   CHF_FRA ( slopes[1] ),
                   CHF_FRA ( slopes[2] ),
                   CHF_CONST_FRA ( a_data ),
                   CHF_CONST_FRA1 ( dx[0], 0 ),
                   CHF_CONST_FRA1 ( dx[1], 0 ),
                   CHF_CONST_FRA1 ( dx[2], 0 ),                   
                   CHF_BOX ( b_mod ),
                   CHF_BOX ( neighborBox ),
                   CHF_BOX (domBox)                   
                   );  
  
  RealVect cc;
  m_csh->getCellCenter(cc,a_iv,m_level);
  
  for (int icomp = 0; icomp < ncomp; ++icomp)
  {
    a_value[icomp] = a_data.get(a_iv,icomp);
  }
  
  for (dir=0; dir<SpaceDim; dir++)
  {   
    Real interp_coef = a_rv[dir] - cc[dir];
                                             
    for (int icomp = 0; icomp < ncomp; ++icomp)
      {
        Real slope = interp_coef * slopes[dir](a_iv,icomp);                          
        a_value[icomp] += slope;
      }
  }
                         
  
}  

/*
// Advance by one timestep
Real AMRLevelIdealMHD::advance()
{
  if ((s_verbosity >= 2) || ((s_verbosity <= -2)&&(procID() == 0)))
  {
    pout() << "  AMRLevelIdealMHD::advance level " << m_level << " to time " << m_time + m_dt << " (dt=" << m_dt << ")" << endl;
  }

  // Copy the new to the old

  DataIterator dit = m_UNew.dataIterator();
  for( ; dit.ok(); ++dit )
  {
	m_UOld[dit()].copy(m_UNew[dit()]);
  }

  Real newDt = 0.0;

  // Set up arguments to LevelMHDAM::step based on whether there are
  // coarser and finer levels

  // Undefined flux register in case we need it
  LevelFluxRegister dummyFR;

  // Undefined leveldata in case we need it
  const LevelData<FArrayBox> dummyData;

  // Set arguments to dummy values and then fix if real values are available
  LevelFluxRegister* coarserFR = &dummyFR;
  LevelFluxRegister* finerFR   = &dummyFR;

  const LevelData<FArrayBox>* coarserDataOld = &dummyData;
  const LevelData<FArrayBox>* coarserDataNew = &dummyData;
  const LevelData<FArrayBox>* coarserPhi     = &dummyData;

  Real tCoarserOld = 0.0;
  Real tCoarserNew = 0.0;

  if (s_verbosity >= 4)
  {
    pout() << "Finer and coarse levels level "  << endl;
  }

  // A coarser level exists
  if (m_hasCoarser)
  {
    AMRLevelIdealMHD* coarserPtr = getCoarserLevel();

    // Recall that my flux register goes between my level and the next
    // finer level
    coarserFR = &coarserPtr->m_fluxRegister;

    coarserDataOld = &coarserPtr->m_UOld;
    coarserDataNew = &coarserPtr->m_UNew;
    coarserPhi     = &coarserPtr->m_phi;

    tCoarserNew = coarserPtr->m_time;
    tCoarserOld = tCoarserNew - coarserPtr->m_dt;
  }

  // A finer level exists
  if( m_hasFiner == true )
  {
    // Recall that my flux register goes between my level and the next
    // finer level
    finerFR = &m_fluxRegister;

    if( m_CTused == true )
    {
      AMRLevelIdealMHD* pFiner = getFinerLevel();

      for( DataIterator ditF = pFiner->m_EAcc.dataIterator(); ditF.ok(); ++ditF )
      {
        pFiner->m_EAcc[ditF()].setVal( 0.0 );
      }
    }
  }
  

  if (s_verbosity >= 4)
  {
    pout() << "m_levelMHDAM.step "  << endl;
  }  
  
  // Advance the solve one timestep
  if (m_patchMHDAM->timeApproximation()==PatchMHDAM::TA_Hancock)
    newDt = m_levelMHDAM.stepHancock(m_UNew,
                                m_E,
                                *finerFR,
                                *coarserFR,
                                m_divB,
                                m_phi,
                                *coarserPhi,
                                m_solver,
                                m_SCData,
                                *coarserDataOld,
                                tCoarserOld,
                                *coarserDataNew,
                                tCoarserNew,
                                m_time,
                                m_dt);
  
   else      
    newDt = m_levelMHDAM.step(m_UNew,
                                m_E,
                                *finerFR,
                                *coarserFR,
                                m_divB,
                                m_SCData,
                                *coarserDataOld,
                                tCoarserOld,
                                *coarserDataNew,
                                tCoarserNew,
                                m_time,
                                m_dt);
  

  // Update the time and store the new timestep
  m_time += m_dt;
  Real returnDt = m_cfl * newDt;

  m_dtNew = returnDt;

  return returnDt;
}

*/

