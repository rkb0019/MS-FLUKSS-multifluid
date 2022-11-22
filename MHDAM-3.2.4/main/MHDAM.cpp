#include <iostream>
using std::ifstream;
using std::ios;
#include <assert.h>

#ifdef CH_OMPCPP  
  #include<omp.h>
#endif  

#include "CHOMBO_VERSION.H"


#include "ParmParse.H"
#include "CH_HDF5.H"
#include "CH_Attach.H"
#include "parstream.H"

#include "AMR_MHDAM.H"
#include "AMRLevel.H"
#include "AMRLevelIdealMHDFactory.H"
#include "AMRLevelIdealMHD.H"
#include "SourceCalculator.H"

#include "Reconstruction.H"
#include "LGintegrator.H"
#include "RefCriteria.H"
#include "MHDAMDefs.H"

// Coordinate systems
#include "CSHandler.H"
#include "CSHCartesian.H"
#include "CSHSpherical.H"

// Equation systems
#include "EquationSystem.H"
#include "EqSysMHDMF.H"
#include "EqSysEuler.H"

#include "TMBreechEtAl2008.H"
#include "PITwoEquations.H"

// Patches
#include "PatchMHDMF.H"
#include "PatchEuler.H"


// Riemann solvers
#include "Roe8WavesRS.H"
#include "RusanovRS.H"
#include "VanLeerRS.H"
#include "AUSMPlusRS.H"

// Problems
#include "RiemannProblem.H"
#include "OrszagProblem.H"
#include "RotorProblem.H"
#include "CloudProblem.H"
#include "KelvinProblem.H"
#include "BlastProblem.H"
#include "RampProblem.H"

#ifdef SWLISM
// Source calculators
#include "KineticSources1D.H"
#include "KineticSources2D.H"
#include "KineticSources3D.H"

#include "SWLISMProblem.H"
#include "SWLISMProblemPolar.H"

#include "HeliosphericProblem.H"
#include "HelioTILTProblem.H"
#include "HelioRealBCProblem.H"

#include "RiemannProblemMF.H"
#include "RiemannProblem_MHDK.H"

#endif

#include "ObliqueProblem.H"
#include "FastShock.H"
#include "FreeStreamSpherical.H"
#include "ShearFlowProblem.H"
#include "ObliqueShockTube.H"
#include "CurrentSheet.H"


#ifdef USE_ARRAYVIEW
#include "ArrayView.H"
#endif

extern "C" {
//#include <fpu_control.h>
#include <fenv.h>
}
/* IM: Invalid operation mask
 * DM: Denormalized operand mask
 * ZM: Zero-divide mask
 * OM: Overflow mask
 * UM: Underflow mask
 * PM: Precision (inexact result) mask
  ---(pm is kinda stupid)
*/



#ifdef NDEBUG
static void __attribute__ ((constructor)) trapfpe(void)
{
//  fpu_control_t cw =
//    _FPU_DEFAULT & ~(_FPU_MASK_IM | _FPU_MASK_ZM | _FPU_MASK_OM | _FPU_MASK_UM);
//  _FPU_SETCW(cw);
}
#else
#define _GNU_SOURCE 1
//static void __attribute__ ((constructor)) trapfpe(void)
static void trapfpe(void)
{
  
  //feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
  //feenableexcept(FE_ALL_EXCEPT);
  //fenv_t cw;
  //fegetenv(&cw);
  // Exception in all cases
  //feraiseexcept(FE_ALL_EXCEPT);
}
#endif

// Available sample problems
#define PROBLEM_RIEMANN       0
#define PROBLEM_ORSZAG        1
#define PROBLEM_ROTOR         2
#define PROBLEM_CLOUD         3
#define PROBLEM_SWLISM        4
#define PROBLEM_KELVIN        5
#define PROBLEM_BLAST         6
#define PROBLEM_RAMP          9
#define PROBLEM_RIEMANN_MHDK 11
#define PROBLEM_RIEMANN_MF   12
#define PROBLEM_HELIO        13
#define PROBLEM_OBLIQUE      14
#define PROBLEM_HELIOTILT    15
#define PROBLEM_FASTSHOCK    16
#define PROBLEM_FREESTREAMSPHERICAL 17
#define PROBLEM_SWLISMPOLAR  18
#define PROBLEM_HELIOREALBC  19
#define PROBLEM_SHEARFLOW    20
#define PROBLEM_CURRENTSHEET 21
#define PROBLEM_OBLIQUESHOCKTUBE 22


// Possible pressure relationships for the initial condition
#define PRESSURE_ISENTROPIC 0
#define PRESSURE_CONSTANT   1

// MHDAM is a function (as opposed to inline in main()) to get
// around MPI scoping problems
void MHDAM(const char* a_inFile);

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile);

// One more function for MPI
void dumpmemoryatexit();

int main(int a_argc, char* a_argv[])
{
#ifdef CH_MPI
  // Start MPI
  
#ifndef CH_OMPCPP  
  MPI_Init(&a_argc,&a_argv);
#else
  int provided,requested=MPI_THREAD_FUNNELED;
  MPI_Init_thread(&a_argc,&a_argv,MPI_THREAD_FUNNELED,&provided);
  if (procID() == 0)
  {
    pout() << "MPI_Init_thread, requested: " <<requested << ",provided: " << provided  << endl;
    #pragma omp parallel
    {
      #pragma omp master
      {
        pout() << "Number of threads: " << omp_get_num_threads()  << endl; 
      }
    }
    pout().flush();
  }
#endif  
  
#endif

  // Check for an input file
  char* inFile = NULL;
  if (a_argc > 1)
  {
    inFile = a_argv[1];
  }
  else
  {
    pout() << "Usage:  MHDAM...ex <inputfile>" << endl;
    pout() << "No input file specified" << endl;
    return -1;
  }

  
  // Parse the command line and the input file (if any)
  ParmParse pp(a_argc-2,a_argv+2,NULL,inFile);

#ifdef CH_MPI  
  /*MPI_Group new_chombo_group, cworld_group;
  int cur_size;
  MPI_Comm_size(MPI_COMM_WORLD, &cur_size);
  
  printf("qqqq %i\n",cur_size);
  
  ParmParse ppmhdam("mhdam");
  
  if (ppmhdam.contains("ncpus"))
  {
    int ncpus, ierr;
    ppmhdam.query("ncpus",ncpus);
    
    if ((ncpus>0) && (ncpus < cur_size))
    {      
      int * members = new int[ncpus];
      for (int i=0;i<ncpus;++i) members[i]=i;
      ierr = MPI_Comm_group(MPI_COMM_WORLD, &cworld_group);      
      CH_assert(ierr == MPI_SUCCESS);
      ierr = MPI_Group_incl(cworld_group, ncpus, members, &new_chombo_group);
      CH_assert(ierr == MPI_SUCCESS);
      ierr = MPI_Comm_create(MPI_COMM_WORLD, new_chombo_group, &Chombo_MPI::comm);
      CH_assert(ierr == MPI_SUCCESS);
      delete[] members;
      
      int nsize;
      
      MPI_Comm_size(Chombo_MPI::comm, &nsize);
      
      printf("new comm created %i\n",nsize);
    }
  } */ 
  
  setChomboMPIErrorHandler();
  
  if (procID() == 0)
  {
    pout() << "Chombo MPI tasks: " << numProc() << endl;
    pout().flush();
  }
  
#endif

    
#ifdef USE_ARRAYVIEW
  trapfpe();
#endif

#ifdef CH_USE_HDF5
  H5open();
#endif
  

  // Run MHDAM, i.e., do the computation
  MHDAM(inFile);

#ifdef CH_USE_MEMORY_TRACKING  
  dumpmemoryatexit();
#endif  

  
#ifdef CH_USE_HDF5  
//printf("closing c hdf5...\n");
  //H5close();
  //printf("closed\n");
#endif
  
#ifdef CH_MPI
  // Exit MPI
  MPI_Finalize();
#endif
}

void MHDAM(const char* a_inFile)
{
  
  //CH_TIME("MHDAM");
  
  //char* timerEnv = getenv("CH_TIMER");
  //if (timerEnv!=NULL) pout() << std::string(timerEnv) << endl;
  
  // Read inputs that are prefixed with "mhdam."
  ParmParse ppmhdam("mhdam");

  // This determines the amount of diagnositic output generated
  int verbosity = 3;
  ppmhdam.query("verbosity",verbosity);  

  // Determine the sample problem specified
  int problem = -1;
  std::string problemString;

  PhysProblem* pProblem = NULL;
  SourceCalculator * pSourCal  = NULL;

  if (ppmhdam.contains("problem"))
    {
      ppmhdam.query("problem",problemString);
      
      // Test problems
      if (problemString == "riemann")
        {
          problem = PROBLEM_RIEMANN;
          pProblem = new RiemannProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "orszag")
        {
          problem = PROBLEM_ORSZAG;
          pProblem = new OrszagProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "rotor")
        {
          problem = PROBLEM_ROTOR;
          pProblem = new RotorProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "cloud")
        {
          problem = PROBLEM_CLOUD;
          pProblem = new CloudProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "kelvin")
        {
          problem = PROBLEM_KELVIN;
          pProblem = new KelvinProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "blast")
        {
          problem = PROBLEM_BLAST;
          pProblem = new BlastProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "ramp")
        {
          problem = PROBLEM_RAMP;
          pProblem = new RampProblem;
          pSourCal = new SourceCalculator;
        }

#ifdef SWLISM        
      // 2D SW-LISM problems  
      // for historical reason we specify number of fluids in "problemString"
      if (problemString == "swlism")
        {
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMProblem(PhysProblem::PP_MHDPM);
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swlism2F")
        {
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMProblem(PhysProblem::PP_2FluidPM);
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swlism4F")
        {
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMProblem(PhysProblem::PP_4FluidPM);
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swlism5F")
        {
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMProblem(PhysProblem::PP_5FluidPM);
          pSourCal = new SourceCalculator;
        }  
      // Commented, use SWLISMProblem with m_subproblem == SW_TURB
      /*if (problemString == "swlismTurb")
        {
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMTurbProblem(PhysProblem::PP_MHDPM);
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swlism2FTurb")
        {
          problem = PROBLEM_SWLISM_2F;
          pProblem = new SWLISMTurbProblem(PhysProblem::PP_2FluidPM);
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swlism4FTurb")
        {
          problem = PROBLEM_SWLISM_4F;
          pProblem = new SWLISMTurbProblem(PhysProblem::PP_4FluidPM);
          pSourCal = new SourceCalculator;
        }  */
      if (problemString == "swlismPolar")
        {
          problem = PROBLEM_SWLISMPOLAR;
          pProblem = new SWLISMProblemPolar(PhysProblem::PP_MHDPM);
          pSourCal = new SourceCalculator;
        }        
      if (problemString == "swlismKinetic")
        {          
          problem = PROBLEM_SWLISM;
          pProblem = new SWLISMProblem(PhysProblem::PP_MHDPM);
          pSourCal = new KineticSources2D;
        }

      // 3D SW-LISM problems
      if (problemString == "heliosph")
        {
          problem = PROBLEM_HELIO;
          pProblem = new HeliosphericProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "swtilt")
        {
          problem  = PROBLEM_HELIOTILT;
          pProblem = new HelioTILTProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "helioRealBC")
        {
          problem = PROBLEM_HELIOREALBC;
          pProblem = new HelioRealBCProblem;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "helioKinetic")
        {
          problem = PROBLEM_HELIO;
          pProblem = new HeliosphericProblem;
          pSourCal = new KineticSources3D;
        }  
      
        
        
      // 1D kinetic problems  
      if (problemString == "riemann_MHDK")
        {
          problem = PROBLEM_RIEMANN_MHDK;
          pProblem = new RiemannProblem_MHDK;
          pSourCal = new KineticSources1D;
        }
        
      if ((problemString == "riemann_2F") || (problemString == "riemann_3F"))
        {
          problem = PROBLEM_RIEMANN_MF;
          pProblem = new RiemannProblemMF;
          pSourCal = new SourceCalculator;
        }      

#endif        
      // Miscellaneous problems           
      if (problemString == "oblique")
        {
          problem  = PROBLEM_OBLIQUE;
          pProblem = new ObliqueProblem;
          pSourCal = new SourceCalculator;
        }
      
      if (problemString == "fastshock")
        {
          problem  = PROBLEM_FASTSHOCK;
          pProblem = new FastShock;
          pSourCal = new SourceCalculator;
        }
      if (problemString == "FreeStreamSpherical")
        {
          problem = PROBLEM_FREESTREAMSPHERICAL;
          pProblem = new FreeStreamSpherical();
          pSourCal = new SourceCalculator;
        }

      if (problemString == "shearflow")
        {
          problem = PROBLEM_SHEARFLOW;
          pProblem = new ShearFlowProblem();
          pSourCal = new SourceCalculator;
        }

      if (problemString == "CurrentSheet")
        {
          problem = PROBLEM_CURRENTSHEET;
          pProblem = new CurrentSheet();
          pSourCal = new SourceCalculator;
        }

      if (problemString == "obliqueshocktube")
        {
          problem = PROBLEM_OBLIQUESHOCKTUBE;
          pProblem = new ObliqueShockTube();
          pSourCal = new SourceCalculator;
        }
                
        
      // The sample problem given isn't valid
      if (problem == -1)
        {
          pout() << "Invalid problem, \"" << problemString << "\", specified in input file" << endl << endl;
          return;
        }
    }
  else
    {
      // A sample problem must be specified
      pout() << "\"mhdam.problem\" not specified in input file" << endl << endl;
      return;
    }

  // Input parameters and initialise
  pProblem->input( ppmhdam, verbosity );

  int i;

  // Stop after this number of steps
  int nstop = 0;
  ppmhdam.query("max_step",nstop);

  // Stop when the simulation time get here
  Real stopTime = 0.0;
  ppmhdam.query("max_time",stopTime);
  
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  ppmhdam.queryarr("num_cells",numCells,0,SpaceDim);

  CH_assert(D_TERM(   (numCells[0] > 0),
                && (numCells[1] > 0),
                && (numCells[2] > 0)));
  CH_assert(D_TERM(   (numCells[0] % 2 == 0),
                && (numCells[1] % 2 == 0),
                && (numCells[2] % 2 == 0)));
                
  // Determine which spatial directions are periodic
  vector<int> isPeriodica(SpaceDim,0);
  bool isPeriodic[SpaceDim];

  ppmhdam.queryarr("is_periodic",isPeriodica,0,SpaceDim);
  // convert periodic from int->bool
  for (int dim=0; dim<SpaceDim; dim++)
  {
    isPeriodic[dim] = (isPeriodica[dim] == 1);
    if (isPeriodic[dim] && verbosity >= 2 && procID() == 0)
      pout() << "Using Periodic BCs in direction: " << dim << endl;
  }

  ProblemDomain probDomain (IntVect::Zero,
                            IntVect(D_DECL(numCells[0]-1,
                                           numCells[1]-1,
                                           numCells[2]-1)),
                            isPeriodic);

  // Set the physical size of the longest dimension of the domain
  RealVect domainLength = RealVect(D_DECL(1.0,1.0,1.0));
  Vector<Real> domainBox;domainBox.resize(2*SpaceDim);
  int num_comp = ppmhdam.countval("domain_length");
  if (num_comp == 1)
  {    
    Real longestDim;
    ppmhdam.query("domain_length",longestDim);
    Real dx = longestDim / probDomain.domainBox().longside();
    domainLength = RealVect(D_DECL(
      dx*probDomain.size(0),
      dx*probDomain.size(1),
      dx*probDomain.size(2)));
    for (i = 0; i < SpaceDim; ++i)
    {
      domainBox[i] = 0.0;      
      domainBox[SpaceDim+i] = domainBox[i] + domainLength[i];
    }
  } else
  {
    if (num_comp != 2*SpaceDim)
      MayDay::Error("Number of elements in \"domain_length\" should be 1 or SpaceDim");
      
    ppmhdam.queryarr("domain_length", domainBox, 0, 2*CH_SPACEDIM);          
    for (i = 0; i < SpaceDim; ++i)
      domainLength[i] = domainBox[SpaceDim+i]-domainBox[i];
  }

  // Maximum AMR level limit
  int maxLevel = 0;
  ppmhdam.query("max_level",maxLevel);
  int numReadLevels = Max(maxLevel,1);

  // Refinement ratios between levels
  std::vector<int> refRatios;
  // Note: this requires a refRatio to be defined for the finest level
  // (even though it will never be used)
  ppmhdam.queryarr("ref_ratio",refRatios,0,numReadLevels+1);
  

  // How far to extend refinement from cells newly tagged for refinement
  int tagBufferSize = 3;
  ppmhdam.query("tag_buffer_size",tagBufferSize);

  // Threshold that triggers refinement
  Real refineDensityThresh = -1.0;
  ppmhdam.query ("refine_density_thresh",refineDensityThresh);
  
  int density_maxlevel = maxLevel;
  ppmhdam.query ("refine_density_maxlevel", density_maxlevel);

  Real refineBThresh = -1.0;
  ppmhdam.query ("refine_B_thresh",refineBThresh);
  
  int B_maxlevel = maxLevel;
  ppmhdam.query ("refine_B_maxlevel", B_maxlevel);

  // Minimum dimension of a grid
  int blockFactor = 1;
  ppmhdam.query("block_factor",blockFactor);

  // Maximum dimension of a grid
  int maxGridSize = 32;
  ppmhdam.query("max_grid_size",maxGridSize);

  // CFL multiplier
  Real cfl = 0.8;
  ppmhdam.query("cfl",cfl);

  // Initial CFL multiplier
  Real initialCFL = 0.1;
  ppmhdam.query("initial_cfl",initialCFL);

  // Determine if a fixed or variable time step will be used
  Real fixedDt = -1;
  ppmhdam.query("fixed_dt",fixedDt);

  //int PoissonInterval = 0;
  //ppmhdam.query("poisson_interval",PoissonInterval);  
  //bool bSolvePoisson = (PoissonInterval > 0 ? true : false);


  limiter1D::eLimiters eLimRho = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimU   = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimV   = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimW   = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimP   = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimBX  = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimBY  = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimBZ  = limiter1D::eVanLeer;
  limiter1D::eLimiters eLimLS  = limiter1D::eMonotonizedCD;

  int iLimiter = 1;
  ppmhdam.query( "limiter_RHO", iLimiter ); eLimRho = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_VX",  iLimiter ); eLimU   = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_VY",  iLimiter ); eLimV   = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_VZ",  iLimiter ); eLimW   = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_P",   iLimiter ); eLimP   = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_BX",  iLimiter ); eLimBX  = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_BY",  iLimiter ); eLimBY  = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_BZ",  iLimiter ); eLimBZ  = (limiter1D::eLimiters)iLimiter;
  ppmhdam.query( "limiter_LS",  iLimiter ); eLimLS  = (limiter1D::eLimiters)iLimiter;

  // Print the parameters
  if (( verbosity >= 3 ) && (procID() == 0))
    {
      pout() << "maximum step = " << nstop << endl;
      pout() << "maximum time = " << stopTime << endl;

      pout() << "number of cells = " << D_TERM(numCells[0] << "  " <<,
                                               numCells[1] << "  " <<,
                                               numCells[2] << ) endl;

      pout() << "maximum level = " << maxLevel << endl;

      pout() << "refinement ratio = ";
    for (size_t i = 0; i < refRatios.size(); ++i) pout() << refRatios[i] << " ";
      pout() << endl;


      pout() << "refinement threshold for density = " << refineDensityThresh << endl;
      pout() << "refinement threshold for B vector = " << refineBThresh << endl;

      pout() << "blocking factor = " << blockFactor << endl;
      pout() << "max grid size = " << maxGridSize << endl;

      
      pout() << "CFL = " << cfl << endl;
      pout() << "initial CFL = " << initialCFL << endl;
      if (fixedDt > 0)
        {
          pout() << "fixed dt = " << fixedDt << endl;
        }

      pout() << "density    limiter = " << eLimRho << endl;
      pout() << "X-velocity limiter = " << eLimU   << endl;
      pout() << "Y-velocity limiter = " << eLimV   << endl;
      pout() << "Z-velocity limiter = " << eLimW   << endl;
      pout() << "pressure   limiter = " << eLimP   << endl;
      pout() << "BX         limiter = " << eLimBX  << endl;
      pout() << "BY         limiter = " << eLimBY  << endl;
      pout() << "BZ         limiter = " << eLimBZ  << endl;
    }

  int nTrackingSurfaces = 0;
  ppmhdam.query("num_tracking_surfaces",nTrackingSurfaces);  
  nTrackingSurfaces = MAX(0,nTrackingSurfaces);
  
  if (nTrackingSurfaces == 0)
  {
    int nRegionTracer = -1;
    ppmhdam.query("region_tracer",nRegionTracer);  
    if (nRegionTracer > 0)  nTrackingSurfaces++; 
  }
  
  int ilsOnly = false;
  ppmhdam.query("lsonly",ilsOnly); 

  TurbulenceModel * pTurbMod = NULL;
  int nModel = 0;
  ppmhdam.query( "turbModel", nModel );
  if( nModel == 1 )
  {
    pTurbMod = new TMBreechEtAl2008();
    pTurbMod->input( ppmhdam, verbosity );
  }

  if (( verbosity >= 3 ) && (procID() == 0))
  {
    pout() << "Turbulence model = " << nModel  << endl;
  }

  nModel = 0;
  PickupIons * pPIModel = NULL;
  ppmhdam.query( "pickupIonsModel", nModel );
  if( nModel == 1 )
  {
    pPIModel = new PITwoEquations();
    pPIModel->input( ppmhdam, verbosity );
  }

  if (( verbosity >= 3 ) && (procID() == 0))
  {
    pout() << "Pickup Ions model = " << nModel  << endl;
  }

  AMRLevelIdealMHDFactory amrMHDFact;

  RefCountedPtr<PatchMHDAM> patchPtr(NULL);
  
  int iDivCleaning = 0;
  if( ppmhdam.contains( "divB_cleaning" ) )
  {
    ppmhdam.query( "divB_cleaning", iDivCleaning );
  }
  
  EquationSystem * eqSys;
  
  int iDedner = 0; Real factorCp = 1.0,factorCh = 1.0;
  if (iDivCleaning == 5)
  {
    iDedner = 1;    
    ppmhdam.query("factorCh",factorCh); 
    ppmhdam.query("factorCp",factorCp); 
  }
  
  switch( pProblem->physicalModel() ) {
  case PhysProblem::PP_EulerPM :
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchEuler()));
    eqSys    = new EqSysEuler(nTrackingSurfaces);
    break;
  case PhysProblem::PP_MHDPM : {
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchMHDMF(1)));
    EqSysMHDMF * eqSysMF    = new EqSysMHDMF(iDedner,1,nTrackingSurfaces);
    eqSysMF->setCorrectionPotentialParams(iDedner,factorCh,factorCp);
    eqSys = eqSysMF;
    break; }
  case PhysProblem::PP_2FluidPM : {
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchMHDMF(2)));
    EqSysMHDMF * eqSysMF    = new EqSysMHDMF(iDedner,2,nTrackingSurfaces);
    eqSysMF->setCorrectionPotentialParams(iDedner,factorCh,factorCp);

    VanLeerRS * pVLeer = new VanLeerRS();

    patchPtr->setRiemannSolverGD( pVLeer );
    pProblem->setRiemannSolverGD( pVLeer );
    
    eqSys = eqSysMF;
    break;
    }
  case PhysProblem::PP_3FluidPM : {
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchMHDMF(3)));
    EqSysMHDMF * eqSysMF    = new EqSysMHDMF(iDedner,3,nTrackingSurfaces);
    eqSysMF->setCorrectionPotentialParams(iDedner,factorCh,factorCp);

    VanLeerRS * pVLeer = new VanLeerRS();

    patchPtr->setRiemannSolverGD( pVLeer );
    pProblem->setRiemannSolverGD( pVLeer );
    
    eqSys = eqSysMF;
    break;
    }
  case PhysProblem::PP_4FluidPM : {
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchMHDMF(4)));
    EqSysMHDMF * eqSysMF    = new EqSysMHDMF(iDedner,4,nTrackingSurfaces);
    eqSysMF->setCorrectionPotentialParams(iDedner,factorCh,factorCp);

    VanLeerRS * pVLeer = new VanLeerRS();

    patchPtr->setRiemannSolverGD( pVLeer );
    pProblem->setRiemannSolverGD( pVLeer );
    
    eqSys = eqSysMF;
    
    break;
    }
  case PhysProblem::PP_5FluidPM : {
    patchPtr = RefCountedPtr<PatchMHDAM>(static_cast<PatchMHDAM*>(new PatchMHDMF(5)));
    EqSysMHDMF * eqSysMF    = new EqSysMHDMF(iDedner,5,nTrackingSurfaces);
    eqSysMF->setCorrectionPotentialParams(iDedner,factorCh,factorCp);

    VanLeerRS * pVLeer = new VanLeerRS();

    patchPtr->setRiemannSolverGD( pVLeer );
    pProblem->setRiemannSolverGD( pVLeer );
    
    eqSys = eqSysMF;
    
    break;
    }
  }

  if( pTurbMod != NULL )
  {
    eqSys->setTurbulenceModel( pTurbMod );
  }

  if( pPIModel != NULL )
  {
    eqSys->setPickupIons( pPIModel );
  }

  pProblem->setEquationSystem( eqSys );
  pProblem->setSourceCalculator( pSourCal );

  pSourCal->input( ppmhdam, verbosity );
  patchPtr->input( ppmhdam, verbosity );

  int iCoordSys = 0;
  ppmhdam.query( "coord_sys", iCoordSys );
  CoordinateSystemHandler::eCoordinateSystem eCS = (CoordinateSystemHandler::eCoordinateSystem)(iCoordSys);
  CoordinateSystemHandler* CSH = NULL;  
  if ( (eCS == CoordinateSystemHandler::CS_Cartesian) || 
       (eCS == CoordinateSystemHandler::CS_Axisymmetric) || 
       (eCS == CoordinateSystemHandler::CS_Cylindrical))
  {
    CSH = new CSHCartesian;
    CSH->input( ppmhdam, verbosity );
    static_cast<CSHCartesian*>(CSH)->define(eCS, maxLevel, refRatios, probDomain.domainBox(), domainLength, *eqSys);    
  }     
#ifdef _CSHSPHERICAL_H_
  if ( (eCS == CoordinateSystemHandler::CS_Polar) || 
       (eCS == CoordinateSystemHandler::CS_PolarAxisym) || 
       (eCS == CoordinateSystemHandler::CS_Spherical))
  {
    if  ( ((eCS == CoordinateSystemHandler::CS_Polar) || 
           (eCS == CoordinateSystemHandler::CS_PolarAxisym)) && (SpaceDim != 2) )
      MayDay::Error("error in mhdam.coord_sys (polar coordinates requires SpaceDim == 2)");
    
    if  ( (eCS == CoordinateSystemHandler::CS_Spherical) && (SpaceDim != 3) )
      MayDay::Error("error in mhdam.coord_sys (spherical coordinates requires SpaceDim == 3)");
    
    if (domainBox.size() != 2*SpaceDim)
      MayDay::Error("Number of elements in \"domain_length\" should be equal to SpaceDim");
                  
    CSH = new CSHSpherical;
    CSH->input( ppmhdam, verbosity );
    static_cast<CSHSpherical*>(CSH)->define(eCS, maxLevel, refRatios, probDomain.domainBox(), domainLength, domainBox, *eqSys);    
  }     
#endif

  pProblem->setCoordinateSystem( CSH );

  if (!ppmhdam.contains("restart_file"))
  {
    pProblem->defineMesh(probDomain,domainBox);
    
    if (ppmhdam.contains("geom_file"))
    {
      std::string geomFile;
      ppmhdam.query("geom_file",geomFile);
      vector<int>  geomDirs(SpaceDim,0);
      bool        bgeomDirs[SpaceDim];
      ppmhdam.queryarr("geom_dirs",geomDirs,0,SpaceDim);
      for (int dim=0; dim<SpaceDim; dim++) bgeomDirs[dim]=(geomDirs[dim]==1);

#ifdef CH_USE_HDF5
      hid_t h5f = H5Fopen(geomFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
      CSH->readGeomInfo(h5f,bgeomDirs);
      H5Fclose(h5f);
#endif
    }
  }


  // Set up the patch integrator for ideal MHD

  int iFlux = 0;
  ppmhdam.query( "flux", iFlux );

  Real entropy_fix_coeff = 0.3;
  ppmhdam.query( "entropy_fix_coeff", entropy_fix_coeff );
    
  RiemannSolver * pRieman = RiemannSolver::make( (iDivCleaning == 5 ? iFlux+100 : iFlux),
                                                 pProblem->physicalModel(),
                                                 entropy_fix_coeff );
                                                 
  patchPtr->setRiemannSolver( pRieman );
  pProblem->setRiemannSolver( pRieman );  


  Reconstruction * pRec = new Reconstruction( eqSys->numSlopes() );
  pRec->setParametersForAll();
  pRec->setLimiterForAll( limiter1D::eMinmod );

  int iCharRecon  = 0;
  if( ppmhdam.contains( "char_reconstruction" ) )
  {
    ppmhdam.query( "char_reconstruction", iCharRecon );
  }

  if( iCharRecon > 0 )
  {
    pRec->setCharaterictic( iCharRecon );

    int iCheckTVD  = 0;
    if( ppmhdam.contains( "check_tvd" ) )
    {
      ppmhdam.query( "check_tvd", iCheckTVD );

      switch( iCheckTVD )
      {
      case 1 :
        {
          pRec->slope(WRHO )->setTVDCheck( true );
          pRec->slope(WVELX)->setTVDCheck( true );
          pRec->slope(WVELY)->setTVDCheck( true );
          pRec->slope(WVELZ)->setTVDCheck( true );
          pRec->slope(WPRES)->setTVDCheck( true );

          if( pProblem->physicalModel() != PhysProblem::PP_EulerPM )
          {
            pRec->slope(WBX)->setTVDCheck( true );
            pRec->slope(WBY)->setTVDCheck( true );
            pRec->slope(WBZ)->setTVDCheck( true );
          }
        }
      case 2 :
        {
          pRec->slope(WRHO )->setTVDCheck( true );
          pRec->slope(WPRES)->setTVDCheck( true );
        }
      }
    }
  }

  pRec->slope(WRHO )->setLimiter( eLimRho );
  pRec->slope(WVELX)->setLimiter( eLimU );
  pRec->slope(WVELY)->setLimiter( eLimV );
  pRec->slope(WVELZ)->setLimiter( eLimW );
  pRec->slope(WPRES)->setLimiter( eLimP );

  if( ppmhdam.contains( "check_positivity" ) )
  {
    int iCheckPos  = 0;
    ppmhdam.query( "check_positivity", iCheckPos );
    if( iCheckPos != 0 )
    {
      pRec->slope(WRHO )->checkPositivity();
      pRec->slope(WPRES)->checkPositivity();
    }
  }

  if( pProblem->physicalModel() != PhysProblem::PP_EulerPM )
  {
    pRec->slope(WBX)->setLimiter( eLimBX );
    pRec->slope(WBY)->setLimiter( eLimBY );
    pRec->slope(WBZ)->setLimiter( eLimBZ );
  }
  
  if( pPIModel != NULL )
  {
    int iRhoPIW   = pPIModel->primInterval().begin();
    int iPressPIW = iRhoPIW + 1;
    pRec->slope(iRhoPIW)->setLimiter( eLimRho );
    pRec->slope(iPressPIW)->setLimiter( eLimP );
    //consider adding check positivity
  }
  
  for (i=0;i<eqSys->numTrackingSurfaces();++i)
    pRec->slope(eqSys->lsIndexPrim(i))->setLimiter(eLimLS);


  patchPtr->setReconstruction( pRec );

  patchPtr->setPhysProblem(pProblem);
  
  patchPtr->setLSonlyFlag(ilsOnly == 1);


  int iMethodT = 0;
  if( ppmhdam.contains( "time_approximation" ) )
  {
    ppmhdam.query( "time_approximation", iMethodT );
  }
  PatchMHDAM::eTimeApproximation MethodT = (PatchMHDAM::eTimeApproximation)iMethodT;
  patchPtr->setTimeApproximation( MethodT );

  int iDivB  = (pProblem->physicalModel() == PhysProblem::PP_EulerPM) ? 0 : 1;
  if( ppmhdam.contains( "divB_method" ) )
  {
    ppmhdam.query( "divB_method", iDivB );
  }
  patchPtr->setDivBMethod( iDivB );
    
  patchPtr->setDivergenceCleaning( (eDivergenceCleaning)iDivCleaning );
  if ((iDivCleaning == 4) && (iMethodT != 1)) 
    MayDay::Error("divB_cleaning = 4 parameter must be used with time_approximation = 1");
  
  int i8waveUse = 1;
  if( ppmhdam.contains( "8wave_use" ) )
  {
    ppmhdam.query( "8wave_use", i8waveUse );
  }
  patchPtr->set8waveFlag( i8waveUse );
  
  if (iDivCleaning == 5) 
  {
    if (i8waveUse == 1) MayDay::Error("divB_cleaning = 5 can not be used with 8wave_use = 1");
    if (! ((iFlux == 1) || (iFlux == 3))) MayDay::Error("divB_cleaning = 5 must be used with iFlux = 1 or 3");
  }
     

  // Set parameters of refinement. It should be revised later  
  if (refineDensityThresh > 0.0)
  {
    patchPtr->RefineParams[REF_RHO].m_PerformRef = true;
    patchPtr->RefineParams[REF_RHO].m_Threshold  = refineDensityThresh;
    patchPtr->RefineParams[REF_RHO].m_maxlevel   = density_maxlevel;
  }
  if (refineBThresh > 0.0)
  {
    patchPtr->RefineParams[REF_B].m_PerformRef = true;
    patchPtr->RefineParams[REF_B].m_Threshold  = refineBThresh;
    patchPtr->RefineParams[REF_B].m_maxlevel   = B_maxlevel;
  }
  
  
    
  // Set up the AMRLevel... factory
  amrMHDFact.CFL(cfl);
  amrMHDFact.domainLength(domainLength);  
  amrMHDFact.tagBufferSize(tagBufferSize);
  amrMHDFact.patchMHDAM(patchPtr);
  amrMHDFact.verbosity(verbosity);
  amrMHDFact.initialDtMultiplier(initialCFL);
  amrMHDFact.input(ppmhdam, verbosity);

  AMR_MHDAM amr;

  // Set up the AMR object
  amr.define(maxLevel,refRatios,probDomain,&amrMHDFact,pSourCal);

  if (fixedDt > 0)
    {
      amr.fixedDt(fixedDt);
    }

  // Set grid generation parameters
  
  #if CHOMBO_VERSION_MAJOR < 4
  amr.maxGridSize(maxGridSize);
  amr.blockFactor(blockFactor);
  #endif

  amr.input( ppmhdam, verbosity );


  // Set up input files
  if (!ppmhdam.contains("restart_file"))
    {
      if (!ppmhdam.contains("fixed_hierarchy"))
        {
          // initialize from scratch for AMR run
          // initialize hierarchy of levels
          
          #if CHOMBO_VERSION_MAJOR < 4
          amr.setupForNewAMRRun();          
          # else
          amr.setupForNewAMRRun(maxGridSize);          
          #endif
        }
      else
        {
                  
        /* To be rewritten for chombo 4
                  
          std::string gridFile;
          ppmhdam.query("fixed_hierarchy",gridFile);

          // initialize from a list of grids in "gridFile"
          Vector<Vector<Box> > amrGrids(maxLevel+1);
          setupFixedGrids(amrGrids,
                          probDomain,
                          maxLevel,
                          maxGridSize,
                          blockFactor,
                          verbosity,
                          gridFile);
          amr.setupForFixedHierarchyRun(amrGrids,1);*/
        } 
                     
    }
  else
    {
      std::string restartFile;
      ppmhdam.query("restart_file",restartFile);

      if ( verbosity >= 2 )
       pout() << "restart_file = " << restartFile << endl;


#ifdef CH_USE_HDF5
      HDF5Handle handle(restartFile,HDF5Handle::OPEN_RDONLY);
      // read from checkpoint file
      CSH->readGeomInfo(handle.fileID());
      amr.setupForRestart(handle,maxGridSize);         
      handle.close();
#else
      MayDay::Error("MHDAM restart only defined with hdf5");
#endif
    }
  
  pSourCal->initExternalSC(amr.getAMRLevels());
  
  
  amr.setInputFile(a_inFile);  
  amr.setup1DProbeFiles();

  // Run the computation
  amr.run(stopTime,nstop);

  // Output the last plot file and statistics
  amr.conclude();
      
  delete CSH;
  delete eqSys;
  delete pSourCal;

  if (procID() == 0)
    CH_TIMER_REPORT();
}

// setupFixedGrids allows fixed grids to be read in and used in this AMR
// computation example
void setupFixedGrids(Vector<Vector<Box> >& a_amrGrids,
                     const ProblemDomain&  a_domain,
                     int                   a_maxLevel,
                     int                   a_maxGridSize,
                     int                   a_blockFactor,
                     int                   a_verbosity,
                     std::string           a_gridFile)
{
#if CHOMBO_VERSION_MAJOR < 4
  // Run this task on one processor
  if (procID() == uniqueProc(SerialTask::compute))
  {
    a_amrGrids.push_back(Vector<Box>(1,a_domain.domainBox()));

    // Read in predefined grids
    ifstream is(a_gridFile.c_str(), ios::in);

    if (is.fail())
    {
      MayDay::Error("Cannot open grids file");
    }

    // Format of file:
    //   number of levels, then for each level (starting with level 1):
    //   number of grids on level, list of boxes

    int inNumLevels;
    is >> inNumLevels;

    CH_assert (inNumLevels <= a_maxLevel+1);

    if (a_verbosity >= 3)
    {
      pout() << "numLevels = " << inNumLevels << endl;
    }

    while (is.get() != '\n');

    a_amrGrids.resize(inNumLevels);

    // Check to see if coarsest level needs to be broken up
    domainSplit(a_domain,a_amrGrids[0],a_maxGridSize,a_blockFactor);

    if (a_verbosity >= 3)
    {
      pout() << "level 0: ";
      for (size_t n = 0; n < a_amrGrids[0].size(); n++)
      {
        pout() << a_amrGrids[0][0] << endl;
      }
    }

    // Now loop over levels, starting with level 1
    int ngrid;
    for (int lev = 1; lev < inNumLevels; lev++)
    {
      is >> ngrid;

      if (a_verbosity >= 3)
      {
        pout() << "level " << lev << " numGrids = " << ngrid << endl;
        pout() << "Grids: ";
      }

      while (is.get() != '\n');

      a_amrGrids[lev].resize(ngrid);

      for (int i = 0; i < ngrid; i++)
      {
        Box bx;
        is >> bx;

        while (is.get() != '\n');

        // Quick check on box size
        Box bxRef(bx);

        if (bxRef.longside() > a_maxGridSize)
        {
          pout() << "Grid " << bx << " too large" << endl;
          MayDay::Error();
        }

        if (a_verbosity >= 3)
        {
          pout() << bx << endl;
        }

        a_amrGrids[lev][i] = bx;
      } // End loop over boxes on this level
    } // End loop over levels
  }

  // Broadcast results to all the processors
  broadcast(a_amrGrids,uniqueProc(SerialTask::compute));
  #endif
}
