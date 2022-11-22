#ifndef _AMR_LEVEL_IDEAL_MHD_H_
#define _AMR_LEVEL_IDEAL_MHD_H_

#include "CHOMBO_VERSION.H"

#include "FArrayBox.H"
#include "EdgeBox.H"
#include "LevelData.H"
#include "AMRLevel.H"
#include "LevelFluxRegister.H"
#include "Box.H"
#include "IntVectSet.H"
#include "Vector.H"
#include "DisjointBoxLayout.H"
#include <AMRMultiGrid.H> 
#include <BiCGStabSolver.H> 

#if CHOMBO_VERSION_MAJOR < 4
  #include "AMRPoissonOp.H"
#else
  #include "AMRPoissonOpFactory.H"
#endif  

#include "CSHandler.H"
#include "PatchMHDAM.H"
#include "CoarseAverageExt.H"
#include "CoarseAverageEdge.H"
#include "FineInterpBCL.H"
#include "FineInterpMHDAM.H"
#include "PiecewiseLinearFillPatchMHDAM.H"
#include "ProbeFilesInfo.H"


class CoordinateSystemHandler;


/// AMR MHD
/**
 */
class AMRLevelIdealMHD : public AMRLevel
{
  friend class AMRLevelIdealMHDFactory;
public:
  /// Constructor
  /**
   */
  AMRLevelIdealMHD();

  /// Destructor
  /**
   */
  virtual ~AMRLevelIdealMHD();

  /// This instance should never get called - historical
  /**
   */
  virtual void define(AMRLevel*  a_coarserLevelPtr,
                      const Box& a_problemDomain,
                      int        a_level,
                      int        a_refRatio);

  /// Define new AMR level
  /**
   */
  virtual void define(AMRLevel*            a_coarserLevelPtr,
                      const ProblemDomain& a_problemDomain,
                      int                  a_level,
                      int                  a_refRatio);

  /// Things to do before a timestep
  /**
   */
  void preTimeStep();

  /// One time step on current level
  /**
   */
  Real step(IntVect & a_minDtCell);
  
#ifdef CH_OMPCPP
  Real stepOMP(IntVect & a_minDtCell);
#endif
  
  /// One time step on current level using Hancock scheme with divB cleaning via Projection scheme
  /**
   */
  Real stepHancockProjection();
  
  /// Advance by one timestep
  /**
   */
  virtual Real advance();

  /// Things to do after a timestep
  /**
   */
  virtual void postTimeStep();

  /// Create tags for regridding
  /**
   */
  virtual void tagCells(IntVectSet& a_tags) ;

  /// Create tags at initialization
  /**
   */
  virtual void tagCellsInit(IntVectSet& a_tags) ;  

  // Performs any post-regridding operations which are necessary
  /**
   */
  virtual void postRegrid(int a_base_level);
  
#if CHOMBO_VERSION_MAJOR < 4         
  /// Initialize grids
  /**
   */
  virtual void initialGrid(const Vector<Box>& a_newGrids);
  
  /// Set up data on this level after regridding
  /**
   */
  virtual void regrid(const Vector<Box>& a_newGrids);
#else
  virtual void initialGrid(const DisjointBoxLayout & a_grid);
  virtual void regrid(const DisjointBoxLayout & a_grid);
#endif

  /// Initialize data
  /**
   */
  virtual void initialData();

  /// Things to do after initialization
  /**
   */
  virtual void postInitialize();

  /// Create a local PatchMHDAM factory using the argument as a factory
  /**
   */
  virtual void patchMHDAM(const PatchMHDAM* const a_patchMHDAM);

  /// ????
  /**
   */
  PatchMHDAM* getpatchMHDAM();
  

#ifdef CH_USE_HDF5
  /// Write checkpoint header
  /**
   */
  virtual void writeCheckpointHeader(HDF5Handle& a_handle) const;

  /// Write checkpoint data for this level
  /**
   */
  virtual void writeCheckpointLevel(HDF5Handle& a_handle) const;

  /// Read checkpoint header
  /**
   */
  virtual void readCheckpointHeader(HDF5Handle& a_handle);

  /// Read checkpoint data for this level
  /**
   */
  virtual void readCheckpointLevel(HDF5Handle& a_handle);

  /// Write plotfile header
  /**
   */
  virtual void writePlotHeader(HDF5Handle& a_handle) const;

  /// Write plotfile data for this level
  /**
   */
  virtual void writePlotLevel(HDF5Handle& a_handle) const;

#endif

  /// Write on level in parallel mode: one patch = one level
  /**
   */  
#ifdef CH_MPI  
  virtual  void writeLevelforTecPlot(MPI_File tecplot_file);
#else
  virtual  void writeLevelforTecPlot(FILE* tecplot_file);
#endif  

#ifdef CH_MPI  
   void writeLevel0forTecPlot(MPI_File tecplot_file);
#else
   void writeLevel0forTecPlot(FILE* tecplot_file);
#endif  

  /// ????
  /**
   */
#ifdef CH_MPI
  void writeSlice(MPI_File tecplot_file, char plane, Real sValue);
#else
  void writeSlice(FILE* tecplot_file, char plane, Real sValue);
#endif

  /// Extract data along coordinate line
  /**    
   */
  int writeOneLine(FILE * tecplot_file, Box a_box, int a_dirLine, LevelData<FArrayBox>* a_lineData);

  /// Output data to TecPlot format. One level = one zone
  /**
    a_nComp - number of components to output
   */
  void writeLevelforTecPlot_LZ(FILE* tecplot_file, int a_nComp);

  /// Write data along X axis for axisymmetric configuration (Honestly, I do not remember why it is needed)
  /**
   */
  void writeAxisXData(int a_step, Real a_time);
  
  /// Not used, scheduled to be deleted
  void write1DProbe(ProbeFilesInfo& a_pfi);
  
  /// Write curent probe into file
  int writeProbeData(const ProbeFilesInfo::ProbeInfo & a_probe);

  /// Load Pogorelov's fort.55 data into MS-FLUKSS
  /**
   */
  void loadfort55();

  /// ????
  /**
   */
  void writeLevel0SF(int a_step);// write level 0 in special format
  
  void MedvedevOutput(int a_step) const;
  
  void V1V2TSOutput(int a_step) const;
  
  void wrtieLevel0H5() const;

  /// Returns the dt computed earlier for this level
  /**
   */
  virtual Real computeDt();

  /// Compute dt using initial data
  /**
   */
  virtual Real computeInitialDt();
  
  /// Compute dt using for current data
  /**
   */
  Real computeNewDt(Real a_cfl);
  
  

  //virtual void computeUExtra();

  /// Set the CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Set the physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(RealVect a_domainLength);

  /// Set the tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// Set the density gradient output flag
  /**
   */
  void output_density_gradient( bool a_output_density_gradient );

  /// ????
  /**
   */
  void output_B_gradient( bool a_output_B_gradient );

  /// ????
  /**
   */
  bool output_density_gradient();

  /// ????
  /**
   */
  bool output_B_gradient( );
  

  /// ????
  /**
   */
  const LevelData<FArrayBox>& getStateNew() const;

  /// ????
  /**
   */
  const LevelData<FArrayBox>& getStateOld() const;

  /// ????
  /**
   */
  LevelData<FArrayBox>& getSCData();

  /// ????
  /**
   */
  const LevelData<FArrayBox>& getdivB() const;

  /// Get the next coarser level
  /**
   */
  AMRLevelIdealMHD* getCoarserLevel() const;

  /// Get the next finer level
  /**
   */
  AMRLevelIdealMHD* getFinerLevel() const;  

  /// Clear divB by applying projection scheme (solving Poisson equation)
  /**
   */
  void ApplyProjectionScheme(const LevelData<FArrayBox>& a_phi);
  
  /// Reinitialize level set
  /**
   */
  void reinitializeLevelSet();
  
  
  void interpolateData(
        const RealVect & a_rv,
        const IntVect & a_iv,
        const FArrayBox & a_data,  
        Real* a_value);
        
        
  /// Return number of points in grid      
  /**
   */
  virtual int	numPts();
  
  /// Return the grids at this level 
  /**
   */
  virtual const DisjointBoxLayout & boxes();
  
  /// Fill ghost cells for conservative variables
  /**
   */
  void fillGhostCellsCons(Interval & a_interval);

        
    
protected:
  // Create a load-balanced DisjointBoxLayout from a collection of Boxes
  DisjointBoxLayout loadBalance(const Vector<Box>& a_grids);

  // Setup menagerie of data structures
  void levelSetup();

////////////////////////////////////////////////////////////////// Data members
  // Number of ghost cells (in each direction) for m_UOld,m_UNew
  int m_numGhost;

  // Number of converved states
  int m_numStates;

  // Conserverd state, U, at old and new time
  LevelData<FArrayBox> m_UOld, m_UNew;

  // Electric field
  LevelData<EdgeBox>   m_E, m_EAcc;

  // Number of ghost cells (in each direction) for m_SCData
  int m_numSCGhost;

  // Data calculated by external source calculator
  // For kinetic problems this is source terms and neutrals parameters
  LevelData<FArrayBox> m_SCData;    

  // Inversion of volumes. Used for refluxing  
  LevelData<FArrayBox> m_invVolumes;

  // CFL number
  Real m_cfl;

  // Interpolation from fine to coarse level
  //FineInterpBCL m_fineInterp;
  FineInterpMHDAM m_fineInterp;
  //FineInterp m_fineInterp;
  
  // Interpolator for filling in ghost cells from the next coarser level
  PiecewiseLinearFillPatchMHDAM m_patcher;

  // Averaging from coarse to fine level
  CoarseAverageExt m_coarseAverage;
  //CoarseAverage m_coarseAverage;
  CoarseAverageEdge m_coarseAverageEdge;

  // Flux register
  LevelFluxRegister m_fluxRegister;

  // New time step
  Real m_dtNew;

  // Names of conserved states
  Vector<string> m_stateNames;

  // Physical dimension of the domain
  RealVect m_domainLength;
    
  // Patch integrator (factory and object)
  PatchMHDAM* m_patchMHDAMFactory;
  PatchMHDAM* m_patchMHDAM;

  // Tag buffer size
  int m_tagBufferSize;

  // Flag coarser and finer levels
  bool m_hasCoarser;
  bool m_hasFiner;

  // Output density gradient in tecplot format
  bool m_output_density_gradient;
  bool m_output_B_gradient;
  bool m_output_divB;
  
  // true if vectors are written in cartesian coordinates
  bool m_output_vecCS;
  
  bool m_CTused; // Flag whether constrained transport is used
  bool m_CPused; // Flag whether correction potential (Dedner) is used

  DisjointBoxLayout m_grids;
      
  // Coordinate system handler
  CoordinateSystemHandler* m_csh;
  
  // System of equations
  EquationSystem* m_eqSys;
  
  // Members for projection scheme  
  AMRMultiGrid<LevelData<FArrayBox> >     m_solver;    
  BiCGStabSolver<LevelData<FArrayBox> >   m_bottomSolver;    
  LevelData<FArrayBox>                    m_phi;  
  LevelData<FArrayBox>                    m_divB;
  
  // Containers for the boundary conditions along the Z-axis
  //
  // These two containers are used to collect data from m_UNew near Z-axis
  LevelData<FArrayBox> m_lstrips_0;
  LevelData<FArrayBox> m_lstripsPI;  
  //
  // These two containers are used to calculate BC near Z-axis
  LevelData<FArrayBox> m_lstrips;  // data near near z-axis
  LevelData<FArrayBox> m_ladjstr;  // corresponding ghost cells
  
  //PiecewiseLinearFillPatchMHDAM m_patcherZ0;
  //PiecewiseLinearFillPatchMHDAM m_patcherZPI;
    
private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelIdealMHD& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelIdealMHD(const AMRLevelIdealMHD& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
