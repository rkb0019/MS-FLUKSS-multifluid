#ifndef  HELIOSPHERICPROBLEM_H_
#define  HELIOSPHERICPROBLEM_H_

#if CH_LANG_CC
// Boundary condition ans sub problem constants
#endif

#if CH_LANG_CC
// Subproblems
#endif
#define HPBC_DEFAULT         0
#define HPBC_FIRST           1
#define HPBC_SUESS_MIN_FIRST 2
#define HPBC_SUESS_MIN_THIRD 3
#define HPBC_THIRD           4
#define HPBC_CIR             5
#define HPBC_PIZZO           6
#define HPBC_LONGTAIL        7
#define HPBC_HPINSTAB        8
#define HPBC_KINETIC         9
#define HPBC_SOLARCYCLE     10

#if CH_LANG_CC
// HelioRealBCProblem constants

// HPBC_WSOPROBLEM - HCS shape in the inner heliosheath
// HPBC_CIRPROBLEM - CIR data determined on sphere
#endif

#define HPBC_GMIRPROBLEM  100
#define HPBC_SUESSPROBLEM 101
#define HPBC_OMNIPROBLEM  102
#define HPBC_WSOPROBLEM   103
#define HPBC_CIRPROBLEM   104
#define HPBC_HDF5INPUT    105

#if CH_LANG_CC
// HelioProblem adaptation constants
#endif

#define HPAD_CARTESIAN      0
#define HPAD_LONGTAIL       1
#define HPAD_HPINSTAB       2
#define HPAD_TAILDETAIL     3

#if CH_LANG_CC

#include "MultiFluidProblem.H"

/// Base class for 3D heliospheric problem
/**
   Parameters:
 */
class HeliosphericProblem : public MultiFluidProblem
{
public:
  /// Null constructor
  /**
   */
  HeliosphericProblem();

  /// Destructor
  /**
   */
  ~HeliosphericProblem() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( );

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new IBCs made with
     new_PhysProblem() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();
  
  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this object.
   */
  virtual void define(const ProblemDomain& a_domain,                      
                      const int            a_level);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  virtual PhysProblem *new_PhysProblem();
  
  /// Copy method 
  /**
     Copy all data from a_PP to this instance.
   */
  virtual void copy_PhysProblem(const PhysProblem* a_PP);

  /// Set the Equation System flag
  /**
   */
  virtual void setEquationSystem( EquationSystem* a_eqSys );

  /// Fill ghost cells
  /**
   */
  virtual void fillGhostCells(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time ); //added

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(       FArrayBox&      a_F,
                             FArrayBox&      a_Bn,
                       const FArrayBox&      a_W,
                       const FArrayBox&      a_Wextrap,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const Real&           a_time);

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);
  
  virtual void initialize(LevelData<FArrayBox>& a_U,
                                     Interval & a_comp);


  /// Problem specific postprocessing
  /**
   */
  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const Real      & a_time,                               
                               const Box       & a_box);
                               
                               
  /// Things to do before advancing one level by one time step.
  /**
   */
  virtual void preTimeStep(LevelData<FArrayBox>&  a_U, Real a_time);

                               
  /// Things to do after advancing one level by one time step.
  /**
   */
  virtual void postTimeStep(LevelData<FArrayBox>&       a_U);

  /// Creates tagged cells for dynamic mesh refinement
  /**
    Problem specific cells tagging
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags);

  void tagCellsCartesian(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags);
                               
                               
  void  tagCellsLongTail(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags);
							   
  void  tagCellsTailDetail(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags);
                               
  void  tagCellsHPInstab(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags);

  // Tags HCS region (no level set)                              
  void tagCellsTilt(const FArrayBox&  a_U,
                    const Box&        a_box,
                          IntVectSet& a_tags);
 
  // Tags HCS region using level set data
  void tagCellsTiltLS(const FArrayBox&  a_U,
                      const Box&        a_box,
                      IntVectSet& a_tags);
  
  // Tags HCS region between V1 and V2
  void tagCellsTiltLSV1V2(const FArrayBox&  a_U,
                          const Box&        a_box,
                          IntVectSet& a_tags);

// generic tags for spherical grid
  void tagCellsSpherical(const FArrayBox&  a_U,
                          const Box&        a_box,
                          IntVectSet& a_tags);

// Tags for refined IHS, OHS, and HP (Federico)  
  void tagCellsTailHP_Federico(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags);

  /// Check geometrical/problem limitations for grid adaptation
  /**
   */  
  virtual void lockedCellsRegrid( BaseFab<int> & a_flag,
                            const FArrayBox&  a_U,
                            const Box&     a_box);
                            
                            
  /// Compute problem specific dt and returns a cell with the minimum dt. 
  /**
   */                               
  virtual Real computeDt( const FArrayBox& a_U,
                          const FArrayBox& a_dt,
                          const Box&     a_box,
                          IntVect&       a_minDtCell);

  /// Return boundary condition flags for all boundaries
  /**
   */
  virtual void getBCFlags( eBoundaryConditions leftBC,
                           eBoundaryConditions rightBC,
                           eBoundaryConditions bottomBC,
                           eBoundaryConditions topBC,
                           eBoundaryConditions frontBC,
                           eBoundaryConditions behindBC );

  /// Return the boundary condition flag for given boundary
  /**
   */
  virtual eBoundaryConditions getBCFlags( int a_dir, Side::LoHiSide a_sd );

  /// Define regions for charge exchange
  /**
   */
  virtual void defineRegions( const FArrayBox    & a_W,
                                    FArrayBox    & a_S,
                                    BaseFab<int> & a_R,
                              const Box          & a_box );
                              
  /// Number additional variables for writing to plot file
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPlotVars();
  
  /// Names of the additional variables for writing to plot file
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual Vector<std::string> plotNames();
  
  /// Calculates variables for plotting using primitive variables
  /**
  */
  virtual void calcPlotVars(FArrayBox&      a_vars,
                           int              a_comp,
                           const FArrayBox& a_W,
                           const Box&       a_box);
                              
  /// Scale primitive variables for plotting
  virtual void primForPlot(      FArrayBox& a_W,
                           const Box&       a_box);
                           
  /// Modifies curvilinear mesh (applies to curvilinear meshes only)
  /**
   */                               
  virtual void defineMesh(const ProblemDomain & a_prob_domain,
                          const Vector<Real>  & a_domainBox);
                          
                          
  //void  meshSWLISMHayashi(const ProblemDomain & a_prob_domain,
  //                        const Vector<Real>  & a_domainBox);
                          
  /// For a tracking surface 'a_s' returns first component of the velocity field for the level set method
  /**
   */       
  virtual int lsIndexField(int a_s);
  
  
  /// Things to do before after calculations are complete
  /**
   */
  virtual void conclude(const LevelData<FArrayBox>& a_U,
                        Real a_time,
                        int  a_cur_step,                        
                        const std::string & a_inputfile);
                        
                        
  /// Calculate problem specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                                  BaseFab<int> & a_REG,
                               const Real      & a_dt,                               
                               const Box       & a_box);
                               
                               
                               
  void prepareForH5Writing();
  
  void writeSphericalSlice(const LevelData<FArrayBox>& a_U, Real a_time);                               
  void shockBC( const FArrayBox & a_W,
                      FArrayBox & a_U,
                      BaseFab<int> & a_REG,
                      const int   a_level);



protected:
                                // True if the Fortran common block has been set
  bool m_isFortranCommonSet;
  
  // Solar wind parameters    
  Real  m_sunN;
  Real  m_sunV;
  Real  m_sunT;
  Real  m_sunB;
  Real  m_sunTILT;  
  Real  m_initR;
  Real  m_dirBrNorth;
  Real  m_sunIntBCRadius;
  int   m_sunBmonopolar;

  // LISM parameters
  Real  m_lismN;
  Real  m_lismM;
  Real  m_lismV;
  Real  m_lismUX;
  Real  m_lismUY;
  Real  m_lismUZ;
  Real  m_lismT;
  Real  m_lismB;
  Real  m_lismBX;
  Real  m_lismBY;
  Real  m_lismBZ;
  Real  m_netN;
  Real  m_Bref;
 
  Real  m_gamma;
  
  Real  m_TMLIM; // Temperature threshold between 1 and 2 regions

  // Turbulence model parameters
  Real  m_sunZ2;
  Real  m_sunSigmaC;
  Real  m_sunLambda;
  Real  m_lismZ2;
  Real  m_lismSigmaC;
  Real  m_lismLambda;

  // Pickup Ions model parameters
  Real  m_sunN_PI;
  Real  m_sunT_PI;
  Real  m_lismN_PI;
  Real  m_lismT_PI;
  
  bool m_photoionize; 
  bool m_const_H;
  bool m_region_tracer;


  /////////////////////////////////////////////////////////
  ////// Parameters for Casrtesian coordinate system //////
  
  Real  m_R0;

  Real  m_adaptRmin;
// If m_adaptRmax > 0 then any adaptation is prohibited outside r < m_adaptRmax region
//                                                           !!! Only for Cartesian CS
  Real  m_adaptRmax;

//                                      Selector for specialized adaptation procedures
//                                                           !!! Only for Cartesian CS
  int   m_helioAdapt;

//   If m_adaptXL1 < 0 then any adaptation is prohibited for x < m_adaptXL1 on level 1
//                                                           !!! Only for Cartesian CS
  Real  m_adaptXL1;

//   If m_adaptXL2 < 0 then any adaptation is prohibited for x < m_adaptXL2 on level 2
//                                                           !!! Only for Cartesian CS
  Real  m_adaptXL2;
  
//       If a point lies inside of m_R0dtIgnore it doesn't diminish time step of level
  Real m_R0dtIgnore;

// Second and third regions can be defined only in strip [m_RegSafeZBot,m_regSafeZTop]
//                                                           !!! Only for Cartesian CS
  Real m_RegSafeZTop;
  Real m_RegSafeZBot;

  // Node coordinates of the Sun. 
  // Sun is located at the node between
  // m_sunIJK and m_sunIJK-1 cells  
  IntVect  m_sunIJK;    
  RealVect m_sunXYZ;  
  
  // Grid must refined within R0
  Box   m_R0Box;  
  
  // Box where the next kinetic level lies (grid is refined within thos Box)
  Box   m_kineticBox;
    
  ////// Parameters for Casrtesian coordinate system //////
  /////////////////////////////////////////////////////////
  
    
  /////////////////////////////////////////////////////////
  ////// Parameters for spherical coordinate system //////
  
  // Data for path-through BC for neutrals
  LevelData<FArrayBox> m_lstrips;
  LevelData<FArrayBox> m_ladjstr;  
    
  
  ////// Parameters for spherical coordinate system //////
  /////////////////////////////////////////////////////////
  
  int m_subproblem;
  
  // true if take into account the Sun gravity
  bool m_bSunGravity;
  
  // true if take into account heating/acceleration of the solar wind
  bool m_bSunHeating;
  
  // true, if we need to plot region
  bool m_output_region;
  
  // true, if we need output vectors in spherical coordinates when we use cartesian coordinates
  bool m_output_vecSPH;
  
  
  // Data and methods for writing spherical distributions to a file
  bool m_writeH5;
  std::string m_writeH5File;
  Real m_writeT0;        // dimensionless time to start  
  Real m_writeDt;        // dt of the data in days
  Real m_lastWrittenTime;
  int  m_lastInd;        // index of the last dataset that present in the file
  int  m_concludeI;      // index in 'r' array of 'R0' parameter
  Real m_concludeR0;     // Radius where output data
  bool m_writeRotating;  // true if we write in rotating frame
  
  // Level set indices
  struct {
  int  m_iHCS;         // index of HCS tracer
  int  m_iHCSb;        // index for boundary HCS tracer
  int  m_iRegTr;       // index for region tracer
  } m_ls_indices;


private:
                                         // Disallowed for all the usual reasons
  void operator=(const HeliosphericProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  HeliosphericProblem(const HeliosphericProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif

#endif /*  HELIOSPHERICPROBLEM_H_ */
