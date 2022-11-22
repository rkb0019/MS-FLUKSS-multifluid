#ifndef  SWLISMProblem_H_
#define  SWLISMProblem_H_


#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "MultiFluidProblem.H"

#if CH_LANG_CC
// Boundary condition constants
#endif

#define SW_TURB     1
#define SW_LONGTAIL 2
#define SW_KINETIC  3


class PatchMHDAM;

/// SW and LISM interaction problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class SWLISMProblem : public MultiFluidProblem
{ 
public:
  /// Null constructor
  /**
   */
  SWLISMProblem();
  
  SWLISMProblem(ePhysicalModel ePM);
        
  /// Destructor
  /**
   */
  ~SWLISMProblem() {};
    
  
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );
  
  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon();
  
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


      
  /// Fill ghost cells
  /**
   */  
  virtual void fillGhostCells(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  /// Fill ghost cells (default problem)
  /**
   */
  void fillGhostCellsDefault(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
                               
  /// Fill ghost cells
  /**
   */
  void fillGhostCellsZankTurb2007(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(       FArrayBox&      a_F,
                             FArrayBox&      a_Bn,
                       const FArrayBox&      a_WMinus,
                       const FArrayBox&      a_WPlus,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const Real&           a_time);
                       
  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);
  
  /// Set up initial conditions
  /**
   */
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
                               
  /// Things to do after advancing one level by one time step.
  /**
   */
  virtual void postTimeStep(LevelData<FArrayBox>&       a_U);                             

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
                               
  virtual void defineRegions( const FArrayBox    & a_W,
                                    FArrayBox    & a_S,
                                    BaseFab<int> & a_R,
                              const Box          & a_box);
                              
  /// Calculate problem specific source terms
  /**
   */                            
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                                  BaseFab<int> & a_REG,
                               const Real      & a_dt,                               
                               const Box       & a_box);
                              
                              
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
                           
  /// Creates tagged cells for dynamic mesh refinement
  /**
    Problem specific cells tagging
   */
  virtual void tagCells(const FArrayBox&  a_U,
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
            


protected:
  // True if the Fortran common block has been set
  bool m_isFortranCommonSet;


  // Solar wind parameters  
  RealVect m_sunXYZ;  
  Real     m_R0;
  Real     m_sunN;
  Real     m_sunV;
  Real     m_sunT ;
  Real     m_initR;

  // Lism parameters
  Real     m_lismN;
  Real     m_lismM;
  Real     m_lismV;
  Real     m_lismT;
  Real     m_lismB;
  Real     m_netN;  // LISM neutral density
  
  Real     m_gamma;
  
  bool m_photoionize; 
  bool m_const_H;
  
  // If x > stop_ref then any refinement is blocked
  Real     m_stop_ref;  
  // Finest level for the TS.
  int      m_max_levelTS;
  
  Real     m_TMLIM;  // Temperature separating regions 1 and 2  
  
  // If a point lies inside of m_R0dtIgnore it doesn't diminish time step of level
  Real m_R0dtIgnore;
  // Second and third regions can be defined only in strip [-m_RegSafeZ,+m_regSafeZ]
  Real m_RegSafeZ;

  // Node coordinates of the Sun. 
  // Sun is located at the node between
  // m_sunIJK and m_sunIJK-1 cells  
  IntVect m_sunIJK;    
  
  // Grid must refined within R0
  Box   m_R0Box;
    
  Box           m_gridBox;         // This box should be tagged during regridding.
                                   // Attention !!! This box doesn't coincide with 
                                   // any of KineticSources2D.m_gridsBoxes  
  
  int m_subproblem;   

  bool m_region_tracer;
 
  // true, if we need to plot region
  bool m_output_region;
  
  // Data for path-through BC for neutrals
  LevelData<FArrayBox> m_lstrips;
  LevelData<FArrayBox> m_ladjstr;  
  
  //  === Turbulence 2D for Zank, Hawaii 2007 ===    
  
  Real m_dTimeFactor;  // Time in years = m_dTimeFactor * TimeMHDAM

  // Turbulence time steps in dimensionless units
  Real m_dTurbLISMTimeStep;  
  Real m_dTurbSWTimeStep;  
  
  // Diviation from average values
  Real m_dLISMDeviation;
  Real m_dSWDeviation;
  
  Real m_dTime;
  Real m_R0Turb;
    
  FArrayBox m_Phi;
  FArrayBox m_LISMTurb; // Left ghost cells for all patches on this CPU
        
  Real m_dNextTurbLISMTime;
  Real m_dNextTurbSWTime;
  
  void CalculatePhi(const LevelData<FArrayBox>& a_U,
                          PatchMHDAM& a_patch);
  // === Turbulence 2D for Zank, Hawaii 2007 ===    
    
  
};

#endif
