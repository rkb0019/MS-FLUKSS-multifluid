#ifndef _PATCH_MHDAM_H_
#define _PATCH_MHDAM_H_

#include <FluxBox.H> 

#include "EdgeBox.H"
#include "PhysProblem.H"
#include "RiemannSolver.H"
#include "Reconstruction.H"
#include "RefCriteria.H"

#include <string>
using std::string;

class CoordinateSystemHandler;
class EquationSystem;

///
/**
   The base class PatchMHDAM provides an implementation of a second-order,
   unsplit Godunov method acting on a single grid/patch.  PatchMHDAM
   provides an interface to the level integrator, LevelMHDAM, which manages
   the entire level and flux corrections (via flux registers) to the coarser
   and finer levels.  In addition, the physics dependent code is not provided
   in PatchMHDAM but is supplied by the user by subclassing PatchMHDAM and
   implementing the pure virtual functions.  Some parameters can also be
   adjusted to modify the algorithm.  All functions are virtual so any of them
   can be reimplemented by the user.

   There are three types of grid variables that appear in the unsplit Godunov
   method: conserved variables, primitive variables, and fluxes, denoted below
   by U, W, F, respectively.  It is often convenient to have the number of
   primitive variables and fluxes exceed the number of conserved variables.
   In the case of primitive variables, redundant quantities are carried that
   parameterize the equation of state in order to avoid multiple calls to that
   the equation of state function.  In the case of fluxes, it is often
   convenient to split the flux for some variables into multiple components,
   e.g., dividing the momentum flux into advective and pressure terms.  The
   API given here provides the flexibility to support these various
   possibilities.
 */

enum eDivergenceCleaning {
  DC_divBIgnore  = 0,
  DC_CT_BS       = 1,  // Constrained transport, E is calculated by Balsara & Spacer method
  DC_CT_GS0      = 2,  // Constrained transport, E is calculated by Gardiner & Stone method 1
  DC_CT_GS1      = 3,  // Constrained transport, E is calculated by Gardiner & Stone method 2
  DC_ProjHancock = 4,  // Projection scheme with Hancock
  DC_CP          = 5   // Correction potential (Dedner approach)
};

 
class PatchMHDAM
{
public:

  enum eTimeApproximation {
            TA_FirstOrder = 0,
            TA_Hancock    = 1,
            TA_RK2        = 2
       };
                      
  /// Constructor
  /**
   */
  PatchMHDAM();

  /// Destructor
  /**
   */
  virtual ~PatchMHDAM();
  
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Define the object
  /**
   */
  virtual void define(ProblemDomain& a_domain,                      
                      const int      a_level,
                      const int      a_verbosity);

  /// Set the initial and boundary condition object
  /**
   */
  virtual void setPhysProblem(PhysProblem* a_bc);

  /// Get the initial and boundary condition object
  /**
   */
  virtual PhysProblem* getPhysProblem() const;

  /// Copy internal data to another PatchMHDAM
  /**
   */
  virtual void copyTo( PatchMHDAM * pPatch ) const;

  /// Set the MHD Riemann solver object
  /**
   */
  virtual void setRiemannSolver( RiemannSolver * a_RS );

  /// Get the MHD Riemann solver object
  /**
   */
  virtual RiemannSolver * getRiemannSolver( void ) const;

  /// Set the Gas Dynamics Riemann solver object
  /**
   */
  virtual void setRiemannSolverGD( RiemannSolver * a_RS );

  /// Get the Gas Dynamics Riemann solver object
  /**
   */
  virtual RiemannSolver * getRiemannSolverGD( void ) const;

  /// Set the Reconstruction object
  /**
   */
  virtual void setReconstruction( Reconstruction * a_RS );

  /// Get the Reconstruction object
  /**
   */
  virtual Reconstruction * getReconstruction( void ) const;

  /// Set the Time Approximation flag
  /**
   */
  virtual void setTimeApproximation( eTimeApproximation method );

  /// Get the Time Approximation flag
  /**
   */
  virtual eTimeApproximation timeApproximation( void ) const;

  /// Set the Divergence Cleaning flag
  /**
   */
  virtual void setDivergenceCleaning( eDivergenceCleaning method );

  /// Get the Divergence Cleaning flag
  /**
   */
  virtual eDivergenceCleaning getDivergenceCleaning( void ) const;

  /// Set the div(B) calculation method flag
  /**
   */
  virtual void setDivBMethod( int method );

  /// Get the div(B) calculation method flag
  /**
   */
  virtual int getDivBMethod( void ) const;
  
  
  /// Set the div(B) calculation method flag
  /**
   */
  virtual void set8waveFlag( int a_8waveUse );
  
  /// Set m_bLSonly flag
  /**
   */
  virtual void setLSonlyFlag( int a_LSonly );
  
  /// Get m_bLSonly flag
  /**
   */
  virtual bool getLSonlyFlag( void ) const;

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new PatchMHDAM object with its initial and boundary
     condtions, slope parameters, and artificial viscosity information defined.
   */
  virtual PatchMHDAM* new_patchMHDAM() const = 0;

  /// Get the source calculator pointer
  /**
   */
  virtual SourceCalculator * getSourceCalculator( void ) const;

  /// Set the current time before calling updateState()
  /**
   */
  virtual void setCurrentTime(const Real& a_currentTime);

  /// Things to do after advancing this level by one time step.
  /**
   */
  virtual void postTimeStep( void );
  
  /// One time step usinf first order scheme
  /**
   */
  void EulerStep(                   FArrayBox& a_U,
                                    FArrayBox  a_F[CH_SPACEDIM],
                                    EdgeBox  & a_E,    
                                    FArrayBox& a_divB,
                              const FArrayBox& a_scale,
                              const FArrayBox  a_areas[CH_SPACEDIM],
                              const FArrayBox& a_S,                              
                              const Real&      a_dt,
                              const Box&       a_box);
                              
                              
  void predictorHancock(            FArrayBox& a_U,                                    
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
                              const Box&       a_box);                              
                                  
                              
  void correctorHancock(            FArrayBox& a_U,
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
                              const Box&       a_box);
                              


  /// Update the conserved variables and return the final fluxes used for this
  /**
     Update the conserved variables and return the final fluxes that were used
     for this.  Compute the fluxes using a second-order, unsplit Godunov method
     based on the input conserved variables, a_U, and source terms, a_S.  If
     there are no source terms then a_S should be null constructed.  Also
     return the maximum wave speed.
   */
  void updateState(        FArrayBox&       a_U,                           
                           FArrayBox        a_F[CH_SPACEDIM],
                                  EdgeBox   & a_E,
                           Real&            a_newDt,
                           IntVect&         a_minDtCell,
                           const FArrayBox& a_S,
                           FArrayBox&       a_divB,
                           const Real&      a_dt,
                           const Box&       a_box);
  
  
  /// Compute dt and returns a cell with the minimum dt. 
  /**
   */                               
  virtual Real computeDt( const FArrayBox& a_U,
                          const Box&     a_box,
                          IntVect&       a_minDtCell) = 0;
                          
  /// Compute dt and returns a cell with the minimum dt. 
  /**
   */                               
  virtual Real computeDtLevelSet( const FArrayBox& a_U,
                          const Box&     a_box,
                          IntVect&       a_minDtCell);
  
  /// Number of ghost cells
  /**
  */
  virtual int numGhostCells();
  
  /// ??? 
  /**
   */
  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_Uold,
                               const Real      & a_dt,
                               const Box       & a_box) = 0;
  

  /// Creates tagged cells for dynamic mesh refinement
  /**
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags) = 0;

  /// This method returns interval that includes indices of all variables used in adaptation criteria.
  /**
    Before calling tagCells ghost cells must be filled.
    But we don't need all variables. Indeed, we need only variables
    which are used in adaptation criteria. This method returns interval
    that includes all variables used in adaptation criteria.
   */
  virtual void tagCellVarsInterval(Interval& a_interval) = 0;
  
  /// Calculate magnetic field in cell centers using electric field on edges
  /**
   */
  virtual void recalculateMagneticField(       FArrayBox & a_U,
                                         const FArrayBox & a_Uold,
                                         const EdgeBox   & a_E,
                                         const Real      & a_dt,                          
                                         const Box       & a_box ) = 0;
  
  /// Calculate electric field on edges
  /**
   */
  virtual void calculateElectricField( const FArrayBox & a_W,
                                       const FArrayBox & a_Wold,
                                       const FArrayBox   a_F[CH_SPACEDIM],
                                             EdgeBox   & a_E,
                                       const Real      & a_dt,                                       
                                       const Box       & a_box ) = 0;

  /// ??? 
  /**
   */
  RefineCriterion RefineParams[REF_NUM];

protected:

public: // temporarily
  /// Is the object completely defined
  /**
     Return true if the object is completely defined.
   */
  virtual bool isDefined() const;
    
  
  /// Calculate the primitive variables to the cell faces
  /**
     Calculate the primitive variables, a_W, in the minus and plus direction,
     using the slopes a_dW within a_box.
   */
  /*virtual void valuesOnFaces(FArrayBox&       a_WMinus,
                             FArrayBox&       a_WPlus,
                             const FArrayBox& a_W,
                             const FArrayBox& a_dW,
                             const Real&      a_scale,
                             const Box&       a_box) = 0;*/

  /// ??? 
  /**
   */
  virtual void updateValuesOnFaces(       FArrayBox & a_WMinus,
                                          FArrayBox & a_WPlus,
                                    const FArrayBox & a_W,
                                    const FArrayBox & a_WOld,
                                    const Box&        a_box);


  /// Update the conserved variables using fluxes and a scaling factor
  /**
     Given the fluxes, a_F, in a direction, a_dir, and a scaling factor,
     a_scale, update the conserved variables, a_U, within a_box:

     a_U = a_U - a_scale * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
   */
  virtual void updateCons(FArrayBox&       a_U,
                          const FArrayBox& a_F,
                                Real       a_dt,
                          const FArrayBox& a_invvol,
                          const int&       a_dir,
                          const Box&       a_box);
                          
                          
  /// Update the level set variables using pseudo-fluxes and a scaling factor
  /**  
   */
  virtual void updateLevelSet(FArrayBox&   a_U,
                          const FArrayBox& a_W,                          
                          const FArrayBox& a_F,
                                Real       a_dt,
                          const int&       a_dir,
                          const Box&       a_box);


  /// Perform final update using fluxes
  /**
     Given the fluxes, a_F, in a direction, a_dir, and a scaling factor,
     a_scale, update the conserved variables, a_U, within a_box:

     a_U = a_U - a_dt * a_scale * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
   */
  virtual void finalUpdate(FArrayBox&       a_U,                           
                           const FArrayBox& a_W,
                           const FArrayBox& a_F,
                                 Real       a_dt,                           
                           const FArrayBox& a_scale,
                           const int&       a_dir,
                           const Box&       a_box);

  /// ??? 
  /**
   */
  virtual void updateHancock(              FArrayBox & a_U,
                                     const FArrayBox & a_W,
                                     const FArrayBox & a_FMinus,
                                     const FArrayBox & a_FPlus,
                                     Real              a_dt,
                                     const FArrayBox & a_scale,  
                                     const FArrayBox   a_areas[CH_SPACEDIM],                                                                       
                                     const int       & a_dir,
                                     const Box       & a_box);


  /// ??? 
  /**
   */
  virtual void fluxesHancock(         FArrayBox & a_FMinus,
                                      FArrayBox & a_FPlus,
                                const FArrayBox & a_WMinus,
                                const FArrayBox & a_WPlus,
                                const int &       a_dir,
                                const Box &       a_box) = 0;
                         
  /// ??? 
  /**
   */
  virtual void fluxesRP(       FArrayBox & a_F,
                         const FArrayBox & a_WLeft,
                         const FArrayBox & a_WRight,
                         const int &       a_dir,
                         const Box &       a_box) = 0;
                         
  
  virtual void levelSetFluxes(       FArrayBox & a_F,
                         const FArrayBox & a_WLeft,
                         const FArrayBox & a_WRight,
                         const int &       a_dir,
                         const Box &       a_box);

  /// ??? 
  /**
   */
  virtual void addExplicitSources(       FArrayBox    & a_U,
                                   const FArrayBox    & a_W,
                                   const FArrayBox    & a_SOut,
                                         FArrayBox    & a_S,
                                   const FluxBox      & a_Bn,
                                         FArrayBox    & a_divB,
                                         BaseFab<int> & a_REG,
                                   const Real         & a_dt, 
                                   const FArrayBox    & a_scale,  // dt/volume
                                   const Box          & a_box) = 0;

  virtual void addExplicitSources(       FArrayBox    & a_U,
                                   const FArrayBox    & a_W,
                                   const FArrayBox    & a_SOut,
                                         FArrayBox    & a_S,
                                   const FluxBox      & a_Bn,
                                         FArrayBox    & a_divB,
                                   const FluxBox      & a_Un,
                                         BaseFab<int> & a_REG,
                                   const Real         & a_dt, 
                                   const FArrayBox    & a_scale,  // dt/volume
                                   const Box          & a_box) = 0;

  /// ??? 
  /**
   */
  virtual void preprocessing ( const FArrayBox    & a_W,
                                     FArrayBox    & a_S,
                                     BaseFab<int> & a_R,
                               const Box          & a_box) = 0;

  /// Computes face centered magnetic for calculating divB
  /**
   */
  virtual void computeBn(       FArrayBox & a_Bn,
                          const FArrayBox & a_W,
                          const FArrayBox & a_WMinus,
                          const FArrayBox & a_WPlus ,
                          const int       & a_method,
                          const int       & a_dir,
                          const Box       & a_box);
  
  virtual void computeUn(       FArrayBox & a_Un,
                          const FArrayBox & a_W,
                          const FArrayBox & a_WMinus,
                          const FArrayBox & a_WPlus ,
                          const int       & a_method,
                          const int       & a_dir,
                          const Box       & a_box);
                         
  /// Computes face centered magnetic for calculating divB
  /**
   */
  virtual void correctBn(       FArrayBox & a_Bn,
                          const FArrayBox & a_phi,                          
                          const int       & a_dir,
                          const Box       & a_box){};
                                                    
  virtual void computeDivB(     FArrayBox & a_divB,
                          const FluxBox   & a_Bn,                                                    
                          const FArrayBox & a_W,
                          const Box       & a_box){};         

  // Has define() been called
  bool                m_isDefined;

  // Problem domain and grid spacing
  ProblemDomain       m_domain;  
  int                 m_level;

  // Initial and boundary condition object and whether it has been set
  PhysProblem *       m_PhPr;
  bool                m_isPPSet;

  // MHD Riemann solver object and whether it has been set
  RiemannSolver *     m_RS;
  bool                m_isRSSet;

  // Gas Dynamics Riemann solver object and whether it has been set
  RiemannSolver *     m_RSGD;
  bool                m_isRSSetGD;

  // Reconstruction object and has it been set
  Reconstruction *    m_Reconstruction;
  bool                m_isRecSet;

  // Current time and has it been set
  Real                m_currentTime;
  bool                m_isCurrentTimeSet;

  // Time approximation flag
  eTimeApproximation  m_TimeMethod;

  // Divergence cleaning method
  eDivergenceCleaning m_DivergenceCleaning;  

  // Method to calcualate divB  
  int                 m_iDivBMethod;
  
  // true when 8-wave source terms are used
  bool                m_b8waveUse;
      
  // Coordinate system handler
  CoordinateSystemHandler* m_csh;
  
  // System of equations
  EquationSystem* m_eqSys;
  
  int                 m_verbosity;
  
  // true if we need to advect level set without updating fluid
  bool                m_bLSonly;
  
  

private:
  // Disallowed for all the usual reasons
  void operator=(const PatchMHDAM& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PatchMHDAM(const PatchMHDAM& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
