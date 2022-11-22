//  PhysProblem.H
// ============
//
// Virtual base class through which a user specifies the initial and boundary
// conditions for a hyperbolic system of PDEs.
//

#ifndef _PHYSPROBLEM_H
#define _PHYSPROBLEM_H

#include "FArrayBox.H"
#include "ParmParse.H"
#include "REAL.H"
#include "LevelData.H"
#include "ProblemDomain.H"
#include "IntVectSet.H"

#include "CSHandler.H"
#include "EquationSystem.H"

class SourceCalculator;
class RiemannSolver;

/// Physical problem stuff (input, initial and boundary conditions etc)
/**
   Virtual base class through which a user specifies physical problem dependent stuff
   input parameters, set the initial and boundary conditions etc.
 */
class PhysProblem
{
public:
  enum ePhysicalModel {
                        PP_Undefined   = 0,
                        PP_EulerPM     = 1,
                        PP_MHDPM       = 2,
                        PP_2FluidPM    = 3,
                        PP_3FluidPM    = 4,
                        PP_4FluidPM    = 5,
                        PP_5FluidPM    = 6
                      };

  enum eBoundaryConditions {
                        BC_Undefined   = 0,
                        BC_Periodic    = 1,
                        BC_Fixed       = 2,
                        BC_Continuous  = 3,
                        BC_Axis        = 4,
                        BC_Mixed       = 5
                      };
                      
                        
                        

  /// Constructor
  /**
   */
  PhysProblem();

  /// Destructor
  /**
   */
  virtual ~PhysProblem() {};

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this object.
   */
  virtual void define(const ProblemDomain& a_domain,                      
                      const int            a_level);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual PhysProblem* new_PhysProblem() = 0;

  /// Copy method
  /**
     Copy all data from a_PP to this instance.
   */
  virtual void copy_PhysProblem(const PhysProblem* a_PP);

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity ) = 0;

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
                          
  /// Modifies primitives variables for writing them to plot file (e.g. make them dimensional)
  /**
  */
  virtual void primForPlot(      FArrayBox& a_W,
                           const Box&       a_box);
                           
  /// Creates "AUXDATA" fields for tecplot
  /**
    atype = 0 DATASETAUXDATA fields are created
    atype = 1 AUXDATA fields are created (i.e. data are associated with a zone)
  */
  virtual void auxDataTecplot(std::string & a_str,
                              Real          a_time,
                              int           a_type); 



  /// Returns physical model used
  /**
   */
  virtual ePhysicalModel physicalModel( void ) const { return PP_Undefined;}

  /// Fill ghost cells
  /**
   */
  virtual void fillGhostCells(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time ) = 0;

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(       FArrayBox&      a_F,
                             FArrayBox&      a_Bn,
                       const FArrayBox&      a_W,
                       const FArrayBox&      a_Wextrap,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const Real&           a_time ) = 0;

  /// Apply artificial viscosity at the boundary
  /**
     The default implementation of this currently does nothing.
   */
  virtual void artViscBC(       FArrayBox& a_F,
                          const FArrayBox& a_U,
                          const FArrayBox& a_divVel,
                          const int&       a_dir,
                          const Real&      a_time );

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;
  
  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U,
                                     Interval & a_comp);

  /// Calculate problem specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                                  BaseFab<int> & a_REG,
                               const Real      & a_dt,                               
                               const Box       & a_box);

  /// Problem specific postprocessing
  /**
   */
  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const Real      & a_time,                               
                               const Box       & a_box);

  /// Shock boundary conditions
  virtual void shockBC( const FArrayBox & a_W,
                              FArrayBox & a_U,
                           BaseFab<int> & a_REG,
			      const int   a_level);

  /// Set the source calculator pointer
  /**
   */
  virtual void setSourceCalculator( SourceCalculator * a_pSC );

  /// Get the source calculator pointer
  /**
   */
  virtual SourceCalculator * getSourceCalculator( void ) const;

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

  /// Set the Coordinate System flag
  /**
   */
  virtual void setCoordinateSystem( CoordinateSystemHandler* a_CSH );

  /// Get the Coordinate System flag
  /**
   */
  virtual CoordinateSystemHandler * coordinateSystem( void ) const;
  
  /// Set the Equation System flag
  /**
   */
  virtual void setEquationSystem( EquationSystem* a_eqSys );

  /// Get the Equation System flag
  /**
   */
  virtual EquationSystem * equationSystem( void ) const;
  
  /// Set the number of ghost cells
  /**
   */
  void setNumGhostCells(int a_numGhostCells);
  
  
  /// Things to do before advancing one level by one time step.
  /**
   */
  virtual void preTimeStep(LevelData<FArrayBox>&  a_U, Real a_time);

  /// Things to do after advancing this level by one time step.
  /**
   */
  virtual void postTimeStep(LevelData<FArrayBox>& a_U);
  

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

  /// Return boundary condition flags for all boundaries
  /**
  depreciated
   */
  virtual void getBCFlags( eBoundaryConditions leftBC,
                           eBoundaryConditions rightBC,
                           eBoundaryConditions bottomBC,
                           eBoundaryConditions topBC,
                           eBoundaryConditions frontBC,
                           eBoundaryConditions behindBC );

  /// Return the boundary condition flag for given boundary
  /**
  depreciated
   */
  virtual eBoundaryConditions getBCFlags( int a_dir, Side::LoHiSide a_sd );
  
  /// Boundary conditions for divB cleaning using the projection scheme
  virtual void projectBbc(FArrayBox &a_state, const Box &a_valid, const ProblemDomain &a_domain, Real a_dx, bool a_homogeneous);


  
  /// Compute problem specific dt and returns a cell with the minimum dt. 
  /**
   */                               
  virtual Real computeDt( const FArrayBox& a_U,
                          const FArrayBox& a_dt,
                          const Box&     a_box,
                          IntVect&       a_minDtCell);
                          
  
  /// Modifies curvilinear mesh (applies to curvilinear meshes only)
  /**
   */                               
  virtual void defineMesh(const ProblemDomain & a_prob_domain,
                          const Vector<Real>  & a_domainBox);
  
  /// Converts dimensionless time to problem specific time
  /**
   */                               
  virtual Real getPhysTime(Real a_time); 
  
  /// For a tracking surface 'a_s' returns first component of the velocity field for the level set method
  /**
   */       
  virtual int lsIndexField(int a_s);
    
      
  /// Returns file header for 1D probes and number of variables
  /**
   */                        
  virtual void probeFilesHeader(int & a_numVars, std::string& a_header);
  
  
  /// Converts 'probe_dt' value in inputs file to dimensionless units
  /**
   */                        
  virtual Real probeFilesDimensiolessDt(Real a_dt);
                             
                             
  /// Things to do before after calculations are complete
  /**
   */
  virtual void conclude(const LevelData<FArrayBox>& a_U,
                        Real a_time,
                        int  a_cur_step,
                        const std::string & a_inputfile);
    

protected:
  // define() has been called
  bool          m_isDefined;

  // The current level's problem domain
  ProblemDomain m_domain;
  

  // Index of current level
  int           m_level;

  // Calculation of source terms for this physical problem
  SourceCalculator * m_pSourceC;

  // MHD Riemann solver object
  // It is needed to solve Riemann problems on boundaries.
  RiemannSolver * m_RS;

  // Gas Dynamics Riemann solver object and whether it has been set
  // It is needed for multifluid problems.
  RiemannSolver * m_RSGD;

  // Coordinate system handler
  CoordinateSystemHandler* m_csh;
  
  // System of equations
  EquationSystem* m_eqSys;
  
  int m_verbosity;
  
  // Number of ghost cells
  int m_numGhostCells;


private:
  // Disallowed for all the usual reasons
  void operator=(const PhysProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PhysProblem(const PhysProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
