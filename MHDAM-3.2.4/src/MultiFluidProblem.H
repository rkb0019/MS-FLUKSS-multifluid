#ifndef _MULTIFLUIDPROBLEM_H
#define _MULTIFLUIDPROBLEM_H

#include "PhysProblem.H"

/// Physical problem stuff (input, initial and boundary conditions etc)
/**
   Virtual base class through which a user specifies physical problem dependent stuff
   input parameters, set the initial and boundary conditions etc.
 */
class MultiFluidProblem : public PhysProblem
{
public:
  /// Constructor
  /**
   */
  MultiFluidProblem();

  /// Destructor
  /**
   */
  virtual ~MultiFluidProblem() {};

  /// Define the object
  /**
     Set the problem domain index space and the grid spacing for this object.
   */
  virtual void define(const ProblemDomain& a_domain,
                      const int            a_level);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new MultiFluidProblem object with m_isDefined = false (i.e.,
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

  /// Returns physical model used
  /**
   */
  virtual ePhysicalModel physicalModel( void ) const { return m_physModel;}

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(      FArrayBox&      a_F,
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
  virtual void artViscBC(      FArrayBox& a_F,
                         const FArrayBox& a_U,
                         const FArrayBox& a_divVel,
                         const int&       a_dir,
                         const Real&      a_time);

  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U) = 0;  

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

  /// Creates tagged cells for dynamic mesh refinement
  /**
    Problem specific cells tagging
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags);
  
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
                              const Box          & a_box ) = 0;
  
  /// Return number of fluids
  /**
   */
  virtual int nFluids();

protected:
  ePhysicalModel m_physModel;

private:
  // Disallowed for all the usual reasons
  void operator=(const MultiFluidProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  MultiFluidProblem(const MultiFluidProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
