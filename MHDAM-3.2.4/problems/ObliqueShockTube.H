#ifndef  OBLIQUESHOCKTUBE_H_
#define  OBLIQUESHOCKTUBE_H_

#include  <iostream>

#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"

#include "PhysProblem.H"


/// Example of planar Riemann problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class ObliqueShockTube : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  ObliqueShockTube();
  
  /// Destructor
  /**
   */
  ~ObliqueShockTube() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Returns physical model used
  /**
   */
  ePhysicalModel physicalModel( void ) const { return PP_MHDPM;}

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

  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  PhysProblem *new_PhysProblem();

  /// Fill ghost cells
  /**
   */
  virtual void fillGhostCells(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );

  /// Set boundary fluxes
  /**
   */
  void fluxBC(      FArrayBox&      a_F,
                    FArrayBox&      a_Bn,
              const FArrayBox&      a_W,
              const FArrayBox&      a_Wextrap,
              const int&            a_dir,
              const Side::LoHiSide& a_side,
              const Real&           a_time);
  
  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);

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

  
protected:
                                // True if the Fortran common block has been set
  bool m_isFortranCommonSet;

  Real  m_densityL;
  Real  m_densityR;
  Real  m_pressureL;
  Real  m_pressureR;
  Real  m_velxL;
  Real  m_velxR;
  Real  m_velyL;
  Real  m_velyR;
  Real  m_velzL;
  Real  m_velzR;
  Real  m_BxL;
  Real  m_BxR;
  Real  m_ByL;
  Real  m_ByR;
  Real  m_BzL;
  Real  m_BzR;  
  Real  m_gamma;
    
  Real  m_startX;    
  Real  m_startY;
  Real  m_tanangle;
  
  
  int m_subproblem;
  

private:
                                         // Disallowed for all the usual reasons
  void operator=(const ObliqueShockTube& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  ObliqueShockTube(const ObliqueShockTube& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
