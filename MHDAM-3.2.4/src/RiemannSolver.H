#ifndef _RIEMANN_SOLVER_H_
#define _RIEMANN_SOLVER_H_

#include "FArrayBox.H"
#include "Box.H"
#include "PhysProblem.H"

//class PhysProblem;

/// Base class for solution of a Riemann problem and calculation fluxes at the faces.
/**
 */
class RiemannSolver
{
public:

  enum RS_SYSTEM {
                   RS_UNDEFINED = 0,
                   RS_EULER     = 1,
                   RS_MHD       = 2,
                   RS_DEDNER    = 3
                 };

  /// Constructor
  /**
   */
  RiemannSolver( void ){};

  /// Destructor
  /**
   */
  virtual ~RiemannSolver( void ){};

  /// Static method for construction new objects
  /**
   */
  static RiemannSolver * make( int iFlux, PhysProblem::ePhysicalModel eModel, Real dEntFix );

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void ) = 0;

  /// Compute a Riemann problem and generate fluxes at the faces
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate fluxes at the faces within a_box.
   */
  virtual void fluxes(       FArrayBox & a_F,
                       const FArrayBox & a_WLeft,
                       const FArrayBox & a_WRight,
                       const int &       a_dir,
                       const int &       a_iRho,
                       const Box &       a_box) = 0;

  /// Compute a Riemann problem and generate middle states at the faces
  /**
     Given input left and right states in a direction, a_dir, compute a
     Riemann problem and generate middle states at the faces within a_box.
   */
  virtual void statesAtFaces(       FArrayBox & a_F,
                              const FArrayBox & a_WLeft,
                              const FArrayBox & a_WRight,
                              const int &       a_dir,
                              const int &       a_iRho,
                              const Box &       a_box) = 0;

  /// Availability of middle state at the faces
  virtual bool isMiddleStateAvailable( void )
  {
    return false;
  }

  virtual RS_SYSTEM getSystem( void ){ return RS_UNDEFINED; }

private:
  // Disallowed for all the usual reasons
  void operator=(const RiemannSolver& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  RiemannSolver(const RiemannSolver& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _RIEMANN_SOLVER_H_ */

