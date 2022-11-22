#ifndef _RIEMANN_HD_H_
#define _RIEMANN_HD_H_

#include "RiemannSolver.H"

/// Simple Riemann solver for hydrodynamics from CHOMBO.
/**
 */
class RiemannHD : public RiemannSolver
{
public:
  /// Constructor
  /**
   */
  RiemannHD( void );

  /// Destructor
  /**
   */
  virtual ~RiemannHD( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void );

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
                       const Box &       a_box);

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
                              const Box &       a_box);

  /// Availability of middle state at the faces
  virtual bool isMiddleStateAvailable( void )
  {
    return true;
  }

  virtual RS_SYSTEM getSystem( void ){ return RS_EULER; }

private:
  // Disallowed for all the usual reasons
  void operator=(const RiemannHD& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  RiemannHD(const RiemannHD& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _RIEMANN_HD_H_ */

