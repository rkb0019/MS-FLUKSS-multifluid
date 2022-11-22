#ifndef _RUSANOV_RS_H_
#define _RUSANOV_RS_H_

#include "RiemannSolver.H"

/// Rusanov solver for MHD Riemann problem.
/**
 */
class RusanovRS : public RiemannSolver
{
private:
  double m_dSmallB;
  int    m_iAveraging;

public:
  /// Constructor
  /**
   */
  RusanovRS( void );

  /// Destructor
  /**
   */
  virtual ~RusanovRS( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void );

  /// Set parameters
  /**
   */
  void setParameters( double dSmallB = 1.0e-9, int iAver = 0 );

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
                              const Box &       a_box)
  {
    MayDay::Error("This function is not implemented");
  }

  virtual RS_SYSTEM getSystem( void ){ return RS_MHD; }

private:
  // Disallowed for all the usual reasons
  void operator=(const RusanovRS& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  RusanovRS(const RusanovRS& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _RUSANOV_RS_H_ */

