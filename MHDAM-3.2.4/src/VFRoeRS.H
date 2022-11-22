#ifndef _VFROE_RS_H_
#define _VFROE_RS_H_

#include "RiemannSolver.H"

/// 8 waves VFRoe solver for MHD Riemann problem.
/**
 */
class VFRoeRS : public RiemannSolver
{
private:
  double m_dEntropyFix;
  double m_dSmallB;
  int    m_iAveraging;

public:
  /// Constructor
  /**
   */
  VFRoeRS( void );

  /// Destructor
  /**
   */
  virtual ~VFRoeRS( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void );

  /// Set parameters
  /**
   */
  void setParameters( double EntFix, double dSmallB, int iAver = 0 );

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

  virtual RS_SYSTEM getSystem( void ){ return RS_MHD; }

private:
  // Disallowed for all the usual reasons
  void operator=(const VFRoeRS& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  VFRoeRS(const VFRoeRS& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _VFROE_RS_H_ */

