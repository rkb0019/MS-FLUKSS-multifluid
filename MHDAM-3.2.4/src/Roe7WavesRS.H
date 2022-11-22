#ifndef _ROE_7_WAVES_RS_H_
#define _ROE_7_WAVES_RS_H_

#include "RiemannSolver.H"

/// 7 wave Roe solver for MHD Riemann problem.
/**
 */
class Roe7WavesRS : public RiemannSolver
{
private:
  double m_dEntropyFix;
  double m_dSmallB;
  int    m_iAveraging;
  int    m_iLaxFriedrix;

public:
  /// Constructor
  /**
   */
  Roe7WavesRS( void );

  /// Destructor
  /**
   */
  virtual ~Roe7WavesRS( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void );

  /// Set parameters
  /**
   */
  void setParameters( double EntFix, double dSmallB, int iAver = 0, int iLF = 0 );

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
  void operator=(const Roe7WavesRS& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  Roe7WavesRS(const Roe7WavesRS& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _ROE_8_WAVES_RS_H_ */

