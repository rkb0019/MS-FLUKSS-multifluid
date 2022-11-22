#ifndef _DEDNER_RS_H_
#define _DEDNER_RS_H_

#include "RiemannSolver.H"

/// 8 waves Roe solver for MHD Riemann problem.
/**
 */
class DednerRS : public RiemannSolver
{
private:
  Real m_dEntropyFix;
  Real m_dSmallB;
  int    m_iAveraging;
  int    m_iLaxFriedrix;  
public:
  /// Constructor
  /**
   */
  DednerRS( void );

  /// Destructor
  /**
   */
  virtual ~DednerRS( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new RiemannSolver object.
   */
  virtual RiemannSolver * new_RiemannSolver( void );

  /// Set parameters
  /**
   */
  void setParameters( Real EntFix, Real dSmallB, int iAver = 0, int iLF = 0 );

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

  virtual RS_SYSTEM getSystem( void ){ return RS_DEDNER; }

private:
  // Disallowed for all the usual reasons
  void operator=(const DednerRS& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  DednerRS(const DednerRS& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _ROE_8_WAVES_RS_H_ */

