#ifndef  OBLIQUEPROBLEM_H_
#define  OBLIQUEPROBLEM_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysProblem.H"

/// Oblique shock wave problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class ObliqueProblem : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  ObliqueProblem();
  
  /// Destructor
  /**
   */
  ~ObliqueProblem() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Returns physical model used
  /**
   */
  ePhysicalModel physicalModel( void ) const { return PP_EulerPM;}

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real&     a_gamma,
                         const Real&     a_densityL,
                         const Real&     a_pressureL,
                         const Real&     a_velxL,
                         const Real&     a_velyL,
                         const Real&     a_angleSW,
                         const Real&     a_startX,
                         const Real&     m_obliqueAe,
                         const Real&     m_obliqueAv,
                         const Real&     m_obliquePsi,
                         const Real&     m_obliqueK   );

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
  
protected:
                                // True if the Fortran common block has been set
  bool m_isFortranCommonSet;

  Real  m_densityL;
  Real  m_pressureL;
  Real  m_velxL;
  Real  m_velyL;
  Real  m_angleSW;
  Real  m_startX;
  Real  m_obliqueAe;
  Real  m_obliqueAv;
  Real  m_obliquePsi;
  Real  m_obliqueK;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const ObliqueProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  ObliqueProblem(const ObliqueProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
