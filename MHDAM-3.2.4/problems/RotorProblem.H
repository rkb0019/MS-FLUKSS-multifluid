#ifndef  ROTORPROBLEM_H_
#define  ROTORPROBLEM_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysProblem.H"

/// Example of planar Riemann problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class RotorProblem : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  RotorProblem();
  
  /// Destructor
  /**
   */
  ~RotorProblem() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Returns physical model used
  /**
   */
  virtual ePhysicalModel physicalModel( void ) const { return PP_MHDPM;}

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real&     a_gamma,
                         const Real&     a_rotorU,
                         const Real&     a_rotorP,
                         const Real&     a_rotorR,
                         const Real&     a_rotorB,
                         const Real&     a_rotorXC,
                         const Real&     a_rotorYC,
                         const Real&     a_rotorR0,
                         const Real&     a_rotorR1 );
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
  void fluxBC(       FArrayBox&      a_F,
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

  Real  m_rotorU;
  Real  m_rotorP;
  Real  m_rotorR;
  Real  m_rotorB;
  Real  m_rotorXC;
  Real  m_rotorYC;
  Real  m_rotorR0;
  Real  m_rotorR1;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const RotorProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  RotorProblem(const RotorProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
