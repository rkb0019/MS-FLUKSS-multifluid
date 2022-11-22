#ifndef  BLASTPROBLEM_H_
#define  BLASTPROBLEM_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysProblem.H"

/// Blast interaction problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class BlastProblem : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  BlastProblem();
  
  /// Destructor
  /**
   */
  ~BlastProblem() {};

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
                         const Real&     a_Blastp0,
                         const Real&     a_Blastp1,
                         const Real&     a_BlastXc,
                         const Real&     a_BlastYc,             
                         const Real&     a_BlastZc,             
                         const Real&     a_BlastBx,
                         const Real&     a_BlastBy,
                         const Real&     a_Blastr0 );
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
  virtual void initialize(LevelData<FArrayBox>& a_U);

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

  Real  m_Blastp0;
  Real  m_Blastp1;
  Real  m_BlastXc;
  Real  m_BlastYc;
  Real  m_BlastZc;      
  Real  m_BlastBx;
  Real  m_BlastBy;
  Real  m_Blastr0;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const BlastProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  BlastProblem(const BlastProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
