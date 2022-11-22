#ifndef  CLOUDPROBLEM_H_
#define  CLOUDPROBLEM_H_

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
class CloudProblem : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  CloudProblem();
  
  /// Destructor
  /**
   */
  ~CloudProblem() {};

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
                         const Real&     a_cloudRL,
                         const Real&     a_cloudUL,
                         const Real&     a_cloudVL,
                         const Real&     a_cloudWL,
                         const Real&     a_cloudPL,
                         const Real&     a_cloudBXL,
                         const Real&     a_cloudBYL,
                         const Real&     a_cloudBZL,
                         const Real&     a_cloudRR,
                         const Real&     a_cloudUR,
                         const Real&     a_cloudPR,
                         const Real&     a_cloudBYR,
                         const Real&     a_cloudBZR,
                         const Real&     a_cloudXS,
                         const Real&     a_cloudXC,
                         const Real&     a_cloudYC,
                         const Real&     a_cloudZC,
                         const Real&     a_cloudR0,
                         const Real&     a_cloudRho );

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
               const FArrayBox&      a_WMinus,
               const FArrayBox&      a_WPlus,
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
  Real  m_densityR;
  Real  m_pressureL;
  Real  m_pressureR;
  Real  m_velxL;
  Real  m_velxR;
  Real  m_velyL;
  Real  m_velzL;
  Real  m_BxL;
  Real  m_ByL;
  Real  m_ByR;
  Real  m_BzL;
  Real  m_BzR;
  Real  m_startX;

  Real  m_cloudRho;
  Real  m_XC;
  Real  m_YC;
  Real  m_ZC;
  Real  m_R0;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const CloudProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  CloudProblem(const CloudProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
