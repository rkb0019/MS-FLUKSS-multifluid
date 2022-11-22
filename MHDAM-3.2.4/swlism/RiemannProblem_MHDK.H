#ifndef  RIEMANNPROBLEM_MHDK_H_
#define  RIEMANNPROBLEM_MHDK_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"


#include "SolidF_F.H"
#include "PhysProblem.H"

/// Example of planar Riemann problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class RiemannProblem_MHDK : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  RiemannProblem_MHDK();
  
  /// Constructor which defines parameters used by Fortran routines
  /**
   */
  RiemannProblem_MHDK( const Real&     a_gamma,
                       const Real&     a_numdenL,
                       const Real&     a_temperatureL,
                       const Real&     a_velxL,
                       const Real&     a_BxL,          
                       const Real&     a_netnumL,
                       const Real&     a_startX                                  );
  /// Destructor
  /**
   */
  ~RiemannProblem_MHDK() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Returns physical model used
  /**
   */
  ePhysicalModel physicalModel( void ) const { return PP_MHDPM;}

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real& a_gamma,   const Real& a_startX,
                         const Real& a_netnumL,  
                         const Real& a_numdenL, const Real& a_temperatureL, 
                         const Real& a_velxL,   const Real& a_velyL, 
                         const Real& a_BxL,     const Real& a_ByL,
                         const Real& a_numdenR, const Real& a_temperatureR,
                         const Real& a_velxR,   const Real& a_velyR, 
                         const Real& a_BxR,     const Real& a_ByR    );

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
  
  /// Set up initial conditions
  /**
   */
  void initialize(LevelData<FArrayBox>& a_U);
  
 
protected:
                                // True if the Fortran common block has been set
  bool m_isFortranCommonSet;

  Real  m_numdenL;
  Real  m_temperatureL;
  Real  m_velxL;
  Real  m_velyL;
  Real  m_BxL;
  Real  m_ByL;

  Real  m_numdenR;
  Real  m_temperatureR;
  Real  m_velxR;
  Real  m_velyR;
  Real  m_BxR;
  Real  m_ByR;

  Real  m_startX;
  Real  m_netnumL;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const RiemannProblem_MHDK& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  RiemannProblem_MHDK(const RiemannProblem_MHDK& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
