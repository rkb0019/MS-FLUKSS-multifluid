#ifndef  FASTSHOCK_H_
#define  FASTSHOCK_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "PhysProblem.H"

/// Instability of fast MHD shock
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class FastShock : public PhysProblem
{
public:
  /// Null constructor
  /**
   */
  FastShock();
  
  /// Destructor
  /**
   */
  ~FastShock() {};

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
  void setFortranCommon( const Real&     a_gamma,
                         const Real&     a_VelModuleL,
                         const Real&     a_VelAngleL,
                         const Real&     a_MagModuleL,
                         const Real&     a_MagAngleL,
                         const Real&     a_PressRatioL,
                         const Real&     a_FastMachL,
                         const Real&     a_startX,
                         const Real&     a_AlfvenMag,
                         const Real&     a_AlfvenK       );

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

  Real  m_VelModuleL;
  Real  m_VelAngleL;
  Real  m_MagModuleL;
  Real  m_MagAngleL;
  Real  m_PressRatioL;
  Real  m_FastMachL;
  Real  m_startX;
  Real  m_AlfvenMag;
  Real  m_AlfvenK;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const FastShock& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  FastShock(const FastShock& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
