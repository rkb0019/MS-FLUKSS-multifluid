#ifndef  RiemannProblemMF_H_
#define  RiemannProblemMF_H_

#include  <iostream>

#include "LevelData.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "RealVect.H"
#include "AMRIO.H"

#include "MultiFluidProblem.H"
#include "SolidF_F.H"

/// Example of planar Riemann problem
/**
   Parameters:
     a_gamma          - Gamma for polytropic, gamma-law gas
 */
class RiemannProblemMF : public MultiFluidProblem
{
public:
  /// Null constructor
  /**
   */
  RiemannProblemMF();
    
  /// Destructor
  /**
   */
  ~RiemannProblemMF() {};

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Returns physical model used
  /**
   */
  ePhysicalModel physicalModel( void ) const { return PP_3FluidPM;}

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real&     a_gamma,
                         const Real&     a_numdenL,
                         const Real&     a_temperatureL,
                         const Real&     a_velxL,
                         const Real&     a_BxL,
                         const Real&     a_netnumL,
                         const Real&     a_startX,
                         const int&      a_uniform_neutrals                                    );
                         
                         
  void setFortranCommonLR();

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new IBCs made with
     new_PhysProblem() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();
  
  /// Copy method 
  /**
     Copy all data from a_PP to this instance.
   */
  virtual void copy_PhysProblem(const PhysProblem* a_PP);


  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  virtual PhysProblem *new_PhysProblem();

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

  /// Define regions for charge exchange
  /**
   */
  virtual void defineRegions( const FArrayBox    & a_W,
                                    FArrayBox    & a_S,
                                    BaseFab<int> & a_R,
                              const Box          & a_box);
 
protected:
                                // True if the Fortran common block has been set
  bool m_isFortranCommonSet;

  Real  m_densityL;
  Real  m_temperatureL;  
  Real  m_velxL;
  Real  m_velyL;
  Real  m_velzL;
  Real  m_BxL;
  Real  m_ByL;
  Real  m_BzL;
  
  Real  m_densityR;
  Real  m_temperatureR;  
  Real  m_velxR;
  Real  m_velyR;
  Real  m_velzR;
  Real  m_BxR;
  Real  m_ByR;
  Real  m_BzR;
    
  Real  m_startX;    
    
  // m_NInitDistr determines initial distribution of neutrals
  // 0 : | m_netnumL,      m_velxL <SHOCK> m_netnumL    ,  m_velxL |
  // 1 : | m_netnumL,      m_velxL <SHOCK> netFluid2Dens,  0       |
  // 2 : | netFluid2Dens,  0       <SHOCK> netFluid2Dens,  0       | 
  int   m_NInitDistr;
  Real  m_netDenL;
  Real  m_netTemperatureL;
  Real  m_netVelx;
  Real  m_netVely;
  Real  m_netVelz;
  
  Real  m_TMLIM;

  Real  m_gamma;

private:
                                         // Disallowed for all the usual reasons
  void operator=(const RiemannProblemMF& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  RiemannProblemMF(const RiemannProblemMF& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
