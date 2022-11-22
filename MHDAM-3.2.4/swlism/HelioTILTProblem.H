#ifndef  HELIOTILTPROBLEM_H_
#define  HELIOTILTPROBLEM_H_

#include "HeliosphericProblem.H"

/// Base class for 3D heliospheric problem
/**
   Parameters:
 */
class HelioTILTProblem : public HeliosphericProblem
{
public:
  /// Null constructor
  /**
   */
  HelioTILTProblem();

  /// Destructor
  /**
   */
  ~HelioTILTProblem();
  


  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  virtual PhysProblem *new_PhysProblem();
  
  /// Copy method 
  /**
     Copy all data from a_PP to this instance.
   */
  virtual void copy_PhysProblem(const PhysProblem* a_PP);

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
  virtual void fluxBC(       FArrayBox&      a_F,
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
  

  /// Creates tagged cells for dynamic mesh refinement
  /**
    Problem specific cells tagging
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags);

  /// Check geometrical/problem limitations for grid adaptation
  /**
   */  
  virtual void lockedCellsRegrid( BaseFab<int> & a_flag,
                            const FArrayBox&  a_U,
                            const Box&     a_box);

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
  
                                
  
  /// Modifies curvilinear mesh (applies to curvilinear meshes only)
  /**
   */                               
  virtual void defineMesh(const ProblemDomain & a_prob_domain,
                          const Vector<Real>  & a_domainBox);

protected:
                                    

private:
                                         // Disallowed for all the usual reasons
  void operator=(const HelioTILTProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  HelioTILTProblem(const HelioTILTProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /*  HelioTILTProblem_H_ */

