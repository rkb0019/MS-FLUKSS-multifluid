#ifndef _PATCH_EULER_H_
#define _PATCH_EULER_H_

#include "PatchMHDAM.H"
#include "PatchMHDAMF_F.H"

/// An version PatchMHDAM for polytropic, gamma-law gas dynamics.
/**
 */
class PatchEuler : public PatchMHDAM
{
public:
  ///
  /**
   */
  PatchEuler();

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new PatchMHDAM object with its initial and boundary
     condtions, slope parameters, and artificial viscosity information defined.
   */
  virtual PatchMHDAM* new_patchMHDAM() const;

  /// Copy internal data to another PatchMHDAM
  /**
   */
  virtual void copyTo( PatchMHDAM * pPatch ) const;

                               
  /// Compute dt and returns a cell with the minimum dt. 
  /**
   */                               
  virtual Real computeDt( const FArrayBox& a_U,
                          const Box&     a_box,
                          IntVect&       a_minDtCell);    
    
                          
  /// Creates tagged cells for dynamic mesh refinement
  /**
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags);
  
  /// This method returns interval that  includes indices of all variables used in adaptation criteria.
  /**
    Before calling tagCells ghost cells must be filled.    
    But we don't need all variables. Indeed, we need only variables 
    which are used in adaptation criteria. This method returns interval
    that includes all variables used in adaptation criteria.     
   */  
  virtual void tagCellVarsInterval(Interval& a_interval);

  /// Calculate magnetic field in cell centers using electric field on edges
  /**
     Empty function for hydrodynamics
   */
  virtual void recalculateMagneticField(       FArrayBox & a_U,
                                         const FArrayBox & a_Uold,
                                         const EdgeBox   & a_E,
                                         const Real      & a_dt,                                         
                                         const Box       & a_box ){};

  /// Calculate electric field on edges
  /**
     Empty function for hydrodynamics
   */
  virtual void calculateElectricField( const FArrayBox & a_W,
                                       const FArrayBox & a_Wold,
                                       const FArrayBox   a_F[CH_SPACEDIM],
                                             EdgeBox   & a_E,
                                       const Real      & a_dt,                                       
                                       const Box       & a_box ){};
                                       
    
protected:
  
  
                              
  virtual void fluxesHancock(         FArrayBox & a_FMinus,
                                      FArrayBox & a_FPlus,
                                const FArrayBox & a_WMinus,
                                const FArrayBox & a_WPlus,
                                const int &       a_dir,
                                const Box &       a_box);

  virtual void fluxesRP(       FArrayBox & a_F,
                         const FArrayBox & a_WLeft,
                         const FArrayBox & a_WRight,
                         const int &       a_dir,
                         const Box &       a_box);

  virtual void addExplicitSources(       FArrayBox    & a_U,
                                   const FArrayBox    & a_W,
                                   const FArrayBox    & a_SOut,
                                         FArrayBox    & a_S,
                                   const FluxBox      & a_Bn,
                                         FArrayBox    & a_divB,
                                         BaseFab<int> & a_REG,
                                   const Real         & a_dt,
                                   const FArrayBox    & a_scale,  
                                   const Box          & a_box);

  virtual void addExplicitSources(       FArrayBox    & a_U,
                                   const FArrayBox    & a_W,
                                   const FArrayBox    & a_SOut,
                                         FArrayBox    & a_S,
                                   const FluxBox      & a_Bn,
                                         FArrayBox    & a_divB,
                                   const FluxBox      & a_Un,
                                         BaseFab<int> & a_REG,
                                   const Real         & a_dt,
                                   const FArrayBox    & a_scale,  
                                   const Box          & a_box);

  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_Uold,
                               const Real      & a_dt,
                               const Box       & a_box);

  virtual void preprocessing ( const FArrayBox    & a_W,
                                     FArrayBox    & a_S,
                                     BaseFab<int> & a_R,
                               const Box          & a_box);
};

#endif
