#ifndef _PATCH_MHDMF_H_
#define _PATCH_MHDMF_H_

#include "PatchMHDAM.H"

/// Patch for MHD and multifluid calculations
/**
 */
class PatchMHDMF : public PatchMHDAM
{
public:
  ///
  /**
   */
  PatchMHDMF(int a_nFluids);

  /// Factory method - this object is its own factory
  /**
     Return a pointer to new PatchMHDAM object with its initial and boundary
     condtions, slope parameters, and artificial viscosity information defined.
   */
  virtual PatchMHDAM* new_patchMHDAM() const;
  
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

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
  
  /// This method returns interval that includes indices of all variables used in adaptation criteria.
  /**
    Before calling tagCells ghost cells must be filled.    
    But we don't need all variables. Indeed, we need only variables 
    which are used in adaptation criteria. This method returns interval
    that includes all variables used in adaptation criteria.     
   */  
  virtual void tagCellVarsInterval(Interval& a_interval);

  /// Things to do after advancing this level by one time step.
  /**
   */
  virtual void postTimeStep( void );
  
  /// Calculate magnetic field in cell centers using electric field on edges
  /**
   */
  virtual void recalculateMagneticField(       FArrayBox & a_U,
                                         const FArrayBox & a_Uold,
                                         const EdgeBox   & a_E,
                                         const Real      & a_dt,                                         
                                         const Box       & a_box);

  /// Calculate electric field on edges
  /**
   */
  virtual void calculateElectricField( const FArrayBox & a_W,
                                       const FArrayBox & a_Wold,
                                       const FArrayBox   a_F[CH_SPACEDIM],
                                             EdgeBox   & a_E,
                                       const Real      & a_dt,                                       
                                       const Box       & a_box );
                                       
  /// Calculate divB from conservtive variables
  /**
   */                                  
  virtual void computeDivB(          FArrayBox& a_divB,                              
                              const  FArrayBox& a_U,                                                                 
                              const  Box&       a_box);
                                       
                                    


protected:

  /// Number of fluids
  /**
   */
  int m_nFluids;    
        
                              
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
                                   const FArrayBox    & a_scale,  // dt/volume
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
                                   const FArrayBox    & a_scale,  // dt/volume
                                   const Box          & a_box);

  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_Uold,
                               const Real      & a_dt,
                               const Box       & a_box);

  virtual void preprocessing ( const FArrayBox    & a_W,
                                     FArrayBox    & a_S,
                                     BaseFab<int> & a_R,
                               const Box          & a_box);
                                                          
  virtual void computeDivB(     FArrayBox & a_divB,
                          const FluxBox   & a_Bn,                                                    
                          const FArrayBox & a_W,
                          const Box       & a_box);                                     

                             
  virtual void computeBn(       FArrayBox & a_Bn,
                          const FArrayBox & a_W,
                          const FArrayBox & a_WMinus,
                          const FArrayBox & a_WPlus ,
                          const int       & a_method,
                          const int       & a_dir,
                          const Box       & a_box);       

  virtual void computeUn(       FArrayBox & a_Un,
                          const FArrayBox & a_W,
                          const FArrayBox & a_WMinus,
                          const FArrayBox & a_WPlus ,
                          const int       & a_method,
                          const int       & a_dir,
                          const Box       & a_box);       

  virtual void computeDivU(     FArrayBox & a_divU,
                          const FluxBox   & a_Un,                                                    
                          const FArrayBox & a_W,
                          const Box       & a_box);                                     
                        
  virtual void correctBn(       FArrayBox & a_Bn,
                          const FArrayBox & a_phi,                          
                          const int       & a_dir,
                          const Box       & a_box);
                      
};

#endif
