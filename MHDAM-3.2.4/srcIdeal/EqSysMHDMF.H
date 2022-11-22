#ifndef _EQSYSMHDMF_H_
#define _EQSYSMHDMF_H_

#include "EquationSystem.H"

///
/**
   The base class EquationSystem descrbes the system of hyperbolic equations 
   solved in the package
 */
  
class EqSysMHDMF : public EquationSystem
{
public:

  /// Constructor
  /**
   */
  EqSysMHDMF(int a_CP, int a_nFluids, int a_ts);

  /// Destructor
  /**
   */
  virtual ~EqSysMHDMF();

  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
   */
  virtual Vector<string> stateNames();
    
  /// Names of the primitive variables
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual Vector<string> primitiveNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes();

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes();


  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void stateToPrim(     FArrayBox& a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box);
                          
  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToState(     FArrayBox& a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box);
  
  /// Returns density index for the fluid 'iFluid' in an array of conservative variables
  /**
   */
  int densityIndexCons(int iFluid);
  
  /// Returns density index for the fluid 'iFluid' in an array of primitive variables
  /**
   */
  int densityIndexPrim(int iFluid);

  /// Transform a_dWLeft and a_dWRight from primitive to characteristic variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the primitive variables.
     On output, it contains the increments in the characteristic variables.
  */
  virtual void charAnalysis(       FArrayBox & a_dWLeft,
                                   FArrayBox & a_dWRight,
                             const FArrayBox & a_W,
                             const int &       a_dir,
                             const Box &       a_box);

  /// Transform a_dWLeft and a_dWRight from characteristic to primitive variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the characteristic variables.
     On output, it contains the increments in the primitive variables.
  */
  virtual void charSynthesis(       FArrayBox & a_dWLeft,
                                    FArrayBox & a_dWRight,
                              const FArrayBox & a_W,
                              const int &       a_dir,
                              const Box &       a_box);
                              
  virtual void vectorVars(Vector<int> & a_vars) const;
  
  /// Returns first components of velocities in the array of conservative/primitive variables
  /**
   */
  virtual void velocityVars(Vector<int> & a_vars) const;

  /// Calculate equation system specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const int       & a_level,
                               const Box       & a_box);

  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                            const BaseFab<int> & a_REG,
                               const FArrayBox & a_divU,
                               const Real      & a_dt,
                               const int       & a_level,
                               const Box       & a_box);

  /// Compute equation system specific dt and returns a cell with the minimum dt.
  /**
   */
  virtual Real computeDt( const FArrayBox& a_U,
                                FArrayBox& a_dt,
                          const Box&       a_box,
                          const int      & a_level,
                                IntVect&   a_minDtCell);

  /// Set the turbulence model object
  /**
   */
  virtual void setTurbulenceModel( TurbulenceModel* a_pTM );

  /// Set the pickup ions object
  /**
   */
  virtual void setPickupIons( PickupIons* a_pPI );

  /// Set usage of the correction potential approach by Dedner et al.
  /**
   */
  void setCorrectionPotentialParams( int a_iUse, Real a_factorCh, Real a_factorCp );
  
  
  /// Set usage of the correction potential approach by Dedner et al.
  /**
   */
  void getCorrectionPotentialParams( Real & a_factorCh, Real & a_factorCp );
  
  /// Set usage of the correction potential approach by Dedner et al.
  /**
   */
  void setCorrectionPotential( Real a_dCh, Real a_dCp );
  
  /// Returns index for correction potential variable, or -1 if the correction potential not used
  /**
   */
  int  correctionPotentialIndex (void);

  /// Is the correction potential approach by Dedner et al used.
  /**
   */
  int isUsedCorrectionPotential( void );


protected:

  int    m_nFluids;
  
  // correction potential data
  int    m_iCorrectionPotential; // correction potential used when positive value
  int    m_UCP;
  int    m_WCP;
  Real   m_dfactorCh;
  Real   m_dfactorCp;
  Real   m_dCh;
  Real   m_dCp;

private:
  
  
  // Disallowed for all the usual reasons
  void operator=(const EqSysMHDMF& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  EqSysMHDMF(const EqSysMHDMF& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
