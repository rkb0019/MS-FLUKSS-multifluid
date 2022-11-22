#ifndef _EQUATIONSYSTEM_H_
#define _EQUATIONSYSTEM_H_

#include "FArrayBox.H"
#include "TurbulenceModel.H"
#include "PickupIons.H"

///
/**
   The base class EquationSystem descrbes the system of hyperbolic equations 
   solved in the package
      
 */
  
class EquationSystem
{
public:  

  /// Constructor
  /**
   */
  EquationSystem();

  /// Destructor
  /**
   */
  virtual ~EquationSystem();  
  
  /// Number of state variables
  /**
     Return the number of state variables.
   */
  virtual int numStates();
  
  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
   */
  virtual Vector<string> stateNames();

  /// Number of conserved variables
  /**
     Return the number of conserved variables.
   */
  inline int numConserved(){return m_consStateInterval.size();}
  
  /// Interval of conserved variables
  /**
     Return the number of conserved variables.
   */
  inline const Interval & consInterval(){return m_consStateInterval;}
  
  /// Interval of level set method variables in state array
  /**
     Interval of level set method variables.
   */
  inline const Interval & lvlsStateInterval(){return m_lvlsStateInterval;}

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  inline int numPrimitives(){return m_nPrim;}
  
  /// Names of the primitive variables
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual Vector<string> primitiveNames();
  
  /// Interval of level set method variables in state array
  /**
     Interval of level set method variables.
   */
  inline const Interval & lvlsPrimInterval(){return m_lvlsPrimInterval;}

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes() = 0;

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes() = 0;

  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void stateToPrim(FArrayBox&      a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box) = 0;
  
  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToState(     FArrayBox& a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box) = 0;


  /// Transform a_dWLeft and a_dWRight from primitive to characteristic variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the primitive variables.
     On output, it contains the increments in the characteristic variables.
  */
  virtual void charAnalysis(       FArrayBox & a_dWLeft,
                                   FArrayBox & a_dWRight,
                             const FArrayBox & a_W,
                             const int &       a_dir,
                             const Box &       a_box) = 0;

  /// Transform a_dWLeft and a_dWRight from characteristic to primitive variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the characteristic variables.
     On output, it contains the increments in the primitive variables.
  */
  virtual void charSynthesis(       FArrayBox & a_dWLeft,
                                    FArrayBox & a_dWRight,
                              const FArrayBox & a_W,
                              const int &       a_dir,
                              const Box &       a_box) = 0;
                              
  /// Number of tracking surfaces by the level set method
  /**
   */
  inline int numTrackingSurfaces() {return m_nTrackingSurfaces;}
      
  
  /// Returns index for a tracking surface 'a_s' in an array of conservative variables
  /**
   */
  inline int lsIndexCons(int a_s) const {CH_assert(a_s<m_nTrackingSurfaces); return m_lvlsStateInterval.begin()+a_s;}
  
  /// Returns index for a tracking surface 'a_s' in an array of primitive variables
  /**
   */
  inline int lsIndexPrim(int a_s) const {CH_assert(a_s<m_nTrackingSurfaces); return m_lvlsPrimInterval.begin()+a_s;}
  
  /// Returns first components of vectors in the array of conservative/primitive variables
  /**
   */
  virtual void vectorVars(Vector<int> & a_vars) const = 0;
  
  /// Returns first components of velocities in the array of conservative/primitive variables
  /**
   */
  virtual void velocityVars(Vector<int> & a_vars) const = 0;

  /// Calculate equation system specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
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

  /// Get the turbulence model object
  /**
   */
  TurbulenceModel* getTurbulenceModel( void ) const;

  /// Get the pickup ions object
  /**
   */
  PickupIons* getPickupIons( void ) const;

  /// Set the pickup ions object
  /**
   */
  virtual void setPickupIons( PickupIons* a_pPI );

protected:

  /// Define the object
  /**
   */
  void define(const Interval & a_consStateInterval,
              int              a_nPrim,
              int              a_nTrackingSurfaces,
              const Interval & a_lvlsStateInterval,
              Interval       & a_lvlsPrimInterval);

  int m_verbosity;
  
  int      m_nTrackingSurfaces;
  Interval m_consStateInterval;    // interval of conservative variables in an array of state variables
  Interval m_lvlsStateInterval;    // interval of level set    variables in an array of state variables
  int      m_nPrim;
  Interval m_lvlsPrimInterval;

  TurbulenceModel * m_pTrModel;
  int               m_isTrModelSet;
  PickupIons      * m_pPIons;
  int               m_isPIonsSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const EquationSystem& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  EquationSystem(const EquationSystem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
