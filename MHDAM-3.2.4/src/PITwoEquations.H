#ifndef _PI_TWO_EQUATIONS_H_
#define _PI_TWO_EQUATIONS_H_

#include "PickupIons.H"

/// Class for two equation pickup ions model.
/**
 */
class PITwoEquations : public PickupIons
{
public:

  /// Constructor
  /**
   */
  PITwoEquations( void );

  /// Destructor
  /**
   */
  virtual ~PITwoEquations( void );

  /// Get pickup ions model ID
  /**
   */
  virtual PI_MODELS modelID( void ) { return PI_TWO_EQNS;}

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real& a_dSmallVal,
                         const Real& a_dGamma    );

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new PI made with
     new_PickupIons() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PITwoEquations object.
   */
  virtual PickupIons * new_PickupIons( void );

  /// Copy internal data to another PITwoEquations
  /**
   */
  virtual void copyTo( PickupIons * pPIons ) const;

  /// Names of the conservative variables
  /**
     Return the names of the conservative variables.
   */
  virtual void AddConsNames( Vector<string> * pNames );

  /// Names of the primitive variables
  /**
     Return the names of the primitive variables.
   */
  virtual void AddPrimNames( Vector<string> * pNames );

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes( void )
  {
    return 2;
  }

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes( void )
  {
    return 2;
  }

  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void stateToPrim(       FArrayBox & a_W,
                            const FArrayBox & a_U,
                            const Box &       a_box );

  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToState(       FArrayBox & a_U,
                            const FArrayBox & a_W,
                            const Box &       a_box );

  /// Calculate fluxes of pickup ions variables
  /**
   */
  virtual void primToFlux(       FArrayBox & a_F,
                           const FArrayBox & a_W,
                           const int &       a_dir,
                           const Box &       a_box);

  virtual void primToFlux(       FArrayBox & a_F,
                           const FArrayBox & a_W,
                        const BaseFab<int> & a_REG,
                           const int &       a_dir,
                           const Box &       a_box);

  /// Compute an upwind fluxes at the faces for pickup ions variables
  /**
   */
  virtual void upwindFluxes(       FArrayBox & a_F,
                             const FArrayBox & a_WLeft,
                             const FArrayBox & a_WRight,
                             const int &       a_dir,
                             const int &       a_iRho,
                             const Box &       a_box);

  virtual void upwindFluxes(       FArrayBox & a_F,
                             const FArrayBox & a_WLeft,
                             const FArrayBox & a_WRight,
                          const BaseFab<int> & a_REG,
                             const int &       a_dir,
                             const int &       a_iRho,
                             const Box &       a_box);

  /// Calculate pickup ions specific source terms with gradP
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const int       & a_level,
                               const Box       & a_box );
  
  /// Calculate pickup ions specific source terms with divU
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                            const BaseFab<int> & a_REG,
                               const Real      & a_dt,
                               const FArrayBox & a_divU,
                               const int       & a_level,
                               const Box       & a_box );

  /// Compute pickup ion model specific dt in all patch points.
  /**
   */
  virtual Real computeDt( const FArrayBox& a_U,
                                FArrayBox& a_dt,
                          const int      & a_level,
                          const Box      & a_box );

protected:
  double m_dSmallVal;
  double m_dGamma;
                                // True if the Fortran common block has been set
  bool   m_isFortranCommonSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const PITwoEquations& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PITwoEquations(const PITwoEquations& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _PI_TWO_EQUATIONS_H_ */
