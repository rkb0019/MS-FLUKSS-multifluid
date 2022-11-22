#ifndef _PICKUP_IONS_H_
#define _PICKUP_IONS_H_

#include "ParmParse.H"
#include "FArrayBox.H"
#include "Box.H"
#include "CSHandler.H"//FF 

/// Base class for pickup ions.
/**
 */
class PickupIons
{
public:

  enum PI_MODELS {
                   PI_UNDEFINED         = 0,
                   PI_TWO_EQNS          = 1
                 };

  /// Constructor
  /**
   */
  PickupIons( void ){};

  /// Destructor
  /**
   */
  virtual ~PickupIons( void ){};

  /// Get pickup ions model ID
  /**
   */
  virtual PI_MODELS modelID( void ) { return PI_UNDEFINED;}

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Factory method - this object is its own factory
  /**
     Return a point to a new PickupIons object.
   */
  virtual PickupIons * new_PickupIons( void ) = 0;

  /// Copy internal data to another PickupIons
  /**
   */
  virtual void copyTo( PickupIons * pTModel ) const;

  /// Number of conservative variables
  /**
     Return the number of conservative variables.
   */
  int numConservative( void )
  {
    return m_consStateInterval.size();
  }

  /// Interval of conservative variables
  /**
     Return interval of conservative variables.
   */
  const Interval & consInterval( void )
  {
    return m_consStateInterval;
  }

  /// Names of the conservative variables
  /**
     Return the names of the conservative variables.  A default implementation is
     provided that puts in generic names (i.e., "tvariable#" which "#" ranges
     for 0 to numConservative()-1.
   */
  virtual void AddConsNames( Vector<string> * pNames );

  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  int numPrimitives( void )
  {
    return m_primStateInterval.size();
  }

  /// Interval of primitive variables in state array
  /**
     Return interval of primitive variables.
   */
  const Interval & primInterval( void )
  {
    return m_primStateInterval;
  }

  /// Names of the primitive variables
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "tvariable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual void AddPrimNames( Vector<string> * pNames );

  /// Adjust intervals of conservative and primitive variables in state array
  /**
   */
  void adjustIntervals( int iConsFirst, int iPrimFirst );

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes( void ) = 0;

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes( void ) = 0;

  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void stateToPrim(       FArrayBox & a_W,
                            const FArrayBox & a_U,
                            const Box &       a_box ) = 0;

  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToState(       FArrayBox & a_U,
                            const FArrayBox & a_W,
                            const Box &       a_box ) = 0;

  /// Calculate fluxes of pickup ions variables
  /**
   */
  virtual void primToFlux(       FArrayBox & a_F,
                           const FArrayBox & a_W,
                           const int &       a_dir,
                           const Box &       a_box) = 0;

  virtual void primToFlux(       FArrayBox & a_F,
                           const FArrayBox & a_W,
                        const BaseFab<int> & a_REG,
                           const int &       a_dir,
                           const Box &       a_box) = 0;

  /// Compute an upwind fluxes at the faces for pickup ions variables
  /**
   */
  virtual void upwindFluxes(       FArrayBox & a_F,
                             const FArrayBox & a_WLeft,
                             const FArrayBox & a_WRight,
                             const int &       a_dir,
                             const int &       a_iRho,
                             const Box &       a_box) = 0;

  virtual void upwindFluxes(       FArrayBox & a_F,
                             const FArrayBox & a_WLeft,
                             const FArrayBox & a_WRight,
                          const BaseFab<int> & a_REG,
                             const int &       a_dir,
                             const int &       a_iRho,
                             const Box &       a_box) = 0;

  /// Calculate pickup specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const int       & a_level,
                               const Box       & a_box ) = 0;

  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                            const BaseFab<int> & a_REG,
                               const Real      & a_dt,
                               const FArrayBox & a_scale,
                               const int       & a_level,
                               const Box       & a_box ) = 0;

  /// Compute pickup ion model specific dt in all patch points.
  /**
   */
  virtual Real computeDt( const FArrayBox& a_U,
                                FArrayBox& a_dt,
                          const int      & a_level,
                          const Box      & a_box );

/// FF set coordinate system
 virtual void setCoordinateSystem( CoordinateSystemHandler* a_csh )
 {
    m_csh=a_csh;
 };


protected:
  int                 m_verbosity;

  Interval m_consStateInterval;    // interval of conservative variables in an array of state variables
  Interval m_primStateInterval;    // interval of primitive    variables in an array of state variables
// Coordinate system handler
  CoordinateSystemHandler* m_csh;

private:
  // Disallowed for all the usual reasons
  void operator=(const PickupIons& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  PickupIons(const PickupIons& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _PICKUP_IONS_H_ */
