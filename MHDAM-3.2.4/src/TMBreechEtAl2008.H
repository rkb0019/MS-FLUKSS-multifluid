#ifndef _TM_BREECH_2008_H_
#define _TM_BREECH_2008_H_

#include "TurbulenceModel.H"

/// Class for turbulence model by Breech et al 2008.
/**
 */
class TMBreechEtAl2008 : public TurbulenceModel
{
public:

  /// Constructor
  /**
   */
  TMBreechEtAl2008( void );

  /// Destructor
  /**
   */
  virtual ~TMBreechEtAl2008( void );

  /// Get turbulence model ID
  /**
   */
  virtual TM_MODELS modelID( void ) { return TM_BREECH_ET_AL_2008;}

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );

  /// Sets parameters in a common block used by Fortran routines
  /**
   */
  void setFortranCommon( const Real& a_dSmallVal,
                         const Real& a_dAlpha,
                         const Real& a_dBeta,
                         const Real& a_dSigmaD,
                         const Real& a_dfD,
                         const Real& a_dnH,
                         const Real& a_dTion,
                         const Real& a_dLcav,
                         const Real& a_dUr1AU,
                         const Real& a_dVa1AU,
                         const Real& a_dNsw1AU    );

  /// Set the flag m_isFortranCommonSet to true
  /**
     Set the flag m_isFortranCommonSet to true so that new TM made with
     new_TurbulenceModel() will have this flag set without calling setFortranCommon()
     (this is a clumsy design and should be improved).
   */
  void setFortranCommonSet();

  /// Factory method - this object is its own factory
  /**
     Return a point to a new TMBreechEtAl2008 object.
   */
  virtual TurbulenceModel * new_TurbulenceModel( void );

  /// Copy internal data to another TMBreechEtAl2008
  /**
   */
  virtual void copyTo( TurbulenceModel * pTModel ) const;

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
    return 3;
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
    return 3;
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

  /// Calculate fluxes of turbulence model variables
  /**
   */
  virtual void primToFlux(       FArrayBox & a_F,
                           const FArrayBox & a_W,
                           const int &       a_dir,
                           const Box &       a_box);

  /// Compute an upwind fluxes at the faces for turbulence model variables
  /**
   */
  virtual void upwindFluxes(       FArrayBox & a_F,
                             const FArrayBox & a_WLeft,
                             const FArrayBox & a_WRight,
                             const int &       a_dir,
                             const int &       a_iRho,
                             const Box &       a_box);

  /// Calculate turbulence model specific source terms
  /**
   */
  virtual void explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const int       & a_iRhoN1,
                               const int       & a_iRhoPI,
                               const int       & a_level,
                               const Box       & a_box );

  /// Set scales in a common block used by Fortran routines
  /**
   */
  void setScales( const Real& a_dScaleLen,
                  const Real& a_dScaleVel ) const;

protected:
  double m_dSmallVal;
  double m_dAlpha;
  double m_dBeta;
  double m_dSigmaD;
  double m_dfD;
  double m_dnH;
  double m_dTion;
  double m_dLcav;
  double m_dUr1AU;
  double m_dVa1AU;
  double m_dNsw1AU;
                                // True if the Fortran common block has been set
  bool   m_isFortranCommonSet;

private:
  // Disallowed for all the usual reasons
  void operator=(const TMBreechEtAl2008& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  TMBreechEtAl2008(const TMBreechEtAl2008& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _TM_BREECH_2008_H_ */
