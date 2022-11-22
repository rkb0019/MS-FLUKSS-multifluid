#ifndef _RECONSTRUCTION_H_
#define _RECONSTRUCTION_H_

#include "ProblemDomain.H"
#include "slopes.H"

/// Control class for reconstruction.
/**
 */
 
class CoordinateSystemHandler;
class EquationSystem;

class Reconstruction
{
protected:

  slopes ** m_pSlopes;
  int       m_iNumOfSlopes;
  int       m_iCharacteristicReconstruction;

public:
  /// Constructor
  /**
   */
  Reconstruction( void );

  Reconstruction( int iSlopes );

  /// Destructor
  /**
   */
  virtual ~Reconstruction( void );

  /// Factory method - this object is its own factory
  /**
     Return a point to a new Reconstruction object.
   */
  virtual Reconstruction * new_Reconstruction( void );

  /// Number of variables for which slopes are computed
  /**
   */
  int numSlopes( void );

  /// Set slopes information for some variable
  /**
   */
  void setSlopes( int index, slopes * pSl );

  /// Get slopes information for some variable
  /**
   */
  slopes * slope( int index );

  /// Set charateristic or not
  /**
   */
  void setCharaterictic( int iCharRec );

  /// Is reconstruction charateristic or not ?
  /**
   */
  int isCharactistic( void );

  /// Fill all slopes with same limiter
  /**
   */
  virtual void setLimiterForAll( limiter1D::eLimiters eLim);

  /// Fill all slopes with same slope parameters
  /**
   */
  virtual void setParametersForAll( bool a_limitSlopes = true,
                                    bool a_checkPositivity = false,
                                    bool a_reduceOrder = true,
                                    bool a_checkTVD = false );

  /// Apply a slope limiter to computed cell face values
  /**
   */
  virtual void faceValues( const FArrayBox & a_W,
                                 FArrayBox & a_WMinus,
                                 FArrayBox & a_WPlus,
                           const int       & a_level,
                           const int       & a_dir,
                           const Box       & a_box,
                                 EquationSystem          * a_pEqSys,
                                 CoordinateSystemHandler * a_csh);

  /// Check positivity of computed face values
  /**
   */
  virtual void checkPositivity(       FArrayBox & a_WMinus,
                                      FArrayBox & a_WPlus,
                                const FArrayBox & a_W,
                                const Box       & a_box);

  /// Check TVD condition for computed face values
  /**
   */
  virtual void checkTVDCondition(       FArrayBox & a_WMinus,
                                        FArrayBox & a_WPlus,
                                  const FArrayBox & a_W,
                                  const int       & a_dir,
                                  const Box       & a_box);

private:
  // Disallowed for all the usual reasons
  void operator=(const Reconstruction& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  Reconstruction(const Reconstruction& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _RECONSTRUCTION_H_ */

