#ifndef _SLOPES_H_
#define _SLOPES_H_

#include "FArrayBox.H"
#include "Box.H"
#include "limiter1D.H"

class CoordinateSystemHandler;

////////////////////////////////////////////////////////////////////////////////
/// Base class for slope computations.
/**
 */
class slopes
{
protected:
  bool        m_limitSlopes;

  bool        m_checkTVD;
  bool        m_checkPositivity;
  bool        m_reduceOrder;

  limiter1D * m_pLimiter;

public:
  /// Constructor
  /**
   */
  slopes( void );

  slopes( limiter1D::eLimiters eLim,
          bool a_limitSlopes = true,
          bool a_checkPositivity = false,
          bool a_reduceOrder = true,
          bool a_checkTVD = false );

  /// Destructor
  /**
   */
  virtual ~slopes( void );

  /// Factory method - this object is its own factory
  /**
     Return a point to a new slopes object.
   */
  virtual slopes * new_slopes( void );

  /// Set parameters for slope computations
  /**
   */
  void setSlopeParameters( bool a_limitSlopes = true,
                           bool a_checkPositivity = false,
                           bool a_reduceOrder = true,
                           bool a_checkTVD = false );

  /// Set limiter for slope limiting
  /**
   */
  void setLimiter( limiter1D::eLimiters eLim );

  /// Slope limiter query
  /**
     Return true if slope limiting is used 
   */
  bool limitSlopes( void );

  /// Check positivity query
  /**
     Return true if positivity check is used
   */
  bool checkPositivity( void );

  /// Check TVD query
  /**
     Return true if TVD check is used
   */
  bool checkTVD( void );

  /// Set TVD check
  /**
     Set true/false for TVD check
   */
  void setTVDCheck( bool bCheckTVD );

  /// Order reduction query
  /**
     Return true if order has to reduce in case of positivity violation
   */
  virtual bool orderReduction( void );

  /// Apply a slope limiter to computed cell face values
  /**
   */
  virtual void faceValues( const FArrayBox               & a_W,
                                 FArrayBox               & a_WLeft,
                                 FArrayBox               & a_WRight,
                           const int                     & a_var,
                           const int                     & a_level,
                           const int                     & a_dir,
                           const Box                     & a_box,
                           const FArrayBox               & a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply a slope limiter to computed cell face slopes
  /**
   */
  virtual void faceSlopes(       FArrayBox               & a_SLeft,
                                 FArrayBox               & a_SRight,
                           const int                     & a_var,
                           const int                     & a_dir,
                           const Box                     & a_box,
                           const FArrayBox               & a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Check positivity
  /**
   */
  virtual void checkPositivity(       BaseFab<int> & a_negative,
                                const FArrayBox    & a_WLeft,
                                const FArrayBox    & a_WRight,
                                const int          & a_var,
                                const Box          & a_box );

  /// Reduce order
  /**
   */
  virtual void reduceOrder(       FArrayBox    & a_WLeft,
                                  FArrayBox    & a_WRight,
                            const FArrayBox    & a_W,
                            const BaseFab<int> & a_negative,
                            const int          & a_var,
                            const Box          & a_box );

  /// Check TVD condition
  /**
   */
  virtual void checkTVDCondition(       FArrayBox    & a_SLeft,
                                        FArrayBox    & a_SRight,
                                  const FArrayBox    & a_W,
                                  const int          & a_var,
                                  const int          & a_dir,
                                  const Box          & a_box );

private:
                                         // Disallowed for all the usual reasons
  void operator=(const slopes& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  slopes(const slopes& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _SLOPES_H_ */

