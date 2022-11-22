#ifndef _LIMITER_1D_H_
#define _LIMITER_1D_H_

#include "FArrayBox.H"
#include "Box.H"

class CoordinateSystemHandler;

////////////////////////////////////////////////////////////////////////////////
/// Base class for 1D slope limiting.
/**
 */
  
class limiter1D
{
public:
  enum eLimiters{ eUnlimitedCD         = 0,
                  eMinmod              = 1,
                  eSuperbee            = 2,
                  eVanAlbada           = 3,
                  eMonotonizedCD       = 4,
                  eHarmonic            = 5,
                  eHyperbee            = 6,
                  eKoren               = 7,
                  eVenkatakrishnan     = 8,
                  eVanLeer             = 9,
                  eWENO3               = 10,
                  eWENO3YC             = 11 };


public:
  /// Constructor
  /**
   */
  limiter1D( void );

  /// Destructor
  /**
   */
  virtual ~limiter1D( void );

  /// Static method for construction objects
  /**
   */
  static limiter1D * make( eLimiters lim);

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void ) = 0;

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh ) = 0;

  /// Apply the slope limiter to computed left and right slopes on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh ) = 0;

  virtual eLimiters name( void ) const = 0;
  
private:
                                         // Disallowed for all the usual reasons
  void operator=(const limiter1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  limiter1D(const limiter1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// minmod 1D slope limiter.
/**
 */
class minmodSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  minmodSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~minmodSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eMinmod; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const minmodSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  minmodSL1D(const minmodSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// superbee 1D slope limiter.
/**
 */
class superbeeSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  superbeeSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~superbeeSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eSuperbee; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const superbeeSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  superbeeSL1D(const superbeeSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D slope limiting by van Albada et al.
/**
 */
class vanAlbadaSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  vanAlbadaSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~vanAlbadaSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eVanAlbada; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const vanAlbadaSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  vanAlbadaSL1D(const vanAlbadaSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D motonized central differences
/**
 */
class MCSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  MCSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~MCSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eMonotonizedCD; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const MCSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  MCSL1D(const MCSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D harmonic limiter by van Leer.
/**
 */
class harmonicSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  harmonicSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~harmonicSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eHarmonic; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const harmonicSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  harmonicSL1D(const harmonicSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D HYPERBEE limiter
/**
 */
class hyperbeeSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  hyperbeeSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~hyperbeeSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eHyperbee; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const hyperbeeSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  hyperbeeSL1D(const hyperbeeSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D k-approximation by Koren with k = 1/3.
/**
 */
class KorenSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  KorenSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~KorenSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eKoren; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const KorenSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  KorenSL1D(const KorenSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D slope limiter by Ventakakrishnan.
/**
 */
class VentakakrishnanSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  VentakakrishnanSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~VentakakrishnanSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eVenkatakrishnan; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const VentakakrishnanSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  VentakakrishnanSL1D(const VentakakrishnanSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// 1D k-approximation by van Leer.
/**
 */
class vanLeerSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  vanLeerSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~vanLeerSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eVanLeer; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const vanLeerSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  vanLeerSL1D(const vanLeerSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// Unlimited Central Difference approximation.
/**
 */
class CDSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  CDSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~CDSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eUnlimitedCD; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const CDSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  CDSL1D(const CDSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// WENO3 approximation.
/**
 */
class WENO3SL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  WENO3SL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~WENO3SL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eWENO3; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const WENO3SL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  WENO3SL1D(const WENO3SL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

////////////////////////////////////////////////////////////////////////////////
/// Modified WENO3 approximation by Yamaleev & Carpenter
/**
 */
class WENO3YCSL1D : public limiter1D
{
public:
  /// Constructor
  /**
   */
  WENO3YCSL1D( void ){};

  /// Destructor
  /**
   */
  virtual ~WENO3YCSL1D( void ){};

  /// Factory method - this object is its own factory
  /**
     Return a point to a new limiter1D object.
   */
  virtual limiter1D * new_limiter1D( void );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given the cell-averaged values, a_W, for variable, a_var, apply a slope
     limiter within a_box to generate left and right values on cell faces.
   */
  virtual void faceValues( const FArrayBox &               a_W,
                                 FArrayBox &               a_WLeft,
                                 FArrayBox &               a_WRight,
                           const int &                     a_var,
                           const int &                     a_level,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  /// Apply the slope limiter to computed left and right values on cell faces
  /**
     Given one-side slopes, a_WLeft and a_WRight, for variable, a_var, apply a
     slope limiter within a_box to generate left and right slopes on cell faces
   */
  virtual void faceSlopes(       FArrayBox &               a_SLeft,
                                 FArrayBox &               a_SRight,
                           const int &                     a_var,
                           const int &                     a_dir,
                           const Box &                     a_box,
                           const FArrayBox &               a_dx,
                                 CoordinateSystemHandler * a_csh );

  virtual eLimiters name( void ) const { return eWENO3YC; }

private:
                                         // Disallowed for all the usual reasons
  void operator=(const WENO3YCSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  WENO3YCSL1D(const WENO3YCSL1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /* _LIMITER_1D_H_ */

