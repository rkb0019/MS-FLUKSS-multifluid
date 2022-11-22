#include "limiter1D.H"
#include "reconstructionF_F.H"

#include "CSHandler.H"

/// Constructor
/**
 */
limiter1D::limiter1D( void )
{
  
}

/// Destructor
/**
 */
limiter1D::~limiter1D( void )
{
}

limiter1D * limiter1D::make( eLimiters lim )
{
  limiter1D * pLim = NULL;
  switch( lim ){
  case eUnlimitedCD :
                    pLim = new CDSL1D();
                    break;
  case eMinmod :
                    pLim = new minmodSL1D();
                    break;
  case eSuperbee :
                    pLim = new superbeeSL1D();
                    break;
  case eVanAlbada :
                    pLim = new vanAlbadaSL1D();
                    break;
  case eMonotonizedCD :
                    pLim = new MCSL1D();
                    break;
  case eHarmonic :
                    pLim = new harmonicSL1D();
                    break;
  case eHyperbee :
                    pLim = new hyperbeeSL1D();
                    break;
  case eKoren :
                    pLim = new KorenSL1D();
                    break;
  case eVenkatakrishnan :
                    pLim = new VentakakrishnanSL1D();
                    break;
  case eVanLeer :
                    pLim = new vanLeerSL1D();
                    break;
  case eWENO3 :
                    pLim = new WENO3SL1D();
                    break;
  case eWENO3YC :
                    pLim = new WENO3YCSL1D();
                    break;
  }
    
  return pLim;
}

////////////////////////////////////////////////////////////////////////////////
/// minmod 1D slope limiter.
limiter1D * minmodSL1D::new_limiter1D( void )
{
  minmodSL1D * retval = new minmodSL1D();

  return static_cast<limiter1D*>(retval);
}

void minmodSL1D::faceValues( const FArrayBox & a_W,
                                   FArrayBox & a_WLeft,
                                   FArrayBox & a_WRight,
                             const int &       a_var,
                             const int &       a_level,
                             const int &       a_dir,
                             const Box &       a_box,
                             const FArrayBox & a_dx,
                                   CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
    FORT_MINMODFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                           CHF_FRA1(a_WLeft,a_var),
                           CHF_FRA1(a_WRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
  else
  {        
    FORT_MINMODFACEVALUES_V( CHF_CONST_FRA1(a_W,a_var),
                             CHF_FRA1(a_WLeft,a_var),
                             CHF_FRA1(a_WRight,a_var),
                             CHF_CONST_FRA1(a_dx,0),
                             CHF_INT(a_dir),
                             CHF_BOX(a_box));
  }
}

void minmodSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                   FArrayBox & a_SRight,
                             const int &       a_var,
                             const int &       a_dir,
                             const Box &       a_box,
                             const FArrayBox & a_dx,
                                   CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
    FORT_MINMODFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                           CHF_FRA1(a_SRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
  else
  {        
    FORT_MINMODFACESLOPES_V( CHF_FRA1(a_SLeft,a_var),
                             CHF_FRA1(a_SRight,a_var),
                             CHF_CONST_FRA1(a_dx,0),
                             CHF_INT(a_dir),
                             CHF_BOX(a_box));
  }
}

////////////////////////////////////////////////////////////////////////////////
/// superbee 1D slope limiter.
limiter1D * superbeeSL1D::new_limiter1D( void )
{
  superbeeSL1D * retval = new superbeeSL1D();

  return static_cast<limiter1D*>(retval);
}

void superbeeSL1D::faceValues( const FArrayBox & a_W,
                                     FArrayBox & a_WLeft,
                                     FArrayBox & a_WRight,
                               const int &       a_var,
                               const int &       a_level,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                                     CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
  FORT_SUPERBEEFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                           CHF_FRA1(a_WLeft,a_var),
                           CHF_FRA1(a_WRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
  else
  FORT_SUPERBEEFACEVALUES_V( CHF_CONST_FRA1(a_W,a_var),
                           CHF_FRA1(a_WLeft,a_var),
                           CHF_FRA1(a_WRight,a_var),
                           CHF_CONST_FRA1(a_dx,0),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
  
}

void superbeeSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                     FArrayBox & a_SRight,
                               const int &       a_var,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                                     CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
  FORT_SUPERBEEFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                           CHF_FRA1(a_SRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
  else
  FORT_SUPERBEEFACESLOPES_V( CHF_FRA1(a_SLeft,a_var),
                           CHF_FRA1(a_SRight,a_var),
                           CHF_CONST_FRA1(a_dx,0),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D slope limiting by van Albada et al.
limiter1D * vanAlbadaSL1D::new_limiter1D( void )
{
  vanAlbadaSL1D * retval = new vanAlbadaSL1D();

  return static_cast<limiter1D*>(retval);
}

void vanAlbadaSL1D::faceValues( const FArrayBox & a_W,
                                      FArrayBox & a_WLeft,
                                      FArrayBox & a_WRight,
                                const int &       a_var,
                                const int &       a_level,
                                const int &       a_dir,
                                const Box &       a_box,
                                const FArrayBox & a_dx,
                                      CoordinateSystemHandler * a_csh )
{
  FORT_VANALBADAFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                            CHF_FRA1(a_WLeft,a_var),
                            CHF_FRA1(a_WRight,a_var),
                            CHF_INT(a_dir),
                            CHF_BOX(a_box));
}

void vanAlbadaSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                      FArrayBox & a_SRight,
                                const int &       a_var,
                                const int &       a_dir,
                                const Box &       a_box,
                                const FArrayBox & a_dx,
                                      CoordinateSystemHandler * a_csh )
{
  FORT_VANALBADAFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                            CHF_FRA1(a_SRight,a_var),
                            CHF_INT(a_dir),
                            CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D motonized central differences
limiter1D * MCSL1D::new_limiter1D( void )
{
  MCSL1D * retval = new MCSL1D();

  return static_cast<limiter1D*>(retval);
}

void MCSL1D::faceValues( const FArrayBox & a_W,
                               FArrayBox & a_WLeft,
                               FArrayBox & a_WRight,
                         const int &       a_var,
                         const int &       a_level,
                         const int &       a_dir,
                         const Box &       a_box,
                         const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
    FORT_MCFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                       CHF_FRA1(a_WLeft,a_var),
                       CHF_FRA1(a_WRight,a_var),
                       CHF_INT(a_dir),
                       CHF_BOX(a_box));
  else                   
    FORT_MCFACEVALUES_V( CHF_CONST_FRA1(a_W,a_var),
                         CHF_FRA1(a_WLeft,a_var),
                         CHF_FRA1(a_WRight,a_var),
                         CHF_CONST_FRA1(a_dx,0),
                         CHF_INT(a_dir),
                         CHF_BOX(a_box));
}

void MCSL1D::faceSlopes(       FArrayBox & a_SLeft,
                               FArrayBox & a_SRight,
                         const int &       a_var,
                         const int &       a_dir,
                         const Box &       a_box,
                         const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  if( a_csh->constStep( a_dir ) == true )
    FORT_MCFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                       CHF_FRA1(a_SRight,a_var),
                       CHF_INT(a_dir),
                       CHF_BOX(a_box));
  else                   
    FORT_MCFACESLOPES_V( CHF_FRA1(a_SLeft,a_var),
                         CHF_FRA1(a_SRight,a_var),
                         CHF_CONST_FRA1(a_dx,0),
                         CHF_INT(a_dir),
                         CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D harmonic limiter by van Leer.
limiter1D * harmonicSL1D::new_limiter1D( void )
{
  harmonicSL1D * retval = new harmonicSL1D();

  return static_cast<limiter1D*>(retval);
}

void harmonicSL1D::faceValues( const FArrayBox & a_W,
                                     FArrayBox & a_WLeft,
                                     FArrayBox & a_WRight,
                               const int &       a_var,
                               const int &       a_level,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_HARMONICFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                           CHF_FRA1(a_WLeft,a_var),
                           CHF_FRA1(a_WRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
}

void harmonicSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                     FArrayBox & a_SRight,
                               const int &       a_var,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_HARMONICFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                           CHF_FRA1(a_SRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D HYPERBEE limiter
limiter1D * hyperbeeSL1D::new_limiter1D( void )
{
  hyperbeeSL1D * retval = new hyperbeeSL1D();

  return static_cast<limiter1D*>(retval);
}

void hyperbeeSL1D::faceValues( const FArrayBox & a_W,
                                     FArrayBox & a_WLeft,
                                     FArrayBox & a_WRight,
                               const int &       a_var,
                               const int &       a_level,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_HYPERBEEFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                           CHF_FRA1(a_WLeft,a_var),
                           CHF_FRA1(a_WRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
}

void hyperbeeSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                     FArrayBox & a_SRight,
                               const int &       a_var,
                               const int &       a_dir,
                               const Box &       a_box,
                               const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_HYPERBEEFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                           CHF_FRA1(a_SRight,a_var),
                           CHF_INT(a_dir),
                           CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D k-approximation by Koren with k = 1/3.
limiter1D * KorenSL1D::new_limiter1D( void )
{
  KorenSL1D * retval = new KorenSL1D();

  return static_cast<limiter1D*>(retval);
}

void KorenSL1D::faceValues( const FArrayBox & a_W,
                                  FArrayBox & a_WLeft,
                                  FArrayBox & a_WRight,
                            const int &       a_var,
                            const int &       a_level,
                            const int &       a_dir,
                            const Box &       a_box,
                            const FArrayBox & a_dx,
                            CoordinateSystemHandler * a_csh )
{
  FORT_KORENFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                        CHF_FRA1(a_WLeft,a_var),
                        CHF_FRA1(a_WRight,a_var),
                        CHF_INT(a_dir),
                        CHF_BOX(a_box));
}

void KorenSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                  FArrayBox & a_SRight,
                            const int &       a_var,
                            const int &       a_dir,
                            const Box &       a_box,
                            const FArrayBox & a_dx,
                            CoordinateSystemHandler * a_csh )
{
  FORT_KORENFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                        CHF_FRA1(a_SRight,a_var),
                        CHF_INT(a_dir),
                        CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D slope limiter by Ventakakrishnan.
limiter1D * VentakakrishnanSL1D::new_limiter1D( void )
{
  VentakakrishnanSL1D * retval = new VentakakrishnanSL1D();

  return static_cast<limiter1D*>(retval);
}

void VentakakrishnanSL1D::faceValues( const FArrayBox & a_W,
                                            FArrayBox & a_WLeft,
                                            FArrayBox & a_WRight,
                                      const int &       a_var,
                                      const int &       a_level,
                                      const int &       a_dir,
                                      const Box &       a_box,
                                      const FArrayBox & a_dx,
                                      CoordinateSystemHandler * a_csh )
{
  FORT_VENKATAKRISHNANFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                                  CHF_FRA1(a_WLeft,a_var),
                                  CHF_FRA1(a_WRight,a_var),
                                  CHF_INT(a_dir),
                                  CHF_BOX(a_box));
}

void VentakakrishnanSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                            FArrayBox & a_SRight,
                                      const int &       a_var,
                                      const int &       a_dir,
                                      const Box &       a_box,
                                      const FArrayBox & a_dx,
                                      CoordinateSystemHandler * a_csh )
{
  FORT_VENKATAKRISHNANFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                                  CHF_FRA1(a_SRight,a_var),
                                  CHF_INT(a_dir),
                                  CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// 1D slope limiting by van Leer.
limiter1D * vanLeerSL1D::new_limiter1D( void )
{
  vanLeerSL1D * retval = new vanLeerSL1D();

  return static_cast<limiter1D*>(retval);
}

void vanLeerSL1D::faceValues( const FArrayBox & a_W,
                                    FArrayBox & a_WLeft,
                                    FArrayBox & a_WRight,
                              const int &       a_var,
                              const int &       a_level,
                              const int &       a_dir,
                              const Box &       a_box,
                              const FArrayBox & a_dx,
                              CoordinateSystemHandler * a_csh )
{
  FORT_VANLEERFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                          CHF_FRA1(a_WLeft,a_var),
                          CHF_FRA1(a_WRight,a_var),
                          CHF_INT(a_dir),
                          CHF_BOX(a_box));
}

void vanLeerSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                    FArrayBox & a_SRight,
                              const int &       a_var,
                              const int &       a_dir,
                              const Box &       a_box,
                              const FArrayBox & a_dx,
                              CoordinateSystemHandler * a_csh )
{
  FORT_VANLEERFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                          CHF_FRA1(a_SRight,a_var),
                          CHF_INT(a_dir),
                          CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// Unlimited central approximation.
limiter1D * CDSL1D::new_limiter1D( void )
{
  CDSL1D * retval = new CDSL1D();

  return static_cast<limiter1D*>(retval);
}

void CDSL1D::faceValues( const FArrayBox & a_W,
                               FArrayBox & a_WLeft,
                               FArrayBox & a_WRight,
                         const int &       a_var,
                         const int &       a_level,
                         const int &       a_dir,
                         const Box &       a_box,
                         const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_CDFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                     CHF_FRA1(a_WLeft,a_var),
                     CHF_FRA1(a_WRight,a_var),
                     CHF_INT(a_dir),
                     CHF_BOX(a_box));
}

void CDSL1D::faceSlopes(       FArrayBox & a_SLeft,
                               FArrayBox & a_SRight,
                         const int &       a_var,
                         const int &       a_dir,
                         const Box &       a_box,
                         const FArrayBox & a_dx,
                               CoordinateSystemHandler * a_csh )
{
  FORT_CDFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                     CHF_FRA1(a_SRight,a_var),
                     CHF_INT(a_dir),
                     CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// WENO3 approximation.
limiter1D * WENO3SL1D::new_limiter1D( void )
{
  WENO3SL1D * retval = new WENO3SL1D();

  return static_cast<limiter1D*>(retval);
}

void WENO3SL1D::faceValues( const FArrayBox & a_W,
                                  FArrayBox & a_WLeft,
                                  FArrayBox & a_WRight,
                            const int &       a_var,
                            const int &       a_level,
                            const int &       a_dir,
                            const Box &       a_box,
                            const FArrayBox & a_dx,
                            CoordinateSystemHandler * a_csh )
{
  FORT_WENO3FACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                        CHF_FRA1(a_WLeft,a_var),
                        CHF_FRA1(a_WRight,a_var),
                        CHF_INT(a_dir),
                        CHF_BOX(a_box));
}

void WENO3SL1D::faceSlopes(       FArrayBox & a_SLeft,
                                  FArrayBox & a_SRight,
                            const int &       a_var,
                            const int &       a_dir,
                            const Box &       a_box,
                            const FArrayBox & a_dx,
                            CoordinateSystemHandler * a_csh )
{
  FORT_WENO3FACESLOPES( CHF_FRA1(a_SLeft,a_var),
                        CHF_FRA1(a_SRight,a_var),
                        CHF_INT(a_dir),
                        CHF_BOX(a_box));
}

////////////////////////////////////////////////////////////////////////////////
/// Modified WENO3 approximation by Yamaleev & Carpenter
limiter1D * WENO3YCSL1D::new_limiter1D( void )
{
  WENO3YCSL1D * retval = new WENO3YCSL1D();

  return static_cast<limiter1D*>(retval);
}

void WENO3YCSL1D::faceValues( const FArrayBox & a_W,
                                    FArrayBox & a_WLeft,
                                    FArrayBox & a_WRight,
                              const int &       a_var,
                              const int &       a_level,
                              const int &       a_dir,
                              const Box &       a_box,
                              const FArrayBox & a_dx,
                              CoordinateSystemHandler * a_csh )
{
  FORT_WENO3YCFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                          CHF_FRA1(a_WLeft,a_var),
                          CHF_FRA1(a_WRight,a_var),
                          CHF_INT(a_dir),
                          CHF_BOX(a_box));
}

void WENO3YCSL1D::faceSlopes(       FArrayBox & a_SLeft,
                                    FArrayBox & a_SRight,
                              const int &       a_var,
                              const int &       a_dir,
                              const Box &       a_box,
                              const FArrayBox & a_dx,
                              CoordinateSystemHandler * a_csh )
{
  FORT_WENO3YCFACESLOPES( CHF_FRA1(a_SLeft,a_var),
                          CHF_FRA1(a_SRight,a_var),
                          CHF_INT(a_dir),
                          CHF_BOX(a_box));
}
