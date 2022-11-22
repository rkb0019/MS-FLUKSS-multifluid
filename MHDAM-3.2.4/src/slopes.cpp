#include "parstream.H"

#include "slopes.H"
#include "reconstructionF_F.H"

slopes::slopes( void )
{
  m_limitSlopes          = true;

  m_checkTVD             = false;
  m_checkPositivity      = false;
  m_reduceOrder          = true;

  m_pLimiter             = NULL;
}

slopes::slopes( limiter1D::eLimiters eLim,
                bool                 a_limitSlopes,
                bool                 a_checkPositivity,
                bool                 a_reduceOrder,
                bool                 a_checkTVD )
{
  setSlopeParameters( a_limitSlopes, a_checkPositivity, a_reduceOrder, a_checkTVD );

  m_pLimiter = limiter1D::make( eLim );
}

slopes::~slopes( void )
{
  if( m_pLimiter != NULL )
  {
    delete m_pLimiter;
  }
}

slopes * slopes::new_slopes( void )
{
  slopes * retval = NULL;

  if( m_pLimiter == NULL )
  {
    retval = new slopes;
    retval->setSlopeParameters( m_limitSlopes,
                                m_checkPositivity,
                                m_reduceOrder,
                                m_checkTVD );
  } else {
    retval = new slopes( m_pLimiter->name(),
                         m_limitSlopes,
                         m_checkPositivity,
                         m_reduceOrder,
                         m_checkTVD );
  }

  return retval;
}

void slopes::setSlopeParameters( bool a_limitSlopes,
                                 bool a_checkPositivity,
                                 bool a_reduceOrder,
                                 bool a_checkTVD )
{
  m_limitSlopes          = a_limitSlopes;

  m_checkTVD             = a_checkTVD;
  m_checkPositivity      = a_checkPositivity;
  m_reduceOrder          = a_reduceOrder;
}

void slopes::setLimiter( limiter1D::eLimiters eLim )
{
  if( m_pLimiter != NULL )
  {
    delete m_pLimiter;
  }
  m_pLimiter = limiter1D::make( eLim );
}

bool slopes::limitSlopes( void )
{
  return (m_limitSlopes && (m_pLimiter != NULL));
}

bool slopes::checkPositivity( void )
{
  return m_checkPositivity;
}

bool slopes::orderReduction( void )
{
  return m_reduceOrder;
}

bool slopes::checkTVD( void )
{
  return m_checkTVD;
}

void slopes::setTVDCheck( bool bCheckTVD )
{
  m_checkTVD = bCheckTVD;
}

void slopes::faceValues( const FArrayBox & a_W,
                               FArrayBox & a_WLeft,
                               FArrayBox & a_WRight,
                         const int       & a_var,
                         const int       & a_level,
                         const int       & a_dir,
                         const Box       & a_box,
                         const FArrayBox & a_dx,
                         CoordinateSystemHandler * a_csh )
{
  if( limitSlopes() == true )
  {
    m_pLimiter->faceValues( a_W, a_WLeft, a_WRight, a_var, a_level, a_dir, a_box, a_dx, a_csh );
  } else {
    FORT_CENTEREDFACEVALUES( CHF_CONST_FRA1(a_W,a_var),
                             CHF_FRA1(a_WLeft,a_var),
                             CHF_FRA1(a_WRight,a_var),
                             CHF_INT(a_dir),
                             CHF_BOX(a_box));
  }
}

void slopes::faceSlopes(       FArrayBox         & a_SLeft,
                               FArrayBox         & a_SRight,
                         const int               & a_var,
                         const int               & a_dir,
                         const Box               & a_box,
                         const FArrayBox         & a_dx,
                         CoordinateSystemHandler * a_csh )
{
  if( limitSlopes() == true )
  {
    m_pLimiter->faceSlopes( a_SLeft, a_SRight, a_var, a_dir, a_box, a_dx, a_csh );
  }
}

void slopes::checkPositivity(       BaseFab<int> & a_negative,
                              const FArrayBox    & a_WLeft,
                              const FArrayBox    & a_WRight,
                              const int          & a_var,
                              const Box          & a_box )
{
  if( m_checkPositivity == true )
  {
    FORT_CHECKPOSITIVITY( CHF_FIA1(a_negative,0),
                          CHF_CONST_FRA1(a_WLeft, a_var),
                          CHF_CONST_FRA1(a_WRight,a_var),
                          CHF_BOX(a_box) );
  }
}

void slopes::reduceOrder(       FArrayBox    & a_WLeft,
                                FArrayBox    & a_WRight,
                          const FArrayBox    & a_W,
                          const BaseFab<int> & a_negative,
                          const int          & a_var,
                          const Box          & a_box )
{
  if( m_reduceOrder == true )
  {
    FORT_REDUCEORDER( CHF_FRA1(a_WLeft, a_var),
                      CHF_FRA1(a_WRight,a_var),
                      CHF_CONST_FRA1(a_W,a_var),
                      CHF_CONST_FIA1(a_negative,0),
                      CHF_BOX(a_box) );
  }
}
                                                          // Check TVD condition
void slopes::checkTVDCondition(       FArrayBox    & a_SLeft,
                                      FArrayBox    & a_SRight,
                                const FArrayBox    & a_W,
                                const int          & a_var,
                                const int          & a_dir,
                                const Box          & a_box )
{
  if( m_checkTVD == true )
  {
    FORT_CHECKTVD( CHF_CONST_FRA1(a_W,a_var),
                   CHF_FRA1(a_SLeft,a_var),
                   CHF_FRA1(a_SRight,a_var),
                   CHF_INT(a_dir),
                   CHF_BOX(a_box));
  }
}

