#include <iomanip>
#include <iostream>
using std::ifstream;
using std::ios;
using std::endl;

#include "LoHiCenter.H"

#include "Reconstruction.H"
#include "CSHandler.H"
#include "EquationSystem.H"
#include "reconstructionF_F.H"

Reconstruction::Reconstruction( void )
{
  m_pSlopes      = NULL;
  m_iNumOfSlopes = 0;
  m_iCharacteristicReconstruction  = 0;
}

Reconstruction::Reconstruction( int iSlopes )
{
  m_iNumOfSlopes = iSlopes;
  m_pSlopes      = new slopes*[iSlopes];

  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    m_pSlopes[i] = NULL;
  }
  m_iCharacteristicReconstruction  = 0;
}

Reconstruction::~Reconstruction( void )
{
  if( m_pSlopes != NULL )
  {
    for( int i = 0; i < m_iNumOfSlopes; i++ )
    {
      slopes * pSl = m_pSlopes[i];
      if( pSl != NULL ) delete pSl;
    }
    delete m_pSlopes;
  }
}

Reconstruction * Reconstruction::new_Reconstruction( void )
{
  Reconstruction * retval = new Reconstruction( m_iNumOfSlopes );

  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      retval->m_pSlopes[i] = pSl->new_slopes();
    }
  }

  retval->m_iCharacteristicReconstruction  = isCharactistic();

  return retval;
}
                            // Number of variables for which slopes are computed
int Reconstruction::numSlopes( void )
{
  return m_iNumOfSlopes;
}
                                     // Set slopes information for some variable
void Reconstruction::setSlopes( int index, slopes * pSl )
{
  if( index >= m_iNumOfSlopes || index < 0 ) return;

  if( m_pSlopes == NULL ) return;

  slopes * pSlope = m_pSlopes[index];
  if( pSlope != NULL ) delete pSlope;

  m_pSlopes[index] = pSl;
}
                                     // Get slopes information for some variable
slopes * Reconstruction::slope( int index )
{
  if( index >= m_iNumOfSlopes || index < 0 ) return NULL;

  if( m_pSlopes == NULL ) return NULL;

  return m_pSlopes[index];
}

                                                     // Set charateristic or not
void Reconstruction::setCharaterictic( int iCharRec )
{
	m_iCharacteristicReconstruction  = (iCharRec == 0) ? 0 : 1;
}
                                     // Is reconstruction charateristic or not ?
int Reconstruction::isCharactistic( void )
{
  return m_iCharacteristicReconstruction;
}
                                            // Fill all slopes with same limiter
void Reconstruction::setLimiterForAll( limiter1D::eLimiters eLim )
{
  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      pSl->setLimiter( eLim );
    } else {
      m_pSlopes[i] = new slopes( eLim );
    }
  }
}

                                   // Fill all slopes with same slope parameters
void Reconstruction::setParametersForAll( bool a_limitSlopes,
                                          bool a_checkPositivity,
                                          bool a_reduceOrder,
                                          bool a_checkTVD )
{
  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      pSl->setSlopeParameters( a_limitSlopes, a_checkPositivity, a_reduceOrder, a_checkTVD );
    }
  }
}
                            // Apply a slope limiter to computed cell face values
void Reconstruction::faceValues( const FArrayBox & a_W,
                                       FArrayBox & a_WMinus,
                                       FArrayBox & a_WPlus,
                                 const int       & a_level,
                                 const int       & a_dir,
                                 const Box       & a_box,
                                       EquationSystem          * a_pEqSys,
                                       CoordinateSystemHandler * a_csh )
{
                                                  // Number of slopes to compute
  int numSlope = numSlopes(); 

  CH_assert( a_W.nComp() >= numSlope );

  FArrayBox dx_v;

  if( a_csh->constStep(a_dir) == false )
  {
    IntVect iv_off(IntVect::Zero);
    iv_off[a_dir]=1;  

    Box b = a_box;
    b.grow(1);
    IntVect lo = b.smallEnd()*iv_off;
    IntVect hi = b.bigEnd()*iv_off;

    dx_v.define(Box(lo,hi),1);

    a_csh->dx(dx_v, dx_v.box(), a_dir, a_level);
  }

  if( isCharactistic() == 0 )
  {
    for( int i = 0; i < numSlope; i++ )
    {
      slopes * pSl = m_pSlopes[i];
      if( pSl != NULL )
      {
        pSl->faceValues( a_W, a_WMinus, a_WPlus, i, a_level, a_dir, a_box, dx_v, a_csh );
      } else {
        FORT_CONSTANTFACEVALUES( CHF_CONST_FRA1(a_W,i),
                                 CHF_FRA1(a_WMinus,i),
                                 CHF_FRA1(a_WPlus,i),
                                 CHF_INT(a_dir),
                                 CHF_BOX(a_box));
      }
    }
  } else {
    int iBGN = 0;
    int iEND = numSlope - 1;

    FORT_INITIALSLOPES( CHF_CONST_FRA(a_W),
                        CHF_FRA(a_WMinus),
                        CHF_FRA(a_WPlus),
                        CHF_INT(a_dir),
                        CHF_INT(iBGN),
                        CHF_INT(iEND),
                        CHF_BOX(a_box));

    a_pEqSys->charAnalysis( a_WMinus, a_WPlus, a_W, a_dir, a_box );

    for( int i = 0; i < numSlope; i++ )
    {
      slopes * pSl = m_pSlopes[i];
      if( pSl != NULL )
      {
        pSl->faceSlopes( a_WMinus, a_WPlus, i, a_dir, a_box, dx_v, a_csh );
      }
    }

    a_pEqSys->charSynthesis( a_WMinus, a_WPlus, a_W, a_dir, a_box );

    checkTVDCondition( a_WMinus, a_WPlus, a_W, a_dir, a_box );

  }
}

                                     // Check positivity of computed face values
void Reconstruction::checkPositivity(       FArrayBox & a_WMinus,
                                            FArrayBox & a_WPlus,
                                      const FArrayBox & a_W,
                                      const Box       & a_box)
{
  BaseFab<int> negative( a_box, 1 );
  negative.setVal( 0 );

  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      pSl->checkPositivity( negative, a_WMinus, a_WPlus, i, a_box );
    }
  }

  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      pSl->reduceOrder( a_WMinus, a_WPlus, a_W, negative, i, a_box );
    }
  }
}
                                 // Check TVD condition for computed face values
void Reconstruction::checkTVDCondition(       FArrayBox & a_WMinus,
                                              FArrayBox & a_WPlus,
                                        const FArrayBox & a_W,
                                        const int       & a_dir,
                                        const Box       & a_box)
{
  for( int i = 0; i < m_iNumOfSlopes; i++ )
  {
    slopes * pSl = m_pSlopes[i];
    if( pSl != NULL )
    {
      pSl->checkTVDCondition( a_WMinus, a_WPlus, a_W, i, a_dir, a_box );
    }
  }
}
