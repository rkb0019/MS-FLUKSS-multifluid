#include "TMBreechEtAl2008.H"
#include "PatchMHDAMF_F.H"
#include "TModel_F.H"
#include "LGintegrator.H"
#include "parstream.H"

//                                                                   Constructor
TMBreechEtAl2008::TMBreechEtAl2008( void )
{
  m_consStateInterval.define( 0, 2 );
  m_primStateInterval.define( 0, 2 );

  m_iEPIModel  = 1;

  m_dAlpha     = 0.8;
  m_dBeta      = 0.4;
  m_dSmallVal  = 1.0e-8;
  m_dSigmaD    =-1.0/3.0;
  m_dfD        = 0.25;
  m_dnH        = 0.1;
  m_dTion      = 1.0e6;
  m_dLcav      = 8.0*1.5e+13;
  m_dUr1AU     = 750.0e+5;
  m_dVa1AU     = 50.0e+5;
  m_dNsw1AU    = 3.0;

  m_isFortranCommonSet = false;
}

//                                                                    Destructor
TMBreechEtAl2008::~TMBreechEtAl2008( void )
{
}

//                                                              Input parameters
void TMBreechEtAl2008::input( ParmParse & parser, int a_verbosity )
{
  m_verbosity  = a_verbosity;

  m_iEPIModel  = 1;

  m_dAlpha     = 0.8;
  m_dBeta      = 0.4;
  m_dSmallVal  = 1.0e-8;
  m_dSigmaD    =-1.0/3.0;
  m_dfD        = 0.25;
  m_dnH        = 0.1;
  m_dTion      = 1.0e6;
  m_dLcav      = 8.0*1.5e+13;
  m_dUr1AU     = 750.0e+5;
  m_dVa1AU     = 50.0e+5;
  m_dNsw1AU    = 3.0;

  parser.query( "EPIModel",  m_iEPIModel );
  parser.query( "alphaTM",   m_dAlpha    );
  parser.query( "betaTM",    m_dBeta     );
  parser.query( "smallTM",   m_dSmallVal );
  parser.query( "sigmaDTM",  m_dSigmaD   );
  parser.query( "dDEPI",     m_dfD       );
  parser.query( "nHEPI",     m_dnH       );
  parser.query( "tauIonEPI", m_dTion     );
  parser.query( "LcavEPI",   m_dLcav     );
  parser.query( "ur1AUEPI",  m_dUr1AU    );
  parser.query( "va1AYEPI",  m_dVa1AU    );
  parser.query( "nSW1AUEPI", m_dNsw1AU   );

                                                         // Print the parameters
  if( m_verbosity >= 2 )
  {
    pout() << "Turbulence model parametrs input:"     << endl;
    pout() << "EPIModel  = " << m_iEPIModel << endl;
    pout() << "alphaTM   = " << m_dAlpha    << endl;
    pout() << "betaTM    = " << m_dBeta     << endl;
    pout() << "smallTM   = " << m_dSmallVal << endl;
    pout() << "sigmaDTM  = " << m_dSigmaD   << endl;
    pout() << "dDEPI     = " << m_dfD       << endl;
    pout() << "nHEPI     = " << m_dnH       << endl;
    pout() << "tauIonEPI = " << m_dTion     << endl;
    pout() << "LcavEPI   = " << m_dLcav     << endl;
    pout() << "ur1AUEPI  = " << m_dUr1AU    << endl;
    pout() << "va1AYEPI  = " << m_dVa1AU    << endl;
    pout() << "nSW1AUEPI = " << m_dNsw1AU   << endl;
  }

  setFortranCommon( m_dSmallVal,  m_dAlpha,
                    m_dBeta,      m_dSigmaD,
                    m_dfD,        m_dnH,
                    m_dTion,      m_dLcav,
                    m_dUr1AU,     m_dVa1AU,
                    m_dNsw1AU                );
}

//                    Sets parameters in a common block used by Fortran routines
void TMBreechEtAl2008::setFortranCommon( const Real& a_dSmallVal,
                                         const Real& a_dAlpha,
                                         const Real& a_dBeta,
                                         const Real& a_dSigmaD,
                                         const Real& a_dfD,
                                         const Real& a_dnH,
                                         const Real& a_dTion,
                                         const Real& a_dLcav,
                                         const Real& a_dUr1AU,
                                         const Real& a_dVa1AU,
                                         const Real& a_dNsw1AU    )
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_SETCONST_TM( CHF_CONST_REAL( a_dAlpha  ),
                    CHF_CONST_REAL( a_dBeta ),
                    CHF_CONST_REAL( a_dSmallVal ),
                    CHF_CONST_REAL( a_dSigmaD ),
                    CHF_CONST_REAL( a_dfD ),
                    CHF_CONST_REAL( a_dnH ),
                    CHF_CONST_REAL( a_dTion ),
                    CHF_CONST_REAL( a_dLcav ),
                    CHF_CONST_REAL( a_dUr1AU ),
                    CHF_CONST_REAL( a_dVa1AU ),
                    CHF_CONST_REAL( a_dNsw1AU ) );

  m_isFortranCommonSet = true;
}

//                                     Set the flag m_isFortranCommonSet to true
void TMBreechEtAl2008::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

//                               Factory method - this object is its own factory
TurbulenceModel * TMBreechEtAl2008::new_TurbulenceModel( void )
{
  TMBreechEtAl2008 * retval = new TMBreechEtAl2008();

  return static_cast<TurbulenceModel*>(retval);
}

//                                Copy internal data to another TMBreechEtAl2008
void TMBreechEtAl2008::copyTo( TurbulenceModel * pTModel ) const
{
  TurbulenceModel::copyTo( pTModel );

  TMBreechEtAl2008* pTM  = dynamic_cast<TMBreechEtAl2008*>(pTModel);

  pTM->m_iEPIModel = m_iEPIModel;
  pTM->m_dAlpha    = m_dAlpha;
  pTM->m_dBeta     = m_dBeta;
  pTM->m_dSmallVal = m_dSmallVal;
  pTM->m_dSigmaD   = m_dSigmaD;
  pTM->m_dfD       = m_dfD;
  pTM->m_dnH       = m_dnH;
  pTM->m_dTion     = m_dTion;
  pTM->m_dLcav     = m_dLcav;
  pTM->m_dUr1AU    = m_dUr1AU;
  pTM->m_dVa1AU    = m_dVa1AU;
  pTM->m_dNsw1AU   = m_dNsw1AU;
}

//                                 Generate names for the conservative variables
void TMBreechEtAl2008::AddConsNames( Vector<string> * pNames )
{
  pNames->push_back( "rho_Z2" );
  pNames->push_back( "rho_sigma_c" );
  pNames->push_back( "rho_lambda" );
}

//                                    Generate names for the primitive variables
void TMBreechEtAl2008::AddPrimNames( Vector<string> * pNames )
{
  pNames->push_back( "Z2" );
  pNames->push_back( "sigma_c" );
  pNames->push_back( "lambda" );
}

//     Compute the primitive variables from the conserved variables within a_box
void TMBreechEtAl2008::stateToPrim(       FArrayBox & a_W,
                                    const FArrayBox & a_U,
                                    const Box &       a_box )
{
  int iRhoZ2 = m_consStateInterval.begin();
  int iZ2    = m_primStateInterval.begin();
  FORT_CONSTOPRIM_TM( CHF_FRA(a_W),
                      CHF_CONST_FRA(a_U),
                      CHF_CONST_INT(iRhoZ2),
                      CHF_CONST_INT(iZ2),
                      CHF_BOX(a_box));
}

//     Compute the conserved variables from the primitive variables within a_box
void TMBreechEtAl2008::primToState(       FArrayBox & a_U,
                                    const FArrayBox & a_W,
                                    const Box &       a_box )
{
  int iRhoZ2 = m_consStateInterval.begin();
  int iZ2    = m_primStateInterval.begin();
  FORT_PRIMTOCONS_TM( CHF_FRA(a_U),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(iRhoZ2),
                      CHF_CONST_INT(iZ2),
                      CHF_BOX(a_box));
}

//                                Calculate fluxes of turbulence model variables
void TMBreechEtAl2008::primToFlux(       FArrayBox & a_F,
                                   const FArrayBox & a_W,
                                   const int &       a_dir,
                                   const Box &       a_box )
{
  int iRhoZ2 = m_consStateInterval.begin();
  int iZ2    = m_primStateInterval.begin();

  FORT_FLUXE_TM( CHF_FRA(a_F),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_INT(iRhoZ2),
                 CHF_CONST_INT(iZ2),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(a_box) );
}

//          Compute an upwind fluxes at the faces for turbulence model variables
void TMBreechEtAl2008::upwindFluxes(       FArrayBox & a_F,
                                     const FArrayBox & a_WPlus,
                                     const FArrayBox & a_WMinus,
                                     const int &       a_dir,
                                     const int &       a_iRho,
                                     const Box &       a_box )
{
  int iRhoZ2 = m_consStateInterval.begin();
  int iZ2    = m_primStateInterval.begin();

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
  FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

                                   // Riemann solver computes a_F for all edges
                                      // that are not on the physical boundary.
  FORT_UPWINDSCALARFLUXES( CHF_FRA(a_F),
                           CHF_CONST_FRA(shiftWLeft),
                           CHF_CONST_FRA(shiftWRight),
                           CHF_CONST_INT(a_dir),
                           CHF_CONST_INT(a_iRho),
                           CHF_CONST_INT(iRhoZ2),
                           CHF_CONST_INT(iZ2),
                           CHF_BOX(a_box));

                            // Shift the left and right primitive variable boxes
                            // back to their original position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

//                              Calculate turbulence model specific source terms
void TMBreechEtAl2008::explicitSource(       FArrayBox & a_U,
                                             FArrayBox & a_S,
                                       const FArrayBox & a_W,
                                       const Real      & a_dt,
                                       const int       & a_iRhoN1,
                                       const int       & a_iRhoPI,
                                       const int       & a_level,
                                       const Box       & a_box )
{
  int iRhoZ2 = m_consStateInterval.begin();
  int iZ2    = m_primStateInterval.begin();

  FORT_SOURCEBREECHETAL2008( CHF_FRA(a_S),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_INT(iRhoZ2),
                             CHF_CONST_INT(iZ2),
                             CHF_CONST_INT(a_level),
                             CHF_BOX(a_box));

  int iEPI = m_iEPIModel;

  if( a_iRhoN1 < 1 ) iEPI   = 1;
  if( a_iRhoPI < 1 ) iEPI   = 1;

  switch( m_iEPIModel ) {
  case 1 :
    {
      FORT_SOURCEEPIICENBERG( CHF_FRA(a_S),
                              CHF_CONST_FRA(a_W),
                              CHF_CONST_REAL(a_dt),
                              CHF_CONST_INT(iRhoZ2),
                              CHF_CONST_INT(iZ2),
                              CHF_CONST_INT(a_level),
                              CHF_BOX(a_box) );
      break;
    }
    case 2 :
    {                            
      FORT_SOURCEEPIICENBERG_PI( CHF_FRA(a_S),
                                 CHF_CONST_FRA(a_W),
                                 CHF_CONST_REAL(a_dt),
                                 CHF_CONST_INT(a_iRhoN1),
                                 CHF_CONST_INT(iRhoZ2),
                                 CHF_CONST_INT(iZ2),
                                 CHF_CONST_INT(a_iRhoPI),
                                 CHF_BOX(a_box) );
    }
  }

  int iRhoLm = m_consStateInterval.end();

  FORT_ADDSOURCES( CHF_FRA(a_U),
                   CHF_CONST_FRA(a_S),
                   CHF_CONST_INT(iRhoZ2),
                   CHF_CONST_INT(iRhoLm),
                   CHF_BOX(a_box) );

  int iEng = UENG;

  FORT_ADDSOURCES( CHF_FRA(a_U),
                   CHF_CONST_FRA(a_S),
                   CHF_CONST_INT(iEng),
                   CHF_CONST_INT(iEng),
                   CHF_BOX(a_box) );
}

//                         Set scales in a common block used by Fortran routines
void TMBreechEtAl2008::setScales( const Real& a_dScaleLen,
                                  const Real& a_dScaleVel ) const
{
  CH_assert(m_isFortranCommonSet == true);

  CH_assert(a_dScaleLen != 0.0);
  CH_assert(a_dScaleVel != 0.0);

  FORT_SETSCALES_FOR_EPI( CHF_CONST_REAL(a_dScaleLen),
                          CHF_CONST_REAL(a_dScaleVel) );
}
