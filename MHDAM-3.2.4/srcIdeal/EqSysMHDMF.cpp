#include "EqSysMHDMF.H"

#include "PatchMHDMFF_F.H"
#include "PatchMHDAMF_F.H"
#include "PatchIdealMHDF_F.H"
#include "DednerF_F.H"
#include "LGintegrator.H"
#include "reconstructionF_F.H"


EqSysMHDMF::EqSysMHDMF(int a_CP, int a_nFluids, int a_ts)
: EquationSystem()
{
  m_iCorrectionPotential = a_CP;
  m_dfactorCh = 1.0;
  m_dfactorCp = 0.5;
  m_dCh       = 1.0;
  m_dCp       = 1.0;

  m_nFluids   = a_nFluids;  
  
  int nCP = 0; 
  m_UCP   = -1;  
  m_WCP   = -1;
  if (m_iCorrectionPotential > 0)
  {
    nCP = 1;
    m_UCP = UBZ+1;
    m_WCP = WBZ+1;
  }
  
  int nPrim   = WNUM + nCP + (m_nFluids-1)*WNUM_E + a_ts; 

  Interval consStateInterval(0,UNUM + nCP + (m_nFluids-1)*UNUM_E - 1);
  
  Interval lvlsStateInterval, lvlsPrimInterval;
  if (a_ts > 0)
  {
    lvlsStateInterval.define(consStateInterval.end()+1, consStateInterval.end()+a_ts);
    lvlsPrimInterval = lvlsStateInterval;
  }

  define(consStateInterval,    // conservative variables interval,
         nPrim,                // number of primitive variables,      
         a_ts,
         lvlsStateInterval,                
         lvlsPrimInterval);             
}

EqSysMHDMF::~EqSysMHDMF()
{
}

                                             // Names of the conserved variables
Vector<string> EqSysMHDMF::stateNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-momentum");
  retval.push_back("y-momentum");
  retval.push_back("z-momentum");
  retval.push_back("energy-density");
  retval.push_back("x-magnetic field");
  retval.push_back("y-magnetic field");
  retval.push_back("z-magnetic field");

  if( m_iCorrectionPotential != 0 )
  {
    retval.push_back("psi");
  }

  TurbulenceModel * pTM = getTurbulenceModel();

  if( pTM != NULL )
  {
    pTM->AddConsNames( &retval );
  }

  PickupIons * pPI = getPickupIons();

  if( pPI != NULL )
  {
    pPI->AddConsNames( &retval );
  }

  char buffer[40];int i;

  for (i=1; i<m_nFluids; i++)
  { 
    sprintf(buffer,"%s%i","density",       i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","x-momentum",    i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","y-momentum",    i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","z-momentum",    i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","energy-density",i);retval.push_back(buffer);
  }
  
  for (i=0; i<m_nTrackingSurfaces; i++)
  { 
    sprintf(buffer,"%s%i","surface",i);retval.push_back(buffer);    
  }

  return retval;
}
                                                     // Number of flux variables
int EqSysMHDMF::numFluxes()
{
  return numStates();
}

Vector<string> EqSysMHDMF::primitiveNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-vel");
  retval.push_back("y-vel");
  retval.push_back("z-vel");
  retval.push_back("p");
  retval.push_back("x-B");
  retval.push_back("y-B");
  retval.push_back("z-B");

  if( m_iCorrectionPotential != 0 )
  {
    retval.push_back( "psi" );
  }

  TurbulenceModel * pTM = getTurbulenceModel();

  if( pTM != NULL )
  {
    pTM->AddPrimNames( &retval );
  }

  PickupIons * pPI = getPickupIons();

  if( pPI != NULL )
  {
    pPI->AddPrimNames( &retval );
  }

  char buffer[40];int i;
      
  for (i=1; i<m_nFluids; i++)
  { 
    sprintf(buffer,"%s%i","density", i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","x-vel",   i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","y-vel",   i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","z-vel",   i);retval.push_back(buffer);
    sprintf(buffer,"%s%i","p",       i);retval.push_back(buffer);
  }    
  
  for (i=0; i<m_nTrackingSurfaces; i++)
  { 
    sprintf(buffer,"%s%i","surface",i);retval.push_back(buffer);    
  }

  return retval;
}

                 //  Number of primitive variables for which slopes are computed
int EqSysMHDMF::numSlopes()
{  
  return m_nPrim;
}

                 // Compute the primitive variables from the conserved variables
void EqSysMHDMF::stateToPrim(       FArrayBox& a_W,
                              const FArrayBox& a_U,
                              const Box&       a_box )
{  
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
    
  Box b     = a_box;
      
  if (m_nFluids==1)   
     FORT_CONSTOPRIM(CHF_FRA(a_W),
                     CHF_CONST_FRA(a_U),
                     CHF_BOX(b));
  else
  {
    int iRhoN = densityIndexCons(1);
    FORT_CONSTOPRIM_MF(CHF_FRA(a_W),
                     CHF_CONST_FRA(a_U),
                     CHF_CONST_INT(iRhoN),     
                     CHF_CONST_INT(m_nFluids),
                     CHF_BOX(b));
  }

  if( m_iCorrectionPotential > 0 )
  {
    a_W.copy(a_U, m_UCP, m_WCP);
  }

  if( m_isTrModelSet == 1 )
  {
    m_pTrModel->stateToPrim( a_W, a_U, a_box );
  }

  if( m_isPIonsSet == 1 )
  {
    m_pPIons->stateToPrim( a_W, a_U, a_box );
  }

  if (m_nTrackingSurfaces>0)
  {
    a_W.copy(a_U, b, m_lvlsStateInterval.begin(), b, m_lvlsPrimInterval.begin(), m_nTrackingSurfaces);
  }
      
}
                 // Compute the conserved variables from the primitive variables
void EqSysMHDMF::primToState(       FArrayBox& a_U,
                              const FArrayBox& a_W,
                              const Box&       a_box )
{  
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));
  
  if (m_nFluids==1)  
    FORT_PRIMTOCONS(CHF_FRA(a_U),
                     CHF_CONST_FRA(a_W),
                     CHF_BOX(a_box));
  else
  {
    int iRhoN = densityIndexCons(1);
    FORT_PRIMTOCONS_MF(CHF_FRA(a_U),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(iRhoN),     
                     CHF_CONST_INT(m_nFluids),
                     CHF_BOX(a_box));
  }

  if( m_iCorrectionPotential > 0 )
  {
    a_U.copy(a_W, m_WCP, m_UCP);
  }

  if( m_isTrModelSet == 1 )
  {
    m_pTrModel->primToState( a_U, a_W, a_box );
  }

  if( m_isPIonsSet == 1 )
  {
    m_pPIons->primToState( a_U, a_W, a_box );
  }

  if (m_nTrackingSurfaces>0)
  {
    a_U.copy(a_W, a_box, m_lvlsPrimInterval.begin(), a_box, m_lvlsStateInterval.begin(), m_nTrackingSurfaces);
  }
}

// Returns density index for the fluid 'iFluid' in an array of conservative variables
int EqSysMHDMF::densityIndexCons(int iFluid)
{
  if( iFluid == 0 ) return 0;
  
  if (iFluid >= m_nFluids) return -1;

  int index  = UNUM + (iFluid-1)*UNUM_E;

  if( m_iCorrectionPotential > 0 )
  {
    index++;
  }

  if( m_isTrModelSet == 1 )
  {
    index += m_pTrModel->numConservative();
  }

  if( m_isPIonsSet == 1 )
  {
    index += m_pPIons->numConservative();
  }

  return index;
}

// Returns density index for the fluid 'iFluid' in an array of primitive variables
int EqSysMHDMF::densityIndexPrim(int iFluid)
{
  if( iFluid == 0 ) return 0;

  int index  = WNUM + (iFluid-1)*WNUM_E;

  if( m_iCorrectionPotential > 0 )
  {
    index++;
  }

  if( m_isTrModelSet == 1 )
  {
    index += m_pTrModel->numPrimitives();
  }

  if( m_isPIonsSet == 1 )
  {
    index += m_pPIons->numPrimitives();
  }

  return index;
}

  // Transform a_dWLeft and a_dWRight from primitive to characteristic variables
void EqSysMHDMF::charAnalysis(       FArrayBox & a_dWLeft,
                                     FArrayBox & a_dWRight,
                               const FArrayBox & a_W,
                               const int &       a_dir,
                               const Box &       a_box)
{
  int iRho = WRHO;
  FORT_CHARANALYSISF( CHF_FRA(a_dWLeft),
                      CHF_FRA(a_dWRight),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRho),
                      CHF_BOX(a_box));
}

  // Transform a_dWLeft and a_dWRight from characteristic to primitive variables
void EqSysMHDMF::charSynthesis(       FArrayBox & a_dWLeft,
                                      FArrayBox & a_dWRight,
                                const FArrayBox & a_W,
                                const int &       a_dir,
                                const Box &       a_box)
{
  int iRho = WRHO;
  FORT_CHARSYNTHESISF( CHF_FRA(a_dWLeft),
                       CHF_FRA(a_dWRight),
                       CHF_CONST_FRA(a_W),
                       CHF_CONST_INT(a_dir),
                       CHF_CONST_INT(iRho),
                       CHF_BOX(a_box));

  int iBGN = WRHO;
  int iEND = iBGN + numSlopes() - 1;
  FORT_FROMSLOPESTOVALUES( CHF_CONST_FRA(a_W),
                           CHF_FRA(a_dWLeft),
                           CHF_FRA(a_dWRight),
                           CHF_CONST_INT(iBGN),
                           CHF_CONST_INT(iEND),
                           CHF_BOX(a_box));
}

void EqSysMHDMF::vectorVars(Vector<int> & a_vars) const
{
  a_vars.clear();
  a_vars.push_back(UMOMX);
  a_vars.push_back(UBX);

  int index  = UNUM + UMOMX;

  if( m_iCorrectionPotential > 0 )
  {
    index++;
  }

  if( m_isTrModelSet == 1 )
  {
    index += m_pTrModel->numConservative();
  }

  if( m_isPIonsSet == 1 )
  {
    index += m_pPIons->numConservative();
  }

  for( int i=1; i<m_nFluids; i++ )
  {
    int ind = index + UNUM_E*(i-1);
    a_vars.push_back(ind);
  }
}

void EqSysMHDMF::velocityVars(Vector<int> & a_vars) const
{
  a_vars.clear();
  a_vars.push_back(UMOMX);

  int index  = UNUM + UMOMX;

  if( m_iCorrectionPotential > 0 )
  {
    index++;
  }

  if( m_isTrModelSet == 1 )
  {
    index += m_pTrModel->numConservative();
  }

  if( m_isPIonsSet == 1 )
  {
    index += m_pPIons->numConservative();
  }

  for( int i=1; i<m_nFluids; i++ )
  {
    int ind = index + UNUM_E*(i-1);
    a_vars.push_back(ind);
  }
}

                              // Calculate equation system specific source terms
void EqSysMHDMF::explicitSource(       FArrayBox & a_U,
                                       FArrayBox & a_S,
                                 const FArrayBox & a_W,
                                 const Real      & a_dt,
                                 const int       & a_level,
                                 const Box       & a_box )
{
  if( m_isTrModelSet == 1 )
  {
    int iRhoPI = ( m_isPIonsSet == 0 ) ? -1 : m_pPIons->primInterval().begin();
    int iRhoN1 = densityIndexPrim( 1 );

    m_pTrModel->explicitSource( a_U, a_S, a_W, a_dt, iRhoN1, iRhoPI, a_level, a_box );
  }

  if( m_isPIonsSet == 1 )
  {
    m_pPIons->explicitSource( a_U, a_S, a_W, a_dt, a_level, a_box );
  }
    
  if( m_iCorrectionPotential > 0 )
  {
    int iRho = densityIndexCons(0);
    //FORT_DEDNERSOURCETERMS(
    //     CHF_FRA(a_S),
    //     CHF_CONST_FRA(a_W),
    //     CHF_CONST_INT(iRho),
    //     CHF_CONST_REAL(a_dt),     
    //     CHF_BOX(a_box));
    int iBgn = UENG;
    int iEnd =  iRho+8;
    //FORT_ADDSOURCES( CHF_FRA(a_U),
    //             CHF_CONST_FRA(a_S),
    //             CHF_CONST_INT(iBgn),
    //             CHF_CONST_INT(iEnd),
    //             CHF_BOX(a_box) );         
  }
}

void EqSysMHDMF::explicitSource(       FArrayBox & a_U,
                                       FArrayBox & a_S,
                                 const FArrayBox & a_W,
                              const BaseFab<int> & a_REG,
                                 const FArrayBox & a_divU,
                                 const Real      & a_dt,
                                 const int       & a_level,
                                 const Box       & a_box )
{
  if( m_isTrModelSet == 1 )
  {
    int iRhoPI = ( m_isPIonsSet == 0 ) ? -1 : m_pPIons->primInterval().begin();
    int iRhoN1 = densityIndexPrim( 1 );

    m_pTrModel->explicitSource( a_U, a_S, a_W, a_dt, iRhoN1, iRhoPI, a_level, a_box );
  }

  if( m_isPIonsSet == 1 )
  {
    m_pPIons->explicitSource( a_U, a_S, a_W, a_REG, a_dt, a_divU, a_level, a_box );
//    m_pPIons->explicitSource( a_U, a_S, a_W, a_dt, a_level, a_box );
 }
    
  if( m_iCorrectionPotential > 0 )
  {
    int iRho = densityIndexCons(0);
    //FORT_DEDNERSOURCETERMS(
    //     CHF_FRA(a_S),
    //     CHF_CONST_FRA(a_W),
    //     CHF_CONST_INT(iRho),
    //     CHF_CONST_REAL(a_dt),     
    //     CHF_BOX(a_box));
    int iBgn = UENG;
    int iEnd =  iRho+8;
    //FORT_ADDSOURCES( CHF_FRA(a_U),
    //             CHF_CONST_FRA(a_S),
    //             CHF_CONST_INT(iBgn),
    //             CHF_CONST_INT(iEnd),
    //             CHF_BOX(a_box) );         
  }
}
//    Compute equation system specific dt and returns a cell with the minimum dt
Real EqSysMHDMF::computeDt( const FArrayBox& a_U,
                                  FArrayBox& a_dt,
                            const Box      & a_box,
                            const int      & a_level,
                                  IntVect&   a_minDtCell )
{
  if( m_isPIonsSet == 1 )
  {
//take out
//try no separate PI time
//    m_pPIons->computeDt( a_U, a_dt, a_level, a_box );
  }
}

//                                               Set the turbulence model object
void EqSysMHDMF::setTurbulenceModel( TurbulenceModel* a_pTM )
{
  if( m_isTrModelSet == 1 )
  {
    int iTMCons  = m_pTrModel->numConservative();
    int iTMPrim  = m_pTrModel->numPrimitives();

    m_consStateInterval.define( 0, m_consStateInterval.size() - iTMCons - 1 );
    m_nPrim -= iTMPrim;

    if( m_isPIonsSet == 1 )
    {
      m_pPIons->adjustIntervals( m_pPIons->consInterval().begin() - iTMPrim,
                                 m_pPIons->primInterval().begin() - iTMCons );
    }

    if( m_nTrackingSurfaces > 0 )
    {
      m_lvlsStateInterval.define( m_lvlsStateInterval.begin() - iTMCons, m_lvlsStateInterval.end() - iTMCons );
      m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  - iTMPrim, m_lvlsPrimInterval.end()  - iTMPrim );
    }

    m_isTrModelSet = 0;
    m_pTrModel     = NULL;
  }

  if( a_pTM != NULL )
  {
    m_pTrModel     = a_pTM;
    m_isTrModelSet = 1;

    int iTMCons  = m_pTrModel->numConservative();
    int iTMPrim  = m_pTrModel->numPrimitives();

    m_consStateInterval.define( 0, m_consStateInterval.size() + iTMCons - 1 );
    m_nPrim += iTMPrim;

    if( m_isPIonsSet == 1 )
    {
      m_pPIons->adjustIntervals( m_pPIons->consInterval().begin() + iTMPrim,
                                 m_pPIons->primInterval().begin() + iTMCons );
    }

    if( m_nTrackingSurfaces > 0 )
    {
      m_lvlsStateInterval.define( m_lvlsStateInterval.begin() + iTMCons, m_lvlsStateInterval.end() + iTMCons );
      m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  + iTMPrim, m_lvlsPrimInterval.end()  + iTMPrim );
    }

    iTMCons  = UNUM;
    iTMPrim  = WNUM;

    if( m_iCorrectionPotential != 0 )
    {
      iTMCons++;
      iTMPrim++;
    }

    m_pTrModel->adjustIntervals( iTMCons, iTMPrim );
  }
}

//                                                    Set the pickup ions object
void EqSysMHDMF::setPickupIons( PickupIons* a_pPI )
{
  if( m_isPIonsSet == 1 )
  {
    int iPICons  = m_pPIons->numConservative();
    int iPIPrim  = m_pPIons->numPrimitives();

    m_consStateInterval.define( 0, m_consStateInterval.size() - iPICons - 1 );
    m_nPrim -= iPIPrim;

    if( m_nTrackingSurfaces > 0 )
    {
      m_lvlsStateInterval.define( m_lvlsStateInterval.begin() - iPICons, m_lvlsStateInterval.end() - iPICons );
      m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  - iPIPrim, m_lvlsPrimInterval.end()  - iPIPrim );
    }

    m_isPIonsSet   = 0;
    m_pPIons       = NULL;
  }

  if( a_pPI != NULL )
  {
    m_pPIons       = a_pPI;
    m_isPIonsSet   = 1;

    int iPICons  = m_pPIons->numConservative();
    int iPIPrim  = m_pPIons->numPrimitives();

    m_consStateInterval.define( 0, m_consStateInterval.size() + iPICons - 1 );
    m_nPrim += iPIPrim;

    if( m_nTrackingSurfaces > 0 )
    {
      m_lvlsStateInterval.define( m_lvlsStateInterval.begin() + iPICons, m_lvlsStateInterval.end() + iPICons );
      m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  + iPIPrim, m_lvlsPrimInterval.end()  + iPIPrim );
    }

    iPICons  = UNUM;
    iPIPrim  = WNUM;

    if( m_iCorrectionPotential != 0 )
    {
      iPICons++;
      iPIPrim++;
    }

    if( m_isTrModelSet == 1 )
    {
      iPICons += m_pTrModel->numConservative();
      iPIPrim += m_pTrModel->numPrimitives();
    }

    m_pPIons->adjustIntervals( iPICons, iPIPrim );
  }
}

//                Set usage of the correction potential approach by Dedner et al
void EqSysMHDMF::setCorrectionPotentialParams( int a_iUse, Real a_factorCh, Real a_factorCp )
{
  if( m_iCorrectionPotential > 0 )
  {
    if( a_iUse == 0 )
    {
//                                                               Simple cleaning
      m_consStateInterval.define( 0, m_consStateInterval.end() - 1 );
      m_nPrim--;

      if( m_isTrModelSet == 1 )
      {
        m_pTrModel->adjustIntervals( m_pTrModel->consInterval().begin() - 1,
                                     m_pTrModel->primInterval().begin() - 1 );
      }

      if( m_isPIonsSet == 1 )
      {
        m_pPIons->adjustIntervals( m_pPIons->consInterval().begin() - 1,
                                   m_pPIons->primInterval().begin() - 1 );
      }

      if( m_nTrackingSurfaces > 0 )
      {
        m_lvlsStateInterval.define( m_lvlsStateInterval.begin() - 1, m_lvlsStateInterval.end() - 1 );
        m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  - 1, m_lvlsPrimInterval.end()  - 1 );
      }

      m_iCorrectionPotential = 0;
    } else {
//                                                      Reuse with new constants
      m_iCorrectionPotential = a_iUse;
      m_dfactorCh = a_factorCh;
      m_dfactorCp = a_factorCp;      
            
    }
  } else {
    if( a_iUse != 0 )
    {
//                                                              First time usage
      m_consStateInterval.define( 0, m_consStateInterval.end() + 1 );
      m_nPrim++;

      if( m_isTrModelSet == 1 )
      {
        m_pTrModel->adjustIntervals( m_pTrModel->consInterval().begin() + 1,
                                     m_pTrModel->primInterval().begin() + 1 );
      }

      if( m_isPIonsSet == 1 )
      {
        m_pPIons->adjustIntervals( m_pPIons->consInterval().begin() + 1,
                                   m_pPIons->primInterval().begin() + 1 );
      }

      if( m_nTrackingSurfaces > 0 )
      {
        m_lvlsStateInterval.define( m_lvlsStateInterval.begin() + 1, m_lvlsStateInterval.end() + 1 );
        m_lvlsPrimInterval.define ( m_lvlsPrimInterval.begin()  + 1, m_lvlsPrimInterval.end()  + 1 );
      }

      m_iCorrectionPotential = a_iUse;
      m_dfactorCh = a_factorCh;
      m_dfactorCp = a_factorCp;      

    }
  }
}

// Set usage of the correction potential approach by Dedner et al.  
void EqSysMHDMF::getCorrectionPotentialParams( Real & a_factorCh, Real & a_factorCp )
{
  a_factorCh = m_dfactorCh;
  a_factorCp = m_dfactorCp;      
}

void EqSysMHDMF::setCorrectionPotential( double a_dCh, double a_dCp )
{
  m_dCh            = a_dCh;  
  m_dCp            = a_dCp;  
      
  FORT_SETDEDNERCONST( 
       CHF_CONST_REAL( m_dCh  ),
       CHF_CONST_REAL( m_dCp  ));
}


// Returns index for correction potential variable, or -1 if the correction potential is not used
int  EqSysMHDMF::correctionPotentialIndex (void)
{
  if (m_iCorrectionPotential > 0)  return m_UCP; else return -1;  
}

//                     Is the correction potential approach by Dedner et al used
int EqSysMHDMF::isUsedCorrectionPotential( void )
{
  return m_iCorrectionPotential;
}
