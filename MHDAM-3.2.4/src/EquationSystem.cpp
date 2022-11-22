#include "EquationSystem.H"


//                                         Flag everything as not defined or set
EquationSystem::EquationSystem()
{
  m_verbosity = 0;
  m_nTrackingSurfaces = 0;

  m_pTrModel     = NULL;
  m_isTrModelSet = 0;
  m_pPIons       = NULL;
  m_isPIonsSet   = 0;
}

EquationSystem::~EquationSystem()
{
}

//                          Define this object and the boundary condition object
void EquationSystem::define(const Interval & a_consStateInterval,              
              int              a_nPrim,
              int              a_nTrackingSurfaces,
              const Interval & a_lvlsStateInterval,      
              Interval       & a_lvlsPrimInterval)
{      
  CH_assert(a_lvlsStateInterval.size() == a_lvlsPrimInterval.size());
  
  m_nTrackingSurfaces = (a_lvlsStateInterval.size()>0 ? a_lvlsStateInterval.size() : 0);    
  
  m_consStateInterval = a_consStateInterval;
  m_nPrim             = a_nPrim;    
  m_lvlsStateInterval = a_lvlsStateInterval;
  m_lvlsPrimInterval  = a_lvlsPrimInterval;  
}

int EquationSystem::numStates()
{
  return numConserved()+numTrackingSurfaces();
}

//               Generate default names for the conserved variables, "variable#"
Vector<string> EquationSystem::stateNames()
{
  Vector<string> retval;

  int cnum = numStates();

  for (int ivar = 0; ivar < cnum; ivar++)
  {
    char varNameChar[80];
    sprintf(varNameChar,"variable%d",ivar);
    retval.push_back(string(varNameChar));
  }

  return retval;
}

//               Generate default names for the primitive variables, "variable#"
Vector<string> EquationSystem::primitiveNames()
{
  Vector<string> retval;

  int cnum = numPrimitives();

  for (int ivar = 0; ivar < cnum; ivar++)
  {
    char varNameChar[80];
    sprintf(varNameChar,"variable%d",ivar);
    retval.push_back(string(varNameChar));
  }

  return retval; 
}

              // Default implementation of equation system specific source terms
void EquationSystem::explicitSource(       FArrayBox & a_U,
                                           FArrayBox & a_S,
                                     const FArrayBox & a_W,
                                     const Real      & a_dt,
                                     const int       & a_level,
                                     const Box       & a_box )
{
}

//                                               Set the turbulence model object
void EquationSystem::setTurbulenceModel( TurbulenceModel* a_pTM )
{
  if( m_isTrModelSet == 1 )
  {
    m_consStateInterval.define( 0, m_consStateInterval.size() - m_pTrModel->numConservative() - 1 );
    m_nPrim -= m_pTrModel->numPrimitives();

    m_isTrModelSet = 0;
    m_pTrModel     = NULL;
  }

  if( a_pTM != NULL )
  {
    m_pTrModel     = a_pTM;
    m_isTrModelSet = 1;

    m_consStateInterval.define( 0, m_consStateInterval.size() + m_pTrModel->numConservative() - 1 );
    m_nPrim += m_pTrModel->numPrimitives();
  }
}

//    Compute equation system specific dt and returns a cell with the minimum dt
Real EquationSystem::computeDt( const FArrayBox& a_U,
                                      FArrayBox& a_dt,
                                const Box&       a_box,
                                const int      & a_level,
                                      IntVect&   a_minDtCell )
{
}

//                                               Get the turbulence model object
TurbulenceModel* EquationSystem::getTurbulenceModel( void ) const
{
  return m_pTrModel;
}

//                                                    Set the pickup ions object
void EquationSystem::setPickupIons( PickupIons* a_pPI )
{
  if( m_isTrModelSet == 1 )
  {
    m_consStateInterval.define( 0, m_consStateInterval.size() - m_pPIons->numConservative() - 1 );
    m_nPrim -= m_pPIons->numPrimitives();

    m_isPIonsSet   = 0;
    m_pPIons       = NULL;
  }

  if( a_pPI != NULL )
  {
    m_pPIons     = a_pPI;
    m_isPIonsSet = 1;

    m_consStateInterval.define( 0, m_consStateInterval.size() + m_pPIons->numConservative() - 1 );
    m_nPrim += m_pPIons->numPrimitives();
  }
}

//                                                    Get the pickup ions object
PickupIons* EquationSystem::getPickupIons( void ) const
{
  return m_pPIons;
}
