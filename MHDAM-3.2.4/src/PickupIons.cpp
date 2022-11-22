#include "PickupIons.H"

//                                                              Input parameters
void PickupIons::input( ParmParse & parser, int a_verbosity )
{
  m_verbosity  = a_verbosity;
}

//                                 Copy internal data to another PickupIons
void PickupIons::copyTo( PickupIons * pPIons ) const
{
  pPIons->m_verbosity         = m_verbosity;
  pPIons->m_consStateInterval = m_consStateInterval;
  pPIons->m_primStateInterval = m_primStateInterval;
}

//       Adjust intervals of conservative and primitive variables in state array
void PickupIons::adjustIntervals( int iConsFirst, int iPrimFirst )
{
  m_consStateInterval.define( iConsFirst, iConsFirst + m_consStateInterval.size() - 1 );
  m_primStateInterval.define( iPrimFirst, iPrimFirst + m_primStateInterval.size() - 1 );
}

//           Generate default names for the conservative variables, "tvariable#"
void PickupIons::AddConsNames( Vector<string> * pNames )
{
  int cnum = numConservative();

  for( int ivar = 0; ivar < cnum; ivar++ )
  {
    char varNameChar[80];
    sprintf( varNameChar, "tvariable%d", ivar );
    pNames->push_back( string(varNameChar) );
  }
}

//              Generate default names for the primitive variables, "tvariable#"
void PickupIons::AddPrimNames( Vector<string> * pNames )
{
  int cnum = numPrimitives();

  for( int ivar = 0; ivar < cnum; ivar++ )
  {
    char varNameChar[80];
    sprintf( varNameChar, "tvariable%d", ivar );
    pNames->push_back( string(varNameChar) );
  }
}

//                      Compute pickup ion model specific dt in all patch points
Real PickupIons::computeDt( const FArrayBox& a_U,
                                  FArrayBox& a_dt,
                            const int      & a_level,
                            const Box&       a_box )
{
}
