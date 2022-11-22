#include "TurbulenceModel.H"

//                                                              Input parameters
void TurbulenceModel::input( ParmParse & parser, int a_verbosity )
{
  m_verbosity  = a_verbosity;
}

//                                 Copy internal data to another TurbulenceModel
void TurbulenceModel::copyTo( TurbulenceModel * pTModel ) const
{
  pTModel->m_verbosity         = m_verbosity;
  pTModel->m_consStateInterval = m_consStateInterval;
  pTModel->m_primStateInterval = m_primStateInterval;
}

//       Adjust intervals of conservative and primitive variables in state array
void TurbulenceModel::adjustIntervals( int iConsFirst, int iPrimFirst )
{
  m_consStateInterval.define( iConsFirst, iConsFirst + m_consStateInterval.size() - 1 );
  m_primStateInterval.define( iPrimFirst, iPrimFirst + m_primStateInterval.size() - 1 );
}

//           Generate default names for the conservative variables, "tvariable#"
void TurbulenceModel::AddConsNames( Vector<string> * pNames )
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
void TurbulenceModel::AddPrimNames( Vector<string> * pNames )
{
  int cnum = numPrimitives();

  for( int ivar = 0; ivar < cnum; ivar++ )
  {
    char varNameChar[80];
    sprintf( varNameChar, "tvariable%d", ivar );
    pNames->push_back( string(varNameChar) );
  }
}
