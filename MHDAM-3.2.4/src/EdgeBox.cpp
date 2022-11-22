#include "EdgeBox.H"
#include "NamespaceHeader.H"

                                                          // Default constructor
EdgeBox :: EdgeBox()
{
  m_iDim = (SpaceDim == 3) ? 3 : 1;

  m_EdgeData.resize( m_iDim, NULL );

  m_nvar = -1;
}
                    // Constructs EdgeBox on cell-centered box with n components
EdgeBox :: EdgeBox( const Box& a_bx,
                          int  a_nComp )
{
  m_iDim = (SpaceDim == 3) ? 3 : 1;

  m_EdgeData.resize( m_iDim, NULL );

  define( a_bx, a_nComp );
}
                                                                   // Destructor
EdgeBox :: ~EdgeBox()
{
  clear();
}
                                                         // Number of components
int EdgeBox :: nComp() const
{
  return m_nvar;
}
                              // Returns cell-centered box which defines EdgeBox
const Box& EdgeBox :: box() const
{
  return m_bx;
}
                              // Returns edge-centered Data in direction \em dir
FArrayBox & EdgeBox :: getData( const int dir )
{
  CH_assert(m_nvar > 0);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  return *m_EdgeData[iDir];
}
           // Returns const reference to edge-centered data in direction \em dir
const FArrayBox& EdgeBox :: getData( const int dir ) const
{
  CH_assert(m_nvar > 0);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  return *m_EdgeData[iDir];
}
                                           // Returns FArrayBox in direction dir
FArrayBox& EdgeBox :: operator[]( const int dir )
{
  CH_assert(m_nvar > 0);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  return *m_EdgeData[iDir];
}
                         // Returns FArrayBox in direction dir. Constant version
const FArrayBox& EdgeBox :: operator[]( const int dir ) const
{
  CH_assert(m_nvar > 0);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  return *m_EdgeData[iDir];
}
                                   // Returns the EdgeBox to the undefined state
void EdgeBox :: clear()
{
                                                         // first delete storage
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    if( m_EdgeData[dir] != NULL )
    {
      delete m_EdgeData[dir];
      m_EdgeData[dir] = NULL;
    }
  }
                                                // now reset all other variables
  m_nvar = -1;
  m_iDim = 0;
                                              // set the box to the empty box...
  m_bx = Box();
}
                                                              // Define function
void EdgeBox :: define( const Box& a_bx,
                              int  a_nComp )
{
  CH_assert(a_nComp > 0);

  m_bx = a_bx;
  m_nvar = a_nComp;

  if( m_iDim == 0 )
  {
    m_iDim = (SpaceDim == 3) ? 3 : 1;
    m_EdgeData.resize( m_iDim, NULL );
  }

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    if( m_EdgeData[dir] != NULL )
    {
      delete m_EdgeData[dir];
    }

    Box edgeBox( surroundingNodes( m_bx ) );
    if( SpaceDim == 3 )
      edgeBox.enclosedCells( dir );

    m_EdgeData[dir] = new FArrayBox( edgeBox, m_nvar );
  }
}
                                  // Resize EdgeBox similar to BaseFab::resize()
   // should resize fluxes in space (could be faster than re-allocating storage)
void EdgeBox :: resize( const Box& a_bx,
                              int  a_nComp )
{
                 // if this object has not already been defined, call define fn.
  if( m_nvar < 0 )
  {
    define(a_bx, a_nComp);
  }
  else
  {
    CH_assert(a_nComp > 0);

    m_bx   = a_bx;
    m_nvar = a_nComp;

    for( int dir = 0; dir < m_iDim; dir++ )
    {
      Box edgeBox( surroundingNodes( m_bx ) );
      if( SpaceDim == 3 )
        edgeBox.enclosedCells( dir );

      if( m_EdgeData[dir] != NULL )
      {
        m_EdgeData[dir]->resize( edgeBox, m_nvar );
      }
      else
      {
        FArrayBox* newFabPtr = new FArrayBox( edgeBox, m_nvar );
        m_EdgeData[dir] = newFabPtr;
      }
    }
  }
}
                                                     // Set all edge data to val
void EdgeBox :: setVal( const Real val )
{
  CH_assert(m_nvar > 0);

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->setVal( val );
  }
}
                                        // Set edge data in direction dir to val
void EdgeBox :: setVal( const Real val,
                        const int  dir )
{
  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);

  CH_assert(m_EdgeData[iDir] != NULL);
  m_EdgeData[iDir]->setVal( val );
}
                                                         // More specific setVal
void EdgeBox :: setVal( const Real val,
                        const int  dir,
                        const int  startComp,
                        const int  nComp      )
{
  CH_assert(startComp >-1);
  CH_assert(startComp + nComp <= m_nvar);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  for( int comp = startComp; comp < startComp + nComp; comp++ )
  {
    m_EdgeData[iDir]->setVal( val, comp );
  }
}
                          // Sets data on edges surrounding cell-centered box bx
void EdgeBox :: setVal( const Real val,
                        const Box& bx   )
{
  CH_assert(m_bx.contains(bx));

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);
                                   // move cell-centered box to appropriate edge
    Box edgeBox( surroundingNodes( bx ) );
    if( SpaceDim == 3 )
      edgeBox.enclosedCells( dir );

    m_EdgeData[dir]->setVal( val, edgeBox, 0, m_nvar );
  }

}
                                                         // Most specific setVal
void EdgeBox :: setVal( const Real val,
                        const Box& bx,
                        const int  dir,
                        const int  startComp,
                        const int  nComp      )
{
  CH_assert(m_bx.contains(bx));
  CH_assert(startComp > -1);
  CH_assert(startComp + nComp <= m_nvar);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  Box edgeBox( surroundingNodes( bx ) );
  if( SpaceDim == 3 )
    edgeBox.enclosedCells( iDir );

  m_EdgeData[iDir]->setVal( val, edgeBox, startComp, nComp );
}
                     // Copy from src to this EdgeBox -- sizes must be identical
void EdgeBox :: copy( const EdgeBox& src )
{
  CH_assert(src.box() == m_bx);
  CH_assert(src.nComp() == m_nvar);

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->copy( src[dir] );
  }
}
                                          // Copy on overlap, for all directions
void EdgeBox :: copy( const EdgeBox& src,
                      const int      srcComp,
                      const int      destComp,
                      const int      numComp   )
{
                                      // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    const FArrayBox& srcFab = src[dir];
    m_EdgeData[dir]->copy( srcFab, srcComp, destComp, numComp );
  }
}
                           // Copy on overlap of EdgeBoxes, in direction dirvoid
void EdgeBox :: copy( const EdgeBox& src,
                 const int      dir,
                 const int      srcComp,
                 const int      destComp,
                 const int      numComp   )
{
                                      // to ensure that neither comp is negative
  CH_assert(srcComp*destComp > -1);
  CH_assert(srcComp+numComp <= src.nComp());
  CH_assert(destComp+numComp <= m_nvar);

  int iDir = (SpaceDim == 3) ? dir : 0;

  CH_assert(iDir < m_iDim);
  CH_assert(m_EdgeData[iDir] != NULL);

  const FArrayBox& srcFab = src[iDir];
  m_EdgeData[iDir]->copy( srcFab, srcComp, destComp, numComp );
}
                             // Copies from a subsection of one box into another
void EdgeBox :: copy( const Box&      R,
                      const Interval& Cdest,
                      const EdgeBox&  src,
                      const Interval& Csrc  )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    Box Redge( surroundingNodes( R ) );
    if( SpaceDim == 3 )
      Redge.enclosedCells( dir );

    const FArrayBox& srcFab = src[dir];
 // all this intersecting is necessary due to the edge-centered nature of things
    Redge &= srcFab.box();
    Redge &= m_EdgeData[dir]->box();

    if( !Redge.isEmpty() )
    {
                               // this is probably wrong in periodic---dtg (???)
      m_EdgeData[dir]->copy( Redge, Cdest, Redge, srcFab, Csrc );
    }
  }
}
                 // Modifies this EdgeBox by copying the contents of src into it
void EdgeBox :: copy( const Box&      srcbox,
                      const Interval& destcomps,
                      const Box&      destbox,
                      const EdgeBox&  src,
                      const Interval& srccomps  )
{
  IntVect transVec( destbox.smallEnd() - srcbox.smallEnd() );

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    const FArrayBox& srcFab = src[dir];

    Box srcEdgeBox( surroundingNodes( srcbox ) );
    if( SpaceDim == 3 )
      srcEdgeBox.enclosedCells( dir );

    srcEdgeBox &= srcFab.box();

    Box destEdgeBox( surroundingNodes( destbox ) );
    if( SpaceDim == 3 )
      destEdgeBox.enclosedCells( dir );

    destEdgeBox &= (m_EdgeData[dir]->box());

    Box transBox( destEdgeBox );
    transBox.shift(-transVec );

    srcEdgeBox &= transBox;

    m_EdgeData[dir]->copy( srcEdgeBox, srccomps, destEdgeBox, srcFab, destcomps );
  }
}
                                // Modifies this EdgeBox to its additive inverse
EdgeBox& EdgeBox :: negate()
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->negate();
  }
  return *this;
}
                                // Modifies this EdgeBox to its additive inverse
EdgeBox& EdgeBox :: negate( int comp,
                            int numcomp )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->negate( comp, numcomp );
  }
  return *this;
}
                                // Modifies this EdgeBox to its additive inverse
EdgeBox& EdgeBox :: negate( const Box& subbox,
                                  int  comp,
                                  int  numcomp )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    Box edgeBox( surroundingNodes( subbox ) );
    if( SpaceDim == 3 )
      edgeBox.enclosedCells( dir );

    m_EdgeData[dir]->negate( edgeBox, comp, numcomp );
  }
  return *this;
}
              // Modifies this EdgeBox by adding the scalar Real r to all values
EdgeBox& EdgeBox :: operator+=( Real r )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->plus( r );
  }
  return *this;
}
              // Modifies this EdgeBox by incrementing with the argument EdgeBox
EdgeBox& EdgeBox :: operator+=( const EdgeBox& f )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->plus( f[dir] );
  }
  return *this;
}
         // Modifies this EdgeBox by subtracting the scalar Real r to all values
EdgeBox& EdgeBox :: operator-=( Real r )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->plus(-r );
  }
  return *this;
}
              // Modifies this EdgeBox by decrementing with the argument EdgeBox
EdgeBox& EdgeBox :: operator-=( const EdgeBox& f )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->minus( f[dir] );
  }
  return *this;
}
         // Modifies this EdgeBox by multiplying all values by the scalar Real r
EdgeBox& EdgeBox :: operator*=( Real r )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->mult( r );
  }
  return *this;
}
                 // Modifies this EdgeBox by multiplying by the argument EdgeBox
EdgeBox& EdgeBox :: operator*=( const EdgeBox& f )
{
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    m_EdgeData[dir]->mult( f[dir] );
  }
  return *this;
}
                             // Modifies this EdgeBox by shifting its domain box
EdgeBox& EdgeBox :: shift( const IntVect& iv )
{
  m_bx.shift( iv );

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    m_EdgeData[dir]->shift( iv );
  }

  return *this;
}
                                      // Returns size of linearized data over bx
int EdgeBox :: size( const Box&      bx,
                     const Interval& comps ) const
{
  int totalSize = 0;

  FArrayBox tempFab;
  for( int dir = 0; dir < m_iDim; dir++ )
  {
    Box edgeBox( surroundingNodes( bx ) );
    if( SpaceDim == 3 )
      edgeBox.enclosedCells( dir );

    const int dirSize = tempFab.size( edgeBox, comps );
    totalSize += dirSize;
  }

  return totalSize;
}
                               // Writes a linear representation of this EdgeBox
void EdgeBox :: linearOut(       void*     buf,
                           const Box&      R,
                           const Interval& comps ) const
{
  Real* buffer = (Real*) buf;

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    Box dirBox( surroundingNodes( R ) );
    if( SpaceDim == 3 )
      dirBox.enclosedCells( dir );

    int dirSize = m_EdgeData[dir]->size( dirBox, comps );
    m_EdgeData[dir]->linearOut( buffer, dirBox, comps );

    buffer += dirSize/sizeof(Real);
  }
}
                      // Read a linear representation of the data over the Box R
void EdgeBox :: linearIn(       void*     buf,
                          const Box&      R,
                          const Interval& comps )
{
  Real* buffer = (Real*) buf;

  for( int dir = 0; dir < m_iDim; dir++ )
  {
    CH_assert(m_EdgeData[dir] != NULL);

    Box dirBox( surroundingNodes( R ) );
    if( SpaceDim == 3 )
      dirBox.enclosedCells( dir );

    int dirSize = m_EdgeData[dir]->size( dirBox, comps) ;
    m_EdgeData[dir]->linearIn( buffer, dirBox, comps );

    buffer += dirSize/sizeof(Real);
  }
}
#include "NamespaceFooter.H"
