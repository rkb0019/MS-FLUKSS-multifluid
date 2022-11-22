#ifndef _EDGEBOX_H_
#define _EDGEBOX_H_

#include "Box.H"
#include "Vector.H"
#include "FArrayBox.H"
#include "NamespaceHeader.H"

/// A FArrayBox-like container for edge-centered data
/** This is a class to contain edge-centered data on a box.
    The class is a combination in some sence of two Chombo's classes: FluxBox and NodeFArrayBox.

    In the 3D case it is a simple modification of FluxBox class. But in the 2D case edge data
    ususally reduces to node data and can be implemented as a single FArrayBox. It means that
    in the 2D case EdgeBox class reduces to NodeFArrayBox.
 */
class EdgeBox
{

public:
  /// Default constructor
  EdgeBox();

  /// Constructs EdgeBox on cell-centered box with n components
  EdgeBox( const Box & bx,
                 int   n   = 1 );

  /// Destructor
  ~EdgeBox();

  /// Resize EdgeBox similar to BaseFab::resize()
  void resize( const Box& bx,
                     int  n   = 1 );

  /// Define function
  void define( const Box& bx,
                     int  n   = 1 );

  /// Returns the EdgeBox to the undefined state
  void clear();

  //  access functions

  /// Number of components
  int nComp() const;

  /// Returns cell-centered box which defines EdgeBox
  const Box& box() const;

  /// Returns edge-centered Data in direction \em dir
  FArrayBox& getData( const int dir );

  /// Returns const reference to edge-centered data in direction \em dir
  const FArrayBox& getData( const int dir ) const;

  /// Returns FArrayBox in direction dir
  FArrayBox& operator[] ( const int dir );

  /// Constant version
  const FArrayBox& operator[] ( const int dir ) const;

  //  data modification functions

  /// Set all edge data to val
  void setVal( const Real val );

  /// Set edge data in direction dir to val
  void setVal( const Real val,
               const int  dir );

  /// More specific setVal
  void setVal( const Real val,
               const int  dir,
               const int  startComp,
               const int  nComp     );

  /// Sets data on edges surrounding cell-centered box bx
  void setVal( const Real val,
               const Box& bx   );

  /// Most specific setVal
  /**
     Sets data on edges surrounding cell-centered box bx
  */
  void setVal( const Real val,
               const Box& bx,
               const int  dir,
               const int  startComp,
               const int  nComp       );

  /// Copy from src to this EdgeBox -- sizes must be identical
  void copy( const EdgeBox& src );

  /// Copy on overlap, for all directions
  void copy( const EdgeBox& src,
             const int      srcComp,
             const int      destComp,
             const int      numComp   );

  /// Copy on overlap of EdgeBoxes, in direction dir
  void copy( const EdgeBox& src,
             const int      dir,
             const int      srcComp,
             const int      destComp,
             const int      numComp   );

  /// Copies from a subsection of one box into another
  /**
     Assumes the boxes are both in the same index space, and that
     box R is completely contained in both the src and destination
     boxes.  This function required by BoxLayoutData
  */
  void copy( const Box&      R,
             const Interval& Cdest,
             const EdgeBox&  src,
             const Interval& Csrc  );

  /// Modifies this EdgeBox by copying the contents of \em src into it.
  /**
     This, the most general form
     of copy, specifies the contents of any sub-box srcbox
     in `EdgeBox' src may be copied into a (possibly
     different) destbox in the destination `EdgeBox'.  Note
     that although the srcbox and the destbox may be
     disjoint, they must be the same size and shape.  If the
     sizes differ, the copy is undefined and a runtime error
     results.  This copy function is the only one of the
     copy functions to allow a copy between differing
     boxes. The user also specifies how many components are
     copied, starting at component srccomp in src and stored
     starting at component destcomp.  The results are
     UNDEFINED if the src and dest EdgeBoxes are the same and
     the srcbox and destbox overlap.
  */
  void copy( const Box&        srcbox,
             const Interval&   destcomps,
             const Box&        destbox,
             const EdgeBox&    src,
             const Interval&   srccomps);

  /// Modifies this EdgeBox to its additive inverse.
  /**
     Modifies this EdgeBox by replacing each value with its additive
     inverse, for the given range of components and within the given subbox.
     Returns *this.
  */
  EdgeBox& negate( const Box& subbox,
                   int        comp    = 0,
                   int        numcomp = 1 );

  /// Modifies this EdgeBox to its additive inverse.
  /**
     Modifies this EdgeBox by replacing each value with its additive
     inverse, for the given range of components over the whole domain of
     the EdgeBox.  Returns *this.
  */
  EdgeBox& negate( int comp,
                   int numcomp = 1) ;

  /// Modifies this EdgeBox to its additive inverse.
  /**
     Modifies this EdgeBox by replacing each value with its additive
     inverse for all components over the whole domain of the EdgeBox.
     Returns *this.
  */
  EdgeBox& negate();

  /// Modifies this EdgeBox by adding the scalar Real r to all values.
  EdgeBox& operator+=( Real r );

  /// Modifies this EdgeBox by incrementing with the argument EdgeBox.
  /**
     Modifies this EdgeBox by pointwise addition of the values of the
     argument EdgeBox.  You might come to grief if the domains of the
     FArrayBoxes don't match, just as in FArrayBox::plus().
  */
  EdgeBox& operator+=( const EdgeBox& f );

  /// Modifies this EdgeBox by subtracting the scalar Real r to all values.
  EdgeBox& operator-=( Real r );

  /// Modifies this EdgeBox by decrementing with the argument EdgeBox.
  /**
     Modifies this EdgeBox by pointwise subtraction of the values of the
     argument EdgeBox.  You might come to grief if the domains of the
     EdgeBoxes don't match, just as in FArrayBox::minus().
  */
  EdgeBox& operator-=( const EdgeBox& f );

  /// Modifies this EdgeBox by multiplying all values by the scalar Real r
  EdgeBox& operator*=( Real r );

  /// Modifies this EdgeBox by multiplying by the argument EdgeBox.
  /**
     Modifies this EdgeBox by pointwise multiplication of the values by the
     argument EdgeBox.  You might come to grief if the domains of the
     EdgeBoxes don't match, just as in FArrayBox::mult().
  */
  EdgeBox& operator*=( const EdgeBox& f );

  /// Modifies this EdgeBox by shifting its domain box
  EdgeBox& shift( const IntVect& v );

  //  Linearization Functions

  /// Returns size of linearized data over \em bx.
  /**
     Returns size, in number of bytes, of a flat linear
     representation of data in components comps in faces around
     cell-centered box R
  */
  int size( const Box&      bx,
            const Interval& comps ) const;

  /// Writes a linear representation of this EdgeBox
  /**
     Write a linear representaion of the internal data for
     the faces surrounding cell-centered box R.  Assumes that sufficient
     memory for the buffer has already been allocated by the caller
  */
  void linearOut(       void*     buf,
                  const Box&      R,
                  const Interval& comps ) const;

  /// Read a linear representation of the data over the Box \rm R
  /**
      Reads in the output of linearOut
  */
  void linearIn(       void*     buf,
                 const Box&      R,
                 const Interval& comps );

  /// Helper function for linearization
  static int preAllocatable() {return 0;}


protected:
  /**
     Cell-centered Box over which this EdgeBox is defined
  */
  Box m_bx;

  /**
     Number of variables on each face
  */
  int m_nvar;

  /**
     Dimension of Vector<FArrayBox*> m_EdgeData
  */
  int m_iDim;

  /**
     CH_SPACEDIM FArrayBoxes which hold edge data
  */
  Vector<FArrayBox*> m_EdgeData;


private:
  // these are disallowed

  /// Disallowed
  EdgeBox( const EdgeBox& );
  /// Disallowed
  EdgeBox& operator =( const EdgeBox& );

};

#include "NamespaceFooter.H"
#endif
