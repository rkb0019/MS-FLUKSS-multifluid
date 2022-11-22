#ifndef _FINE_INTERP_BCL_H_
#define _FINE_INTERP_BCL_H_

#include "FineInterp.H"

class DisjointBoxLayout;

/// replaces fine level data with interpolation of coarse level data.
/**
   This class replaces data at a fine level of refinement with data
   interpolated from a coarser level of refinement.  Interpolation is
   piecewise bi(tri)linear, with van Leer slopes if there is room for
   the stencil, with lower-order slopes if there isn't.  See the
   design document.

 */
class FineInterpBCL : public FineInterp
{
public:
  ///
  /**
     Default constructor.  User must subsequently call define().
  */
  FineInterpBCL();

  ///
  /**
     Destructor.
  */
  virtual ~FineInterpBCL();

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

  */
  FineInterpBCL(const DisjointBoxLayout& a_fine_domain,
                const int& a_numcomps,
                const int& a_ref_ratio,
                const Box& a_fine_problem_domain);

  ///
  /**
     Defining constructor.  Constructs a valid object.
     Equivalent to default construction followed by define().

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\
  */
  FineInterpBCL(const DisjointBoxLayout& a_fine_domain,
                const int& a_numcomps,
                const int& a_ref_ratio,
                const ProblemDomain& a_fine_problem_domain);

  /// 
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\ 
     ---This object is modified.---

  */
  virtual void 
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
         const Box& a_fine_problem_domain);      // problem domain


  /// 
  /**
     Defines this object.  Existing information is overriden.

     {\bf Arguments:}\\
     a_fine_domain (not modified): the fine level domain.\\
     a_numcomps (not modified): the number of components.\\
     a_ref_ratio (not modified): the refinement ratio.\\
     a_fine_problem_domain (not modified): problem domain at the fine level.\\

     {\bf This:}\\ 
     ---This object is modified.---

  */
  virtual void 
  define(const DisjointBoxLayout& a_fine_domain, // the fine level domain
         const int& a_numcomps,                   // the number of components
         const int& a_ref_ratio,                  // the refinement ratio
         const ProblemDomain& a_fine_problem_domain);


  ///
  /**
     Returns true if this object was created with the defining
     constructor or if define() has been called.

     {\bf This:}\\ 
     This object is not modified.
  */
  virtual bool 
  isDefined() const;

  /// 
  /**
     Replaces a_fine_data with data interpolated from a_coarse_data. It
     is an error to call if not this->isDefined().  The domain of
     a_fine_data should be the same as the fine domain specified in the
     most recent call to define().  It is expected that the coarse and
     fine level's domains are properly nested.  Both a_coarse_data and
     a_fine_data should have the same number of components specified in
     the most recent call to define().

     {\bf Arguments:}\\
     a_fine_data (modified): fine data. \\
     a_coarse_data (not modified): coarse data. \\

     {\bf This:}\\ 
     Well, it's complicated.  As far as the user is concerned, this object
     is not modified.  See the design document if you care for details.

  */
  virtual void 
  interpToFine(LevelData<FArrayBox>& a_fine_data,
               const LevelData<FArrayBox>& a_coarse_data);

protected:
  virtual void
  interpGridData(BaseFab<Real>& a_fine,
                 const BaseFab<Real>& a_coarse,
                 const Box& a_coarsened_fine_box,
                 int a_ref_ratio)
    const;
};

#endif
