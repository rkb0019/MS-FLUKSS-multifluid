/* _______              __
  / ___/ /  ___  __ _  / /  ___
 / /__/ _ \/ _ \/  ' \/ _ \/ _ \
 \___/_//_/\___/_/_/_/_.__/\___/ 
*/
//
// This software is copyright (C) by the Lawrence Berkeley
// National Laboratory.  Permission is granted to reproduce
// this software for non-commercial purposes provided that
// this notice is left intact.
// 
// It is acknowledged that the U.S. Government has rights to
// this software under Contract DE-AC03-765F00098 between
// the U.S.  Department of Energy and the University of
// California.
//
// This software is provided as a professional and academic
// contribution for joint exchange. Thus it is experimental,
// is provided ``as is'', with no warranties of any kind
// whatsoever, no support, no promise of updates, or printed
// documentation. By using this software, you acknowledge
// that the Lawrence Berkeley National Laboratory and
// Regents of the University of California shall have no
// liability with respect to the infringement of other
// copyrights by any part of this software.
//
// DTGraves, Mon, July 19, 1999

#ifndef _POISSONBCMHD_H_
#define _POISSONBCMHD_H_

#include "Box.H"
#include "FArrayBox.H"
#include "Vector.H"
#include "LoHiSide.H"
#include "List.H"
//#include "GhostBC.H"

#include "ProblemDomain.H"

/// Class to enforce homogeneous Neumann boundary conditions 
/**
   BoxGhostBC class implentation to enforce homogeneous neumann 
   boundary conditions 
   Recall that  BoxGhostBC is a class to encapsulate the 
   operations of ghost-cell boundary conditions at a face.  
   If the solution is phi and the face normal direction is x,
   the boundary conditions usually used can be expressed as
   A*phi + B*dphi/dx = C.
   For homogeneous Neumann BCs, A=C=0, B= 1
   These functions are only meant to be called by DomainGhostBC            
 */

/// Class to encapsulate the operations  of ghost-cell boundary conditions.
/**
    BoxGhostBC is a class which encapsulates the
    operations  of ghost-cell boundary conditions at a single face.
    If the solution is phi and the face normal direction is x,
    the boundary conditions usually used can be expressed as
    A*phi + B*dphi/dx = C.
    The derived class provides a function to
    to fill A B and C.
    These functions are only meant to be called by DomainGhostBC.
*/
class BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Null constructor
  BoxGhostBC()
    :
    m_components(-1,-1)
  {
  }

  /// Destructor
  virtual ~BoxGhostBC()
  {
  }

  /// Sets interval to be (0,0)
  BoxGhostBC(int            a_dir,
             Side::LoHiSide a_sd) {};

  /// Full constructor
  BoxGhostBC(int             a_dir,
             Side::LoHiSide  a_sd,
             const Interval& a_comps) {};

protected:
  /**
     Virtual constructor workaround.
  */
  virtual BoxGhostBC* new_boxghostbc() const = 0;

  /**
     Sets interval to be (0,0)
  */
  void define(int            a_dir,
              Side::LoHiSide a_sd) {};

  ///
  /**
     Full define function
   */
  void define(int             a_dir,
              Side::LoHiSide  a_sd,
              const Interval& a_comps) {};

  ///
  /**
     Apply the inhomogeneous form of the boundary conditions to a_state
     along the boundary of a_domain. (Deprecated function -- use
     ProblemDomain instead of a Box to represent the domain)
   */
  virtual void applyInhomogeneousBCs(FArrayBox& a_state,
                                     const Box& a_domain,
                                     Real       a_dx) const {};

  ///
  /**
     Apply the homogeneous form of the boundary conditions to a_state
     along the boundary of a_domain.   (Deprecated function -- use
     ProblemDomain instead of a Box to represent the domain)
   */
  virtual void applyHomogeneousBCs(FArrayBox& a_state,
                                   const Box& a_domain,
                                   Real       a_dx) const {};

  ///
  /**
     Apply the inhomogeneous form of the boundary conditions to a_state
     along the boundary of a_domain.
   */
  virtual void applyInhomogeneousBCs(FArrayBox&           a_state,
                                     const ProblemDomain& a_domain,
                                     Real                 a_dx) const {};

  ///
  /**
     Apply the homogeneous form of the boundary conditions to a_state
     along the boundary of a_domain.
   */
  virtual void applyHomogeneousBCs(FArrayBox&           a_state,
                                   const ProblemDomain& a_domain,
                                   Real                 a_dx) const {};

  ///
  /**
      Deprecated function -- use ProplemDomain version instead
  */
  virtual void fillBCValues(FArrayBox& a_neumfac,
                            FArrayBox& a_dircfac,
                            FArrayBox& a_inhmval,
                            Real       a_dx,
                            const Box& a_domain) const = 0;

  ///
  /**
     Actually applies the boundary conditions defined by the values in
     a_neumfac, a_dircfac, and a_inhmval to a_state in the region
     defined by a_bcBox (which should be the appropriate ghost cells)
  */
  virtual void applyBCs(const Box&       a_bcbox,
                        FArrayBox&       a_state,
                        const FArrayBox& a_neumfac,
                        const FArrayBox& a_dircfac,
                        const FArrayBox& a_inhmval,
                        Real             a_dx) const {};

  ///
  /**
      Pure-virtual Worker function to fill arrays which define the
      boundary condition.
      Essentially, if the boundary condition is represented as
      A*phi + B*dphi/dn = C, then a_dircfac is A, a_neumfac is B, and
      a_inhmval is C.
  */
  virtual void fillBCValues(FArrayBox&           a_neumfac,
                            FArrayBox&           a_dircfac,
                            FArrayBox&           a_inhmval,
                            Real                 a_dx,
                            const ProblemDomain& a_domain) const = 0;

  /**
     Which side this boundary condition is on
  */
  Side::LoHiSide m_side;

  /**
     Which direction this boundary condition is for
  */
  int m_direction;

  /**
     Which components this boundary condition should be applied to
  */
  Interval m_components;

private:
  BoxGhostBC(const BoxGhostBC&)
    :
    m_components(-1,-1)
  {
  }

  virtual void operator=(const BoxGhostBC&)
  {
  }
};

/// Class to enforce ghost-cell based boundary conditions
/**
    This class holds a domain's worth of boxghostbcs,
    at least one for each face.  The user adds BoxGhostBC-derived
    objects for each face to enforce boundary conditions.
    If the solution is phi and the face normal direction is x,
    the boundary conditions at a face may often be expressed as
    A*phi + B*dphi/dx = C.
    The class derived from BoxGhostBC provides a function to
    to fill A B and C.  DomainGhostBC calls this function
    and fills the ghost cells outside the domain appropriately.
    Because BoxGhostBC operates on a component Interval, more than
    one BoxGhostBC per side is supported, to support the case where
    different components have different boundary conditions.
 */
class DomainGhostBC
{
public:
  /// Null constructor
  DomainGhostBC();

  /// Destructor
  ~DomainGhostBC();

  /// Set boundary conditions at a face to those defined in a_ghostBC
  /**
     Adds a_ghostBC to whatever BCs are already defined for the gviven
     side and face.
   */
  void setBoxGhostBC(const BoxGhostBC& a_ghostBC);

  /// Access function -- returns a reference to the BC on the given face.
  const List<BoxGhostBC*>& operator()(int            a_direction,
                                      Side::LoHiSide a_side) const;

  /// Deprecated function -- use ProblemDomain instead of a Box for the domain
  /**
  */
  void applyHomogeneousBCs(FArrayBox& a_state,
                           const Box& a_domain,
                           Real       a_dx) const;

  /// apply boundary conditions on all sides of state which abut domain edge
  /**
  */
  void applyHomogeneousBCs(FArrayBox&           a_state,
                           const ProblemDomain& a_domain,
                           Real                 a_dx) const;

  /// Deprecated function -- use ProblemDomain instead of a Box for the domain
  /**
  */
  void applyInhomogeneousBCs(FArrayBox& a_state,
                             const Box& a_domain,
                             Real       a_dx) const;

  /// Apply boundary conditions on all sides of state which abut the domain edge  /**
  /**
  */
  void applyInhomogeneousBCs(FArrayBox&           a_state,
                             const ProblemDomain& a_domain,
                             Real                 a_dx) const;

  /// Assignment operator
  DomainGhostBC& operator=(const DomainGhostBC& a_dgbcin);

  /// Copy constructor
  DomainGhostBC(const DomainGhostBC& a_dgbcin);

protected:
  /**
     Returns true if this GhostBC is defined
  */
  bool isBCDefined(const int            a_dir,
                   const Side::LoHiSide a_side) const;

  /**
     Resets BC in the given direction and side to an unset state
     (deletes any existing BC's and returns to an empty list)
  */
  void resetBoxGhostBC(const int            a_dir,
                       const Side::LoHiSide a_side);

  /// clear all BoxGhostBC's
  void clear();

  /**
     Contains the lo-side boundary conditions
  */
  Tuple<List<BoxGhostBC*>, SpaceDim> m_loGhostBC;
  /**
     Contains the hi-side boundary conditions
  */
  Tuple<List<BoxGhostBC*>, SpaceDim> m_hiGhostBC;

private:
};

 

class FixedBC : public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Default constructor
  FixedBC();

  /// Basic constructor -- sets Interval to (0,0)
  FixedBC(int dir, Side::LoHiSide sd);

  /// Full constructor
  FixedBC(int dir, Side::LoHiSide sd, 
            const Interval& a_comps);

  /// Destructor
  ~FixedBC();


protected:

  /// 
  /**
     Fills coefficient arrays -- 
     fills neumann factor with 1.0
     and other two fields with zero.
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const Box& a_domain) const;

  ///
  /**
     Fills coefficient arrays --
     fills neumann factor with 1.0
     and other two fields with zero.
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real a_dx,
               const ProblemDomain& a_domain) const;

  ///
  /**
     Virtual constructor workaround.
  */
  BoxGhostBC* new_boxghostbc() const;

private:
  FixedBC(const FixedBC& neuin){;}
  virtual void operator=(const FixedBC& neuin){;}

};

/// Class to enforce homogeneous Dirichlet boundary conditions 
/**
   BoxGhostBC class implentation to enforce homogeneous dirchelet 
   boundary conditions 
   Recall that
   BoxGhostBC is a class to encapsulate the 
   operations  of ghost-cell boundary conditions at a face.  
   If the solution is phi and the face normal direction is x,
   the boundary conditions usually used can be expressed as
   A*phi + B*dphi/dx = C.
   For homogeneous Dirichlet BCs, B=C=0, A= 1
   These functions are only meant to be called by DomainGhostBC
 */
class AxisBC : public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Default constructor
  AxisBC();

  /// Simple constructor (sets Interval to (0,0)
  AxisBC(int dir, Side::LoHiSide sd);

  /// Full constructor
  AxisBC(int dir, Side::LoHiSide sd, 
              const Interval& a_comps);

  /// Destructor
  ~AxisBC();

private:

  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const Box& a_domain) const;


  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const ProblemDomain& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;
private:
  AxisBC(const AxisBC& dircin){;}
  virtual void operator=(const AxisBC& dircin){;}

};

/// Class implements continuous boundary condition 
/**
   BoxGhostBC class implentation to enforce homogeneous dirchelet 
   boundary conditions 
   Recall that
   BoxGhostBC is a class to encapsulate the 
   operations  of ghost-cell boundary conditions at a face.  
   If the solution is phi and the face normal direction is x,
   the boundary conditions usually used can be expressed as
   A*phi + B*dphi/dx = C.
   For continuous BCs,A=B=1 C=0
   These functions are only meant to be called by DomainGhostBC
 */

class ContinuousBC : public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Default constructor
  ContinuousBC();

  /// Simple constructor (sets Interval to (0,0)
  ContinuousBC(int dir, Side::LoHiSide sd);

  /// Full constructor
  ContinuousBC(int dir, Side::LoHiSide sd, 
              const Interval& a_comps);

  /// Destructor
  ~ContinuousBC();

private:

  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const Box& a_domain) const;


  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const ProblemDomain& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;
private:
  ContinuousBC(const ContinuousBC& dircin){;}
  virtual void operator=(const ContinuousBC& dircin){;}

};

/// Class to enforce periodic boundary conditions 
/**
   BoxGhostBC has implementation of periodic BC
   Here we overide all abstract functions with empty functions
   to be able to use BoxGhostBC explicitly 
 */
class PeriodicBC : public BoxGhostBC
{
public:
  friend class DomainGhostBC;

  /// Default constructor
  PeriodicBC();

  /// Simple constructor (sets Interval to (0,0)
  PeriodicBC(int dir, Side::LoHiSide sd);

  /// Full constructor
  PeriodicBC(int dir, Side::LoHiSide sd, 
              const Interval& a_comps);

  /// Destructor
  ~PeriodicBC();

private:

  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const Box& a_domain) const;


  ///
  /**
     fills Dirichlet factor with 1.0
     and other two fields with zero
   */
  virtual void
  fillBCValues(FArrayBox& a_neumfac,
               FArrayBox& a_dircfac,
               FArrayBox& a_inhmval,
               Real dx,
               const ProblemDomain& a_domain) const;

  ///
  BoxGhostBC* new_boxghostbc() const;
private:
  PeriodicBC(const PeriodicBC& dircin){;}
  virtual void operator=(const PeriodicBC& dircin){;}

};


#endif
