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

#include "Box.H"
#include "FArrayBox.H"
#include "REAL.H"
#include "SPACE.H"
#include "Tuple.H"
#include "Vector.H"
#include "LoHiSide.H"

#include "PoissonBCMHD.H"
#include "MayDay.H"


// ================== FixedBC class ====================

FixedBC::FixedBC() : BoxGhostBC()
{
}

FixedBC::~FixedBC() 
{
}

void
FixedBC::fillBCValues(FArrayBox& a_neumfac,
                        FArrayBox& a_dircfac,
                        FArrayBox& a_inhmval,
                        Real a_dx,
                        const Box& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(1.0);
  a_dircfac.setVal(0.0);
  a_inhmval.setVal(0.0);
}


void
FixedBC::fillBCValues(FArrayBox& a_neumfac,
                        FArrayBox& a_dircfac,
                        FArrayBox& a_inhmval,
                        Real a_dx,
                        const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(1.0);
  a_dircfac.setVal(0.0);
  a_inhmval.setVal(0.0);
}

BoxGhostBC* 
FixedBC::new_boxghostbc() const
{
  FixedBC* newop = new FixedBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in FixedBC::new_boxghostbc");
    }
  return static_cast<BoxGhostBC*>(newop);
}

FixedBC::FixedBC(int dir, Side::LoHiSide side, const Interval& a_comps)
    : BoxGhostBC(dir,side,a_comps)
{
}


FixedBC::FixedBC(int dir, Side::LoHiSide side)
  : BoxGhostBC(dir,side)
{
}

// =================== AxisBC class =======================
void
AxisBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const Box& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(0.0);
}

void
AxisBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(0.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(0.0);
}


BoxGhostBC* 
AxisBC::new_boxghostbc() const
{
  AxisBC* newop = new AxisBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in AxisBC::new_boxghostbc");
    }
  return static_cast<BoxGhostBC*>(newop);
}



AxisBC::AxisBC() : BoxGhostBC()
{
}

AxisBC::~AxisBC()
{
}

AxisBC::AxisBC(int dir, Side::LoHiSide side)
    : BoxGhostBC(dir,side)
{
}

AxisBC::AxisBC(int dir, Side::LoHiSide side, const Interval& a_comps)
    : BoxGhostBC(dir,side,a_comps)
{
}


// =================== ContinuousBC class =======================
void
ContinuousBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const Box& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(1.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(0.0);
}

void
ContinuousBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const ProblemDomain& a_domain) const
{
  CH_assert(!a_neumfac.box().isEmpty());
  CH_assert(!a_dircfac.box().isEmpty());
  CH_assert(!a_inhmval.box().isEmpty());
  a_neumfac.setVal(1.0);
  a_dircfac.setVal(1.0);
  a_inhmval.setVal(0.0);
}


BoxGhostBC* 
ContinuousBC::new_boxghostbc() const
{
  ContinuousBC* newop = new ContinuousBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in ContinuousBC::new_boxghostbc");
    }
  return static_cast<BoxGhostBC*>(newop);
}



ContinuousBC::ContinuousBC() : BoxGhostBC()
{
}

ContinuousBC::~ContinuousBC()
{
}

ContinuousBC::ContinuousBC(int dir, Side::LoHiSide side)
    : BoxGhostBC(dir,side)
{
}

ContinuousBC::ContinuousBC(int dir, Side::LoHiSide side, const Interval& a_comps)
    : BoxGhostBC(dir,side,a_comps)
{
}


// =================== PeriodicBC class =======================
void
PeriodicBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const Box& a_domain) const
{
  MayDay::Error("PeriodicBC::fillBCValues must never be called");
}

void
PeriodicBC::fillBCValues(FArrayBox& a_neumfac,
                          FArrayBox& a_dircfac,
                          FArrayBox& a_inhmval,
                          Real a_dx,
                          const ProblemDomain& a_domain) const
{
  MayDay::Error("PeriodicBC::fillBCValues must never be called");
}


BoxGhostBC* 
PeriodicBC::new_boxghostbc() const
{
  PeriodicBC* newop = new PeriodicBC();
  if(newop == NULL)
    {
      MayDay::Error("Out of Memory in PeriodicBC::new_boxghostbc");
    }
  return static_cast<BoxGhostBC*>(newop);
}



PeriodicBC::PeriodicBC() : BoxGhostBC()
{
}

PeriodicBC::~PeriodicBC()
{
}

PeriodicBC::PeriodicBC(int dir, Side::LoHiSide side)
    : BoxGhostBC(dir,side)
{
}

PeriodicBC::PeriodicBC(int dir, Side::LoHiSide side, const Interval& a_comps)
    : BoxGhostBC(dir,side,a_comps)
{
}
