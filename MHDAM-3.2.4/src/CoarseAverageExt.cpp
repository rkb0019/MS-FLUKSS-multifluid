#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*/
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "FArrayBox.H"
#include "DataIterator.H"
#include "LayoutIterator.H"
#include "parstream.H"
#include "CoarseAverageExtF_F.H"

#include "CoarseAverageExt.H"
#include "NamespaceHeader.H"

CoarseAverageExt::CoarseAverageExt()
  :
  is_defined(false)
{
}

CoarseAverageExt::~CoarseAverageExt()
{
}

CoarseAverageExt::CoarseAverageExt(const DisjointBoxLayout& a_fine_domain,
                             int a_numcomps,                             
                             int a_ref_ratio,
                             int a_level,
                             CoordinateSystemHandler* a_csh)
  :
  is_defined(false), m_is_copier_defined(false)
{
  define(a_fine_domain, a_numcomps, a_ref_ratio, a_level, a_csh);
}

CoarseAverageExt::CoarseAverageExt(const DisjointBoxLayout& a_fine_domain,
                             const DisjointBoxLayout& a_crse_domain,
                             int a_numcomps,
                             int a_ref_ratio,
                             int a_level,
                             CoordinateSystemHandler* a_csh)
  :
  is_defined(false), m_is_copier_defined(false)
{
  define(a_fine_domain, a_crse_domain, a_numcomps, a_ref_ratio, a_level, a_csh);
}

void
CoarseAverageExt::define(const DisjointBoxLayout& a_fine_domain,
                      int a_numcomps,                      
                      int a_ref_ratio,
                      int a_level,
                      CoordinateSystemHandler* a_csh)
{
  m_ref_ratio = a_ref_ratio;
  DisjointBoxLayout coarsened_fine_domain;
  coarsen (coarsened_fine_domain, a_fine_domain, m_ref_ratio);
  m_coarsened_fine_data.define (coarsened_fine_domain, a_numcomps);
  
  m_crseDomain = coarsened_fine_domain.physDomain();
  
  m_fineDomain = a_fine_domain.physDomain();
  
  m_csh    = a_csh;
  m_level  = a_level;
  
  m_is_copier_defined = false;
  is_defined = true;
}

void
CoarseAverageExt::define(const DisjointBoxLayout& a_fine_domain,
                      const DisjointBoxLayout& a_crse_domain,
                      int a_numcomps,
                      int a_ref_ratio,
                      int a_level,
                      CoordinateSystemHandler* a_csh)
{
  m_ref_ratio = a_ref_ratio;
  DisjointBoxLayout coarsened_fine_domain;
  coarsen (coarsened_fine_domain, a_fine_domain, m_ref_ratio);
  m_coarsened_fine_data.define (coarsened_fine_domain, a_numcomps);
  
  m_crseDomain = a_crse_domain.physDomain();
  
  m_fineDomain = a_fine_domain.physDomain();

  // also can pre-define copier here
  // note that since CoarseAverageExt only operates on interior
  // data, there is no need to worry about whether the domain
  // is periodic.
  m_copier.define(coarsened_fine_domain, a_crse_domain);
  m_is_copier_defined = true;
  
  m_csh    = a_csh;
  m_level  = a_level;

  is_defined = true;
}

bool
CoarseAverageExt::isDefined() const
{
  return ( is_defined );
}

void
CoarseAverageExt::averageToCoarse(LevelData<FArrayBox>& a_coarse_data,
                               const LevelData<FArrayBox>& a_fine_data)
{
  computeAverages(a_coarse_data, a_fine_data, arithmetic);
}

void 
CoarseAverageExt::averageToCoarseHarmonic(LevelData<FArrayBox>& a_coarse_data,
                                       const LevelData<FArrayBox>& a_fine_data)
{
  computeAverages(a_coarse_data, a_fine_data, harmonic);
}

void
CoarseAverageExt::computeAverages(LevelData<FArrayBox>& a_coarse_data,
                               const LevelData<FArrayBox>& a_fine_data,
                               int a_averageType)
{
  CH_assert(is_defined);
  // it would be nice if this could check for validity of a_fine_data.
  // this could be done with a redundant DisjointBoxLayout
  DataIterator dit = a_fine_data.boxLayout().dataIterator();
  for (dit.begin(); dit.ok(); ++dit)

    {
      BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
      const BaseFab<Real>& fine = a_fine_data[dit()];
      // coarsenGridData coarsens from the entire fine grid onto the entire
      // coarse grid.
      averageGridData(coarsened_fine,
                      fine,
                      m_ref_ratio,
                      a_averageType);
    }

  if (m_is_copier_defined)
    {
      // we can use the pre-defined copier to make things faster
      m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                   a_coarse_data,
                                   a_coarse_data.interval(),
                                   m_copier);
    }
  else
    {
      m_coarsened_fine_data.copyTo(m_coarsened_fine_data.interval(),
                                   a_coarse_data,
                                   a_coarse_data.interval() );
    }
}

void
CoarseAverageExt::averageGridData(BaseFab<Real>& a_coarse,
                               const BaseFab<Real>& a_fine,
                               int a_ref_ratio,
                               int a_averageType)
  const
{
  const Box& b = a_coarse.box();
  Box refbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);
             
  Box boxf(b);
  boxf.refine(a_ref_ratio);
  FArrayBox volf(boxf,1);
    
  m_csh->getCellVolumes(volf,boxf,m_level+1);
  
  FORT_AVERAGEWITHVOLUMES( CHF_FRA(a_coarse),
                    CHF_CONST_FRA(a_fine),                    
                    CHF_CONST_FRA1(volf,0),
                    CHF_BOX(b),
                    CHF_CONST_INT(a_ref_ratio),
                    CHF_BOX(refbox)
                    );               
  
}
#include "NamespaceFooter.H"
