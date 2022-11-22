#ifdef CH_LANG_CC
/*
*      _______              __
*     / ___/ /  ___  __ _  / /  ___
*    / /__/ _ \/ _ \/  V \/ _ \/ _ \
*    \___/_//_/\___/_/_/_/_.__/\___/
*    Please refer to Copyright.txt, in Chombo's root directory.
*
*    Extension of FineInterp class to support different coordinate systems
*/
#endif

#include "DisjointBoxLayout.H"
#include "LevelData.H"
#include "BaseFab.H"
#include "REAL.H"
#include "DataIterator.H"
#include "Tuple.H"
#include "InterpF_F.H"
#include "FineInterpMHDAMF_F.H"

#include "FineInterpMHDAM.H"
#include "CSHandler.H"

#include "NamespaceHeader.H"

FineInterpMHDAM::FineInterpMHDAM()
: m_iBX(-1)
{
  is_defined = false;
}

FineInterpMHDAM::~FineInterpMHDAM()
{         
}

FineInterpMHDAM::FineInterpMHDAM(const DisjointBoxLayout& a_fine_domain,
                       const int&  a_numcomps,                       
                       const int& a_ref_ratio,                       
                       const Box& a_fine_problem_domain,
                       const int& a_level,
                       CoordinateSystemHandler* a_csh)
: m_iBX(-1)
{
  is_defined = false;

  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, fineProbDomain, a_level, a_csh);
}

FineInterpMHDAM::FineInterpMHDAM(const DisjointBoxLayout& a_fine_domain,
                       const int&  a_numcomps,                       
                       const int& a_ref_ratio,                       
                       const ProblemDomain& a_fine_problem_domain,
                       const int& a_level,
                       CoordinateSystemHandler* a_csh)
: m_iBX(-1)
{
  is_defined = false;

  define(a_fine_domain, a_numcomps, a_ref_ratio, a_fine_problem_domain, a_level, a_csh);
}

void
FineInterpMHDAM::define(const DisjointBoxLayout& a_fine_domain,
                   const int& a_numcomps,                   
                   const int& a_ref_ratio,                   
                   const Box& a_fine_problem_domain,
                   const int& a_level,
                   CoordinateSystemHandler* a_csh)
{
  ProblemDomain fineProbDomain(a_fine_problem_domain);
  define(a_fine_domain, a_numcomps, a_ref_ratio, fineProbDomain, a_level, a_csh);
}

void
FineInterpMHDAM::define(const DisjointBoxLayout& a_fine_domain,
                   const int& a_numcomps,                   
                   const int& a_ref_ratio,                   
                   const ProblemDomain& a_fine_problem_domain,
                   const int& a_level,
                   CoordinateSystemHandler* a_csh)
{
  // check for consistency
  CH_assert (a_fine_domain.checkPeriodic(a_fine_problem_domain));
  m_ref_ratio = a_ref_ratio;
  m_crse_problem_domain = coarsen(a_fine_problem_domain, m_ref_ratio);
  
  m_level = a_level;
  m_csh   = a_csh;
  
  m_fine_problem_domain = a_fine_problem_domain;
  //
  // create the work array
  DisjointBoxLayout coarsened_fine_domain;
  coarsen ( coarsened_fine_domain,
            a_fine_domain,
            m_ref_ratio );
  m_coarsened_fine_data.define ( coarsened_fine_domain,
                                 a_numcomps,
                                 IntVect::Unit );
                                   
  is_defined = true;
  
  // Initialize aux data    
  for (int dir = SpaceDim; dir < 3; ++dir)
  {
    // Create fictions arrays
    m_dx[dir].define(Box(IntVect::Zero,IntVect::Zero), 1);
  }        
  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    IntVect iv_off(IntVect::Zero);
    iv_off[dir]=1;
    
    Box dxBox(m_crse_problem_domain.domainBox().smallEnd()*iv_off, 
              m_crse_problem_domain.domainBox().bigEnd()  *iv_off);
    dxBox.grow(dir,1);
    dxBox = dxBox & m_crse_problem_domain;
    
    m_dx[dir].define(dxBox,1);                     
    m_csh->dx(m_dx[dir],m_dx[dir].box(),dir,m_level);          
    
    Box crseCCbox(dxBox);    
    
    m_crseCC[dir].define(crseCCbox,1);    
    m_csh->getCellCenters(m_crseCC[dir],m_crseCC[dir].box(),dir,m_level);
      
    m_fineCC[dir].define(Box(
      m_fine_problem_domain.domainBox().smallEnd()*iv_off,
      m_fine_problem_domain.domainBox().bigEnd()  *iv_off),1);                
    m_csh->getCellCenters(m_fineCC[dir],m_fineCC[dir].box(),dir,m_level+1);
    
  }  
}

bool
FineInterpMHDAM::isDefined() const
{
  return ( is_defined );
}

// interpolate from coarse level to fine level
void
FineInterpMHDAM::interpToFine(LevelData<FArrayBox>& a_fine_data,
                         const LevelData<FArrayBox>& a_coarse_data)
{
 CH_assert(is_defined);
#ifndef NDEBUG
  // debugging check
  {
    DataIterator crseDit = m_coarsened_fine_data.dataIterator();
    for (crseDit.reset(); crseDit.ok(); ++crseDit)
      {
        m_coarsened_fine_data[crseDit()].setVal(1.0e9);
      }
  }
#endif

  // this should handle all the periodic BCs as well,
  // by filling in the ghost cells in an appropriate way
  a_coarse_data.copyTo(a_coarse_data.interval(),
                       m_coarsened_fine_data,
                       m_coarsened_fine_data.interval() );

  const BoxLayout fine_domain = a_fine_data.boxLayout();
  DataIterator dit = fine_domain.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
      const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
      BaseFab<Real>& fine = a_fine_data[dit()];

                       // interpGridData interpolates from an entire coarse grid
                       // onto an entire fine grid.
    interpGridData( fine, coarsened_fine, coarsened_fine_box, m_ref_ratio );
    }
}
void
FineInterpMHDAM::pwcinterpToFine(LevelData<FArrayBox>& a_fine_data,
                         const LevelData<FArrayBox>& a_coarse_data)
{
  CH_assert(is_defined);
  // this should handle all the periodic BCs as well,
  // by filling in the ghost cells in an appropriate way
  a_coarse_data.copyTo(a_coarse_data.interval(),
                       m_coarsened_fine_data,
                       m_coarsened_fine_data.interval() );

  const BoxLayout fine_domain = a_fine_data.boxLayout();
  DataIterator dit = fine_domain.dataIterator();

  for (dit.begin(); dit.ok(); ++dit)
    {
      const BaseFab<Real>& coarsened_fine = m_coarsened_fine_data[dit()];
      const Box& coarsened_fine_box = m_coarsened_fine_data.getBoxes()[dit()];
      BaseFab<Real>& fine = a_fine_data[dit()];
      // interpGridData interpolates from an entire coarse grid onto an
      // entire fine grid.
      pwcinterpGridData(fine,
                        coarsened_fine,
                        coarsened_fine_box,
                        m_ref_ratio);
    }
}

void
FineInterpMHDAM::pwcinterpGridData(BaseFab<Real>& a_fine,
                           const BaseFab<Real>& a_coarse,
                           const Box& a_coarsened_fine_box,
                           int a_ref_ratio) const
{
  // fill fine data with piecewise constant coarse data
  const Box& b = a_coarsened_fine_box;
  Box refbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);

  FORT_INTERPCONSTANT ( CHF_FRA(a_fine),
                        CHF_CONST_FRA(a_coarse),
                        CHF_BOX(b),
                        CHF_CONST_INT(a_ref_ratio),
                        CHF_BOX(refbox)
                        );
}
// interpolate from fine grid to coarse grid.  prerequisite:
// coarsened.box contains coarsen(fine.box).
//
// uses piecewise bilinear interpolation with multidimensional-limited
// slopes.  see design document for details.
void
FineInterpMHDAM::interpGridData(BaseFab<Real>& a_fine,
                           const BaseFab<Real>& a_coarse,
                           const Box& a_coarsened_fine_box,
                           int a_ref_ratio)
  const
{
  // fill fine data with piecewise constant coarse data
  const Box& b = a_coarsened_fine_box;
  const int num_comp = a_fine.nComp ();
  Box refbox(IntVect::Zero,
             (a_ref_ratio-1)*IntVect::Unit);

  FORT_INTERPCONSTANT ( CHF_FRA(a_fine),
                        CHF_CONST_FRA(a_coarse),
                        CHF_BOX(b),
                        CHF_CONST_INT(a_ref_ratio),
                        CHF_BOX(refbox)
                        );
  //  Tuple<BaseFab<Real>, SpaceDim> slopes;
  //  for (int dir = 0; dir < SpaceDim; ++dir)
  // hardwired to 3 due to lack of variable number of arguments in chfpp
  BaseFab<Real> slopes[3];
  for (int dir = 0; dir < 3; ++dir)
    {
      BaseFab<Real>& dir_slope = slopes[dir];
      dir_slope.resize(b, num_comp);
    }
  
  for (int dir = 0; dir < SpaceDim; ++dir)
  {    
    BaseFab<Real>& dir_slope = slopes[dir];        

    const Box bcenter = grow(m_crse_problem_domain,-BASISV(dir)) & b;
    if (!bcenter.isEmpty())
    {              
      if (m_csh->constStep(dir) == true)
      {
        Real dx = m_dx[dir].get(m_dx[dir].smallEnd(),0);
        FORT_INTERPCENTRALDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( bcenter ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( dx )
                        );
      } else
      {          
        FORT_INTERPCENTRALDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( bcenter ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( m_dx[dir], 0 )
                        );                   
      }        
    }
    const Box blo = b & adjCellLo(grow(m_crse_problem_domain,-BASISV(dir)),dir);
    if (!blo.isEmpty())
    {
      if (m_csh->constStep(dir) == true)
      {                    
        Real dx = m_dx[dir].get(m_dx[dir].smallEnd(),0);
        FORT_INTERPHISIDEDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( blo ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( dx )
                        );
      } else
      {          
        FORT_INTERPHISIDEDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( blo ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( m_dx[dir], 0 )
                        );                   
      }                
    }
    const Box bhi = b & adjCellHi(grow(m_crse_problem_domain,-BASISV(dir)),dir);
    if (!bhi.isEmpty())
    {
      if (m_csh->constStep(dir) == true)
      {     
        Real dx = m_dx[dir].get(m_dx[dir].smallEnd(),0);
        FORT_INTERPLOSIDEDERIV_C ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( bhi ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_REAL ( dx )
                        );
      } else
      {          
        FORT_INTERPLOSIDEDERIV_V ( CHF_FRA ( dir_slope ),
                        CHF_CONST_FRA ( a_coarse ),
                        CHF_BOX ( bhi ),
                        CHF_CONST_INT ( dir ),
                        CHF_CONST_FRA1 ( m_dx[dir], 0 )
                        );                   
      }                          
    }    
  }
    

  // to do limits, we need to have a box which includes
  // the neighbors of a given point (to check for the
  // local maximum...
  Box neighborBox(-1*IntVect::Unit,
                  IntVect::Unit);

  // GHM 7/12/01
  // interplimit iterates over box b_mod (was b), but cells within
  // 1 of the physical boundary never enter result (and this
  // wasted calculation may call upon uninitialized memory).
  // DFM 10/8/01
  // note that this turns off slope limiting for cells adjacent to the
  // boundary -- may want to revisit this in the future
  Box b_mod(b);
  b_mod.grow(1);
  b_mod = m_crse_problem_domain & b_mod;
  b_mod.grow(-1);

  // create a box grown big enough to remove periodic BCs from domain
  Box domBox = grow(b, 2);
  domBox = m_crse_problem_domain & domBox;
  
  
  FORT_INTERPLIMIT_V ( CHF_FRA ( slopes[0] ),
                   CHF_FRA ( slopes[1] ),
                   CHF_FRA ( slopes[2] ),
                   CHF_CONST_FRA ( a_coarse ),
                   CHF_CONST_FRA1 ( m_dx[0], 0 ),
                   CHF_CONST_FRA1 ( m_dx[1], 0 ),
                   CHF_CONST_FRA1 ( m_dx[2], 0 ),                   
                   CHF_BOX ( b_mod ),
                   CHF_BOX ( neighborBox ),
                   CHF_BOX (domBox)                   
                   );  

  if( m_iBX > -1 && m_iBX < num_comp )
  {
    FORT_INTERP_DIVB_0( CHF_FRA ( slopes[0] ),
                        CHF_FRA ( slopes[1] ),
                        CHF_FRA ( slopes[2] ),
                        CHF_BOX ( b_mod ),
                        CHF_CONST_INT ( m_iBX )
                      );
  }

  FArrayBox volc,volf;
  volc.define(m_crse_problem_domain & b,1);
  volf.define(m_fine_problem_domain & a_fine.box(),1);
    
  m_csh->getCellVolumes(volc,volc.box(),m_level);
  m_csh->getCellVolumes(volf,volf.box(),m_level+1);

  for (int dir = 0; dir < SpaceDim; ++dir)
  {      
    BaseFab<Real>& dir_slope = slopes[dir];
    
    FORT_INTERPLINEAR_V ( CHF_FRA ( a_fine ),
                          CHF_CONST_FRA ( dir_slope ),
                          CHF_CONST_FRA1 ( m_crseCC[dir], 0 ),
                          CHF_CONST_FRA1 ( m_fineCC[dir], 0 ),                          
                          CHF_CONST_FRA1 ( volc, 0 ),
                          CHF_CONST_FRA1 ( volf, 0 ),
                          CHF_BOX ( b ),
                          CHF_CONST_INT ( dir ),                       
                          CHF_CONST_INT ( a_ref_ratio ),                          
                          CHF_BOX ( refbox )
                          );    
  }

/*
  if( m_iBX > -1 && m_iBX < num_comp )
  {
    for( int dir = 0; dir < SpaceDim; ++dir )
    {
      BaseFab<Real>& dir_slope = slopes[dir];

      b_mod  = b;
      b_mod.growHi( dir, -1 );

      FORT_CORRECT_DIVB( CHF_FRA      ( a_fine ),
                         CHF_CONST_FRA( dir_slope ),
                         CHF_BOX      ( b_mod ),
                         CHF_CONST_INT( dir ),
                         CHF_CONST_INT( a_ref_ratio ),
                         CHF_CONST_INT( m_iBX )
                       );
    }

    int iBX = 5;
    b_mod   = a_fine.box();
    b_mod.grow( -1 );
    FORT_CHECKDIVB( CHF_CONST_FRA(a_fine), CHF_CONST_INT(iBX), CHF_BOX(b_mod) );

  }
*/
}
#include "NamespaceFooter.H"
