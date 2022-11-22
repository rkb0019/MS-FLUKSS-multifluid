#include <cmath>

#include "CHOMBO_VERSION.H"

#include "REAL.H"
#include "IntVect.H"
#include "Box.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "IntVectSet.H"
#include "DisjointBoxLayout.H"
#include "LayoutIterator.H"

#include "InterpF_F.H"
#include "FineInterpMHDAMF_F.H"

#include "CH_Timer.H"
#include "MayDay.H"
using std::cout;
using std::endl;

#include <execinfo.h> 



#include "DebugF_F.H"
#include "CSHandler.H"
#include "PiecewiseLinearFillPatchMHDAM.H"

#include "NamespaceHeader.H"

#ifndef copysign
template <class T>
inline
T
copysign (const T& a,
          const T& b)
{
    return ( b >= 0 ) ? ( ( a >= 0 ) ? a : -a) : ( (a >= 0 ) ? -a : a);
}

#endif

const int PiecewiseLinearFillPatchMHDAM::s_stencil_radius = 1;

PiecewiseLinearFillPatchMHDAM::PiecewiseLinearFillPatchMHDAM()
: m_iBX(-1)
{
  m_is_defined = false;
}

PiecewiseLinearFillPatchMHDAM::~PiecewiseLinearFillPatchMHDAM()
{
}

PiecewiseLinearFillPatchMHDAM::PiecewiseLinearFillPatchMHDAM(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const Box& a_crse_problem_domain,  
  int a_ref_ratio,
  int a_level,
  int a_interp_radius,
  CoordinateSystemHandler* a_csh)
    : m_iBX(-1)
{
  m_is_defined = false;

  ProblemDomain crseProbDomain(a_crse_problem_domain);
  define(a_fine_domain,
         a_coarse_domain,
         a_num_comps,
         crseProbDomain,         
         a_ref_ratio,
         a_level,
         a_interp_radius,
         a_csh);
}

PiecewiseLinearFillPatchMHDAM::PiecewiseLinearFillPatchMHDAM(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const ProblemDomain& a_crse_problem_domain,  
  int a_ref_ratio,
  int a_level,
  int a_interp_radius,
  CoordinateSystemHandler* a_csh)
    : m_iBX(-1)
{
  m_is_defined = false;

  define(a_fine_domain,
         a_coarse_domain,
         a_num_comps,
         a_crse_problem_domain,         
         a_ref_ratio,
         a_level,
         a_interp_radius,
         a_csh);
}

void
PiecewiseLinearFillPatchMHDAM::define(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const Box& a_crse_problem_domain,  
  int a_ref_ratio,
  int a_level,
  int a_interp_radius,
  CoordinateSystemHandler* a_csh  
  )
{
  ProblemDomain crseProbDomain(a_crse_problem_domain);
  define(a_fine_domain, a_coarse_domain, a_num_comps,
         crseProbDomain, a_ref_ratio, a_level, a_interp_radius, a_csh);
}

void
PiecewiseLinearFillPatchMHDAM::define(
  const DisjointBoxLayout& a_fine_domain,
  const DisjointBoxLayout& a_coarse_domain,
  int a_num_comps,
  const ProblemDomain& a_crse_problem_domain,  
  int a_ref_ratio,
  int a_level,
  int a_interp_radius,
  CoordinateSystemHandler* a_csh
  )
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::define");
  m_ref_ratio = a_ref_ratio;
  m_interp_radius = a_interp_radius;
  m_crse_problem_domain = a_crse_problem_domain;  
  m_level  = a_level;
  m_csh    = a_csh;
  

  const ProblemDomain  fine_problem_domain = refine(m_crse_problem_domain,
                                                    m_ref_ratio);

  ///////////////////////////
  /*pout() << "PiecewiseLinearFillPatchMHDAM::define(const DisjointBoxLayout& a_fine_domain,"  << endl;
  a_fine_domain.physDomain().domainBox().p();
  pout() << "size " << a_fine_domain.size() << endl;
  
  
  int j, nptrs;
  const int SIZE = 100;
  void *buffer[SIZE];
  char **strings;

  nptrs = backtrace(buffer, SIZE);

   // The call backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO)
   //    would produce similar output to the following: 

  strings = backtrace_symbols(buffer, nptrs);
  if (strings == NULL) {
      pout() << "backtrace_symbols error" << endl;        
  } else
  {
   for (j = 0; j < nptrs; j++)
        pout() << std::string(strings[j]) << endl;        

   free(strings);
  }*/
   
  ///////////////////////////

  
  
  // quick sanity checks
  CH_assert (a_fine_domain.checkPeriodic(fine_problem_domain));
  
  //pout() << "a_fine_domain.checkPeriodic(fine_problem_domain)) complete " << endl;
  
  if (a_coarse_domain.isClosed())
    {
                
      CH_assert (a_coarse_domain.checkPeriodic(a_crse_problem_domain));

      //
      // create the work array
      DisjointBoxLayout coarsened_fine_domain;
      coarsen ( coarsened_fine_domain,
                a_fine_domain,
                m_ref_ratio );

      const int coarse_slope_radius =
        (m_interp_radius + m_ref_ratio - 1) / m_ref_ratio;
      const int coarse_ghost_radius = coarse_slope_radius + s_stencil_radius;
      const IntVect coarse_slope = coarse_slope_radius * IntVect::Unit;
      // (wasteful) extra storage here, but who cares?
      for(int dir=0; dir<3; dir++)
        m_slopes[dir].define(coarsened_fine_domain,
                             a_num_comps,
                             coarse_slope);
      const IntVect coarse_ghost = coarse_ghost_radius * IntVect::Unit;
      m_coarsened_fine_data.define(coarsened_fine_domain,
                                   a_num_comps,
                                   coarse_ghost);
                                   
                                   
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        IntVect iv_off(IntVect::Zero);
        iv_off[dir]=1;
        
        Box dxBox(m_crse_problem_domain.domainBox().smallEnd()*iv_off, 
                  m_crse_problem_domain.domainBox().bigEnd()  *iv_off);
        dxBox.grow(dir,coarse_slope_radius);
        dxBox = dxBox & m_crse_problem_domain;
        
        m_dx[dir].define(dxBox,1);                     
        m_csh->dx(m_dx[dir],m_dx[dir].box(),dir,m_level);          
        
        Box crseCCbox( m_crse_problem_domain.domainBox().smallEnd()*iv_off,          
                       m_crse_problem_domain.domainBox().bigEnd()  *iv_off);
        crseCCbox.grow(dir,coarse_ghost_radius);
        crseCCbox = crseCCbox & m_crse_problem_domain;
                    
        m_crseCC[dir].define(crseCCbox,1);    
        m_csh->getCellCenters(m_crseCC[dir],m_crseCC[dir].box(),dir,m_level);
        
        Box fineCCbox( fine_problem_domain.domainBox().smallEnd()*iv_off,
                       fine_problem_domain.domainBox().bigEnd()  *iv_off);
        fineCCbox.grow(dir,m_interp_radius);
        fineCCbox = fineCCbox & fine_problem_domain;
                       
        m_fineCC[dir].define(fineCCbox,1);                
        m_csh->getCellCenters(m_fineCC[dir],m_fineCC[dir].box(),dir,m_level+1);        
      }      
      // Initialize aux data    
      for (int dir = SpaceDim; dir < 3; ++dir)
      {
        // Create fictions arrays
        m_dx[dir].define(Box(IntVect::Zero,IntVect::Zero), 1);
      }
      
      // Initialize uninitialized memory - a quick hack around the fact that
      // the ghost cells in coarsened_fine_domain may extend outside the
      // problem domain (PC, 7/21/00).
      {
        DataIterator dit = coarsened_fine_domain.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            m_coarsened_fine_data[dit()].setVal(-666.666);
          }
      }
      // allocate intvectsets
      m_fine_interp.define(a_fine_domain);
      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          m_coarse_centered_interp[dir].define(coarsened_fine_domain);
          m_coarse_lo_interp[dir].define(coarsened_fine_domain);
          m_coarse_hi_interp[dir].define(coarsened_fine_domain);
        }
      // Compute intvectsets. We do this mostly in the coarsened domain, with
      // only an intersection with the fine domain at the end.

      // first, create a box which will determine whether a given box
      // adjoins a periodic boundary
      Box periodicTestBox(m_crse_problem_domain.domainBox());
      if (m_crse_problem_domain.isPeriodic())
        {
          for (int idir=0; idir<SpaceDim; idir++)
            {
              if (m_crse_problem_domain.isPeriodic(idir))
                periodicTestBox.grow(idir,-1);
            }
        }

      // create regions in which to interpolate, then intersect with borders
      // to form regions for centered, one-sided low, and one-sided high
      // differences, per coordinate direction.

      DataIterator dit = coarsened_fine_domain.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          const Box& fine_box = a_fine_domain[dit()];
          Box coarsened_fine_box = m_crse_problem_domain
            & coarsen(grow(fine_box, m_interp_radius),m_ref_ratio);
          IntVectSet coarsened_fine_interp(coarsened_fine_box);

          // Iterate over boxes in coarsened fine domain, and subtract off
          // from the set of coarse cells from which the fine ghost cells
          // will be interpolated.

          LayoutIterator other_lit = coarsened_fine_domain.layoutIterator();
          for (other_lit.begin(); other_lit.ok(); ++other_lit)
            {
              const Box& other_coarsened_box
                = coarsened_fine_domain.get(other_lit());
              coarsened_fine_interp -= other_coarsened_box;
              // also need to remove periodic images from list of cells
              // to be filled, since they will be filled through exchange
              // as well
              if (m_crse_problem_domain.isPeriodic()
                  && !periodicTestBox.contains(other_coarsened_box)
                  && !periodicTestBox.contains(coarsened_fine_box))
                {
                  ShiftIterator shiftIt = m_crse_problem_domain.shiftIterator();
                  IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                  Box shiftedBox(other_coarsened_box);
                  for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                    {
                      IntVect shiftVect = shiftMult*shiftIt();
                      shiftedBox.shift(shiftVect);
                      coarsened_fine_interp -= shiftedBox;
                      shiftedBox.shift(-shiftVect);
                    }
                }
            }

          // Now that we have the coarsened cells required for interpolation,
          // construct IntvectSets specifying the one-sided and centered
          // stencil locations.

          for (int dir = 0; dir < SpaceDim; ++dir)
            {
              IntVectSet& coarse_centered_interp
                = m_coarse_centered_interp[dir][dit()];
              coarse_centered_interp = coarsened_fine_interp;

              IntVectSet& coarse_lo_interp = m_coarse_lo_interp[dir][dit()];
              coarse_lo_interp = coarse_centered_interp;
              coarse_lo_interp.shift(BASISV(dir));

              IntVectSet& coarse_hi_interp = m_coarse_hi_interp[dir][dit()];
              coarse_hi_interp = coarse_centered_interp;
              coarse_hi_interp.shift(-BASISV(dir));

              // We iterate over the coarse grids and subtract them off of the
              // one-sided stencils.
              LayoutIterator coarse_lit = a_coarse_domain.layoutIterator();
              for (coarse_lit.begin();coarse_lit.ok();++coarse_lit)
                {
                  Box bx = a_coarse_domain.get(coarse_lit());
                  coarse_lo_interp -= bx;
                  coarse_hi_interp -= bx;
                  // once again, need to do periodic images, too
                  if (m_crse_problem_domain.isPeriodic()
                      && !periodicTestBox.contains(bx)
                      && !periodicTestBox.contains(coarsened_fine_box))
                    {
                      ShiftIterator shiftIt = m_crse_problem_domain.shiftIterator();
                      IntVect shiftMult(m_crse_problem_domain.domainBox().size());
                      Box shiftedBox(bx);
                      for (shiftIt.begin(); shiftIt.ok(); ++shiftIt)
                        {
                          IntVect shiftVect = shiftMult*shiftIt();
                          shiftedBox.shift(shiftVect);
                          coarse_lo_interp -= shiftedBox;
                          coarse_hi_interp -= shiftedBox;
                          shiftedBox.shift(-shiftVect);
                        }
                    }
                }

              coarse_lo_interp.shift(-BASISV(dir));
              coarse_hi_interp.shift(BASISV(dir));
              coarse_centered_interp -= coarse_lo_interp;
              coarse_centered_interp -= coarse_hi_interp;

            }
          // Finally, we construct the fine cells that are going to be
          // interpolated by intersecting them with the refined version of the
          // coarse IntVectSet.

          IntVectSet& fine_interp = m_fine_interp[dit()];
          fine_interp = refine(coarsened_fine_interp,m_ref_ratio);
          fine_interp &= fine_problem_domain & grow(fine_box, m_interp_radius);

        }
      m_is_defined = true;
    } // end if coarser level is well-defined
}

bool
PiecewiseLinearFillPatchMHDAM::isDefined() const
{
  return ( m_is_defined );
}

// fill the interpolation region of the fine level ghost cells
void
PiecewiseLinearFillPatchMHDAM::fillInterp(
                                     LevelData<FArrayBox>& a_fine_data,
                                     const LevelData<FArrayBox>& a_old_coarse_data,
                                     const LevelData<FArrayBox>& a_new_coarse_data,
                                     Real a_time_interp_coef,
                                     int a_src_comp,
                                     int a_dest_comp,
                                     int a_num_comp
                                     )
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::fillInterp");
  // sanity checks
  CH_assert (m_is_defined);
  CH_assert (a_time_interp_coef >= 0.);
  CH_assert (a_time_interp_coef <= 1.);

  const DisjointBoxLayout oldCrseGrids = a_old_coarse_data.getBoxes();
  const DisjointBoxLayout newCrseGrids = a_new_coarse_data.getBoxes();
  const DisjointBoxLayout fineGrids = a_fine_data.getBoxes();

  CH_assert (oldCrseGrids.checkPeriodic(m_crse_problem_domain));
  CH_assert (newCrseGrids.checkPeriodic(m_crse_problem_domain));
  CH_assert (fineGrids.checkPeriodic(refine(m_crse_problem_domain,
                                           m_ref_ratio)));

  // time interpolation of coarse level data, to coarsened fine level work array
  timeInterp(a_old_coarse_data,
             a_new_coarse_data,
             a_time_interp_coef,
             a_src_comp,
             a_dest_comp,
             a_num_comp);
  
  // piecewise contant interpolation, from coarsened fine level work
  // array, to fine level
  fillConstantInterp(a_fine_data,
                     a_src_comp,
                     a_dest_comp,
                     a_num_comp);
  //
  // increment fine level data with per-direction linear terms
  computeSlopes( a_src_comp, a_num_comp );

  incrementLinearInterp(a_fine_data,
                        a_src_comp,
                        a_dest_comp,
                        a_num_comp);
}

// time interpolation of coarse level data, to coarsened fine level work array
void
PiecewiseLinearFillPatchMHDAM::timeInterp(
                                     const LevelData<FArrayBox>& a_old_coarse_data,
                                     const LevelData<FArrayBox>& a_new_coarse_data,
                                     Real a_time_interp_coef,
                                     int a_src_comp,
                                     int a_dest_comp,
                                     int a_num_comp
                                     )
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::timeInterp");
  Interval src_interval (a_src_comp,  a_src_comp  + a_num_comp - 1);
  Interval dest_interval(a_dest_comp, a_dest_comp + a_num_comp - 1);
  
  bool data_was_transformed = false;

  if ( (a_old_coarse_data.boxLayout().size() == 0)  &&
       (a_new_coarse_data.boxLayout().size() == 0) )
    {
      MayDay::Error ( "PiecewiseLinearFillPatchMHDAM::fillInterp: no old coarse data and no new coarse data" );
    }
  else if ( (a_time_interp_coef == 1.) ||
            (a_old_coarse_data.boxLayout().size() == 0) )
    {
      // old coarse data is absent, or fine time level is the new coarse time level
      a_new_coarse_data.copyTo( src_interval,
                                m_coarsened_fine_data,
                                dest_interval
                                );
    }
  else if ( (a_time_interp_coef == 0.) ||
            (a_new_coarse_data.boxLayout().size() == 0) )
    {
      // new coarse data is absent, or fine time level is the old coarse time level
      a_old_coarse_data.copyTo( src_interval,
                                m_coarsened_fine_data,
                                dest_interval
                                );
    }
  else
    {
      // linearly interpolate between old and new time levels
      a_new_coarse_data.copyTo( src_interval,
                                m_coarsened_fine_data,
                                dest_interval
                                );
      const DisjointBoxLayout&
        coarsened_fine_layout = m_coarsened_fine_data.disjointBoxLayout();
      LevelData<FArrayBox>
        tmp_coarsened_fine_data(coarsened_fine_layout,
                                m_coarsened_fine_data.nComp(),
                                m_coarsened_fine_data.ghostVect());

      // Initialize uninitialized memory - a quick hack around the fact that
      // the ghost cells in coarsened_fine_domain may extend outside the
      // problem domain (PC, 7/21/00).
      {
        DataIterator dit = coarsened_fine_layout.dataIterator();
        for (dit.begin(); dit.ok(); ++dit)
          {
            tmp_coarsened_fine_data[dit()].setVal(-666.666);
          }
      }
      a_old_coarse_data.copyTo( src_interval,
                                tmp_coarsened_fine_data,
                                dest_interval
                                );
      DataIterator dit = coarsened_fine_layout.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& coarsened_fine_fab = m_coarsened_fine_data[dit()];
          FArrayBox& tmp_coarsened_fine_fab = tmp_coarsened_fine_data[dit()];
          coarsened_fine_fab *= a_time_interp_coef;
          tmp_coarsened_fine_fab *= (1.- a_time_interp_coef);
          coarsened_fine_fab += tmp_coarsened_fine_fab;
          
          m_csh->PLFP_ConsToInterpVars(coarsened_fine_fab, coarsened_fine_fab.box(), m_level);
        }
        
      data_was_transformed = true;
    }
    
    if (data_was_transformed == false)
    {
      DataIterator dit = m_coarsened_fine_data.disjointBoxLayout().dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          FArrayBox& coarsened_fine_fab = m_coarsened_fine_data[dit()];                
          m_csh->PLFP_ConsToInterpVars(coarsened_fine_fab, coarsened_fine_fab.box(), m_level);
        }    
    }
}

// fill the fine interpolation region piecewise-constantly
void
PiecewiseLinearFillPatchMHDAM::fillConstantInterp(
                                             LevelData<FArrayBox>& a_fine_data,
                                             int a_src_comp,
                                             int a_dest_comp,
                                             int a_num_comp
                                             )
  const
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::fillConstantInterp");
  DataIterator dit = a_fine_data.boxLayout().dataIterator();
  
  Real tmp;
  int iFineLevel = m_level+1;

  for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& fine_fab = a_fine_data[dit()];
      const FArrayBox& coarse_fab = m_coarsened_fine_data[dit()];
      
      m_csh->PLFP_ConsToInterpVars(fine_fab, fine_fab.box(), iFineLevel);

#ifndef NDEBUG
      FORT_VIEWBOXDATA(
      CHF_FRA(fine_fab));
      
      FORT_VIEWBOXDATACONST(
      CHF_FRA(coarse_fab));
#endif
      
      const IntVectSet& local_fine_interp = m_fine_interp[dit()];
      IVSIterator ivsit(local_fine_interp);
      for (ivsit.begin(); ivsit.ok(); ++ivsit)
        {
          const IntVect& fine_iv = ivsit();
          IntVect coarse_iv = coarsen(fine_iv, m_ref_ratio);
          int coarse_comp = a_src_comp;
          int fine_comp   = a_dest_comp;
          for(; coarse_comp < a_src_comp + a_num_comp; ++fine_comp, ++coarse_comp)
          {
            tmp = coarse_fab(coarse_iv, coarse_comp);
            fine_fab(fine_iv, fine_comp) = tmp;
          }  

#ifndef NDEBUG          
          FORT_VIEWBOXDATA(
          CHF_FRA(fine_fab));
#endif          
        }
#ifndef NDEBUG        
      FORT_VIEWBOXDATA(
      CHF_FRA(fine_fab));  
#endif      
    }
}

// compute slopes at the coarse interpolation sites in the specified direction
void
PiecewiseLinearFillPatchMHDAM::computeSlopes(int a_src_comp,
                                        int a_num_comp)
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::computeSlopes");
  for (int dir=0; dir<3; dir++)
  {
    // Initialize uninitialized memory - a quick hack around the fact that
    // the ghost cells in coarsened_fine_domain may extend outside the
    // problem domain (PC, 7/21/00).
      DataIterator dit = m_slopes[0].boxLayout().dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {
          m_slopes[dir][dit()].setVal(-666.666);
        }
  }
  
  Real dxw, dx, dxe, c1, c2, DXP;
    
  
  DataIterator dit = m_coarsened_fine_data.boxLayout().dataIterator();
  for (int dir=0; dir<SpaceDim; dir++)
  {
      IntVect iv_off(IntVect::Zero);
      iv_off[dir]=1;
      
      //  for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
      //  {
      for (dit.begin(); dit.ok(); ++dit)
        {
          const FArrayBox& data_fab = m_coarsened_fine_data[dit()];
          #ifndef NDEBUG        
            FORT_VIEWBOXDATACONST(
             CHF_FRA(data_fab));  
          #endif      

          FArrayBox& slope_fab = m_slopes[dir][dit()];
          const IntVectSet& local_centered_interp = m_coarse_centered_interp[dir][dit()];
          const IntVectSet& local_lo_interp = m_coarse_lo_interp[dir][dit()];
          const IntVectSet& local_hi_interp = m_coarse_hi_interp[dir][dit()];
                    
          // van leer limited central difference
          IVSIterator centered_ivsit(local_centered_interp);
          if (m_csh->constStep(dir) == true)
          {
            c1 = 0.5/m_csh->dx(dir,m_level);
            for (centered_ivsit.begin(); centered_ivsit.ok(); ++centered_ivsit)
            {
              const IntVect& iv = centered_ivsit();
              const IntVect ivlo = iv - BASISV(dir);
              const IntVect ivhi = iv + BASISV(dir);                            
              
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
                {
                  Real dcenter = c1 * (data_fab(ivhi,comp) - data_fab(ivlo,comp));
                  slope_fab(iv,comp) = dcenter;
                }
            }
          } else
          {
            for (centered_ivsit.begin(); centered_ivsit.ok(); ++centered_ivsit)
            {            
              const IntVect& iv = centered_ivsit();
              const IntVect ivlo = iv - BASISV(dir);
              const IntVect ivhi = iv + BASISV(dir);
                                          
              dx  = m_dx[dir].get(iv*iv_off,0);
              dxe = m_dx[dir].get(ivhi*iv_off,0);
              dxw = m_dx[dir].get(ivlo*iv_off,0);
         
              c1  = 2.0*(0.5*dx+dxw)/(dx+dxe);
              c2  = 2.0*(0.5*dx+dxe)/(dx+dxw);
         
              DXP = 1.0/(dxw + dx + dxe);
        
              c1  = c1*DXP;
              c2  = c2*DXP;
              
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
                {            
                  Real dcenter = c1 * (data_fab(ivhi,comp) - data_fab(iv,  comp))+
                                 c2 * (data_fab(iv,  comp) - data_fab(ivlo,comp));
                                
                  slope_fab(iv,comp) = dcenter;
                }
            }
          } 
          //
          // one-sided difference (low)
          IVSIterator lo_ivsit(local_lo_interp);
          if (m_csh->constStep(dir) == true)
          {
            c1 = 1.0/m_csh->dx(dir,m_level);
            for (lo_ivsit.begin(); lo_ivsit.ok(); ++lo_ivsit)
            {
              const IntVect& iv = lo_ivsit();
              const IntVect ivlo = iv - BASISV(dir);                            
              
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
                {
                  Real dlo = c1*(data_fab(iv,comp) - data_fab(ivlo,comp));
                  slope_fab(iv,comp) = dlo;
                }
            }
          } else  
          {            
            for (lo_ivsit.begin(); lo_ivsit.ok(); ++lo_ivsit)
            {
              const IntVect& iv = lo_ivsit();
              const IntVect ivlo = iv - BASISV(dir);
                                          
              dx  = m_dx[dir].get(iv*iv_off,0);              
              dxw = m_dx[dir].get(ivlo*iv_off,0);
                            
              c1  = 2.0/(dx + dxw);       
              
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
                {
                  Real dlo = c1*(data_fab(iv,comp) - data_fab(ivlo,comp));
                  slope_fab(iv,comp) = dlo;
                }
            }
          }
          //
          // one-sided difference (high)
          IVSIterator hi_ivsit(local_hi_interp);
          if (m_csh->constStep(dir) == true)
          {
            c1 = 1.0/m_csh->dx(dir,m_level);
            for (hi_ivsit.begin(); hi_ivsit.ok(); ++hi_ivsit)
            {
              const IntVect& iv = hi_ivsit();
              const IntVect ivhi = iv + BASISV(dir);
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
                {
                  Real dhi = c1*(data_fab(ivhi,comp) - data_fab(iv,comp));
                  slope_fab(iv,comp) = dhi;
                }
            }
          } else 
          {
            for (hi_ivsit.begin(); hi_ivsit.ok(); ++hi_ivsit)
            {
              const IntVect& iv = hi_ivsit();
              const IntVect ivhi = iv + BASISV(dir);
                            
              
              dx  = m_dx[dir].get(iv*iv_off,0);
              dxe = m_dx[dir].get(ivhi*iv_off,0);                            
              c1  = 2.0/(dx + dxe);              
            
              for (int comp = a_src_comp; comp < a_src_comp + a_num_comp; ++comp)
              {
                Real dhi = c1*(data_fab(ivhi,comp) - data_fab(iv,comp));
                slope_fab(iv,comp) = dhi;
              }
            }
          }  
        } // end loop over boxes
        
    } // end loop over directions for simple slope computation
      
    
      
  
  // now do multidimensional slope limiting
  for (dit.begin(); dit.ok(); ++dit)
    {
      // this is the same stuff that is in FineInterp.cpp
      const Box& b = m_slopes[0][dit()].box();
      Box b_mod(b);
      b_mod.grow(1);
      b_mod = m_crse_problem_domain & b_mod;
      b_mod.grow(-1);

      // create a box big enough to remove periodic BCs from domain
      Box domBox = grow(b,2);
      domBox = m_crse_problem_domain & domBox;

      // to do limits, we need to have a box which includes
      // the neighbors of a given point (to check for the
      // local maximum...
      Box neighborBox(-1*IntVect::Unit,
                      IntVect::Unit);

      const FArrayBox& data_fab = m_coarsened_fine_data[dit()];
      
      FORT_INTERPLIMIT_V ( CHF_FRA(m_slopes[0][dit()]),
                   CHF_FRA(m_slopes[1][dit()]),
                   CHF_FRA(m_slopes[2][dit()]),
                   CHF_CONST_FRA(data_fab),
                   CHF_CONST_FRA1 ( m_dx[0], 0 ),
                   CHF_CONST_FRA1 ( m_dx[1], 0 ),
                   CHF_CONST_FRA1 ( m_dx[2], 0 ),
                   CHF_BOX ( b_mod ),
                   CHF_BOX ( neighborBox ),
                   CHF_BOX (domBox)                   
                   );        
            
    if( m_iBX >= a_src_comp && m_iBX + SpaceDim - 1 < a_src_comp + a_num_comp)
    {
      FORT_INTERP_DIVB_0( CHF_FRA      (m_slopes[0][dit()]),
                          CHF_FRA      (m_slopes[1][dit()]),
                          CHF_FRA      (m_slopes[2][dit()]),
                          CHF_BOX      (b_mod),
                          CHF_CONST_INT(m_iBX)
                        );
    }
    } // end loop over boxes for multiD slope limiting
}

// increment the fine interpolation sites with linear term for the
// specified coordinate direction
void
PiecewiseLinearFillPatchMHDAM::incrementLinearInterp(
                                                LevelData<FArrayBox>& a_fine_data,
                                                int a_src_comp,
                                                int a_dest_comp,
                                                int a_num_comp
                                                )
  const
{
  CH_TIME("PiecewiseLinearFillPatchMHDAM::incrementLinearInterp");
  
  int dir;
  Real refScale, coeff;  
  Real interp_coef;
      
  refScale = 1.0;
  for (dir=0; dir<SpaceDim; dir++) refScale*=m_ref_ratio;    
     
  int iFineLevel = m_level+1;    
  
  ProblemDomain fine_problem_domain = a_fine_data.disjointBoxLayout().physDomain();
     
                
  for (dir=0; dir<SpaceDim; dir++)
  {
    IntVect iv_off(IntVect::Zero);
    iv_off[dir]=1;
    
    DataIterator dit = a_fine_data.boxLayout().dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
      {
        const FArrayBox& slope_fab = m_slopes[dir][dit()];
        FArrayBox& fine_data_fab = a_fine_data[dit()];
        const IntVectSet& fine_interp = m_fine_interp[dit()];
        IVSIterator ivsit(fine_interp);
                        
#ifndef NDEBUG          
        FORT_VIEWBOXDATACONST(
        CHF_FRA(slope_fab));
        
        FORT_VIEWBOXDATA(
        CHF_FRA(fine_data_fab));
#endif          
        
        FArrayBox volc,volf;
        volc.define(m_crse_problem_domain & slope_fab.box()    ,1);
        volf.define(  fine_problem_domain & fine_data_fab.box(),1);
  
        m_csh->getCellVolumes(volc,volc.box(),m_level);
        m_csh->getCellVolumes(volf,volf.box(),iFineLevel);
        
        for (ivsit.begin(); ivsit.ok(); ++ivsit)
          {
            const IntVect& fine_iv = ivsit();
            const IntVect coarse_iv = coarsen(fine_iv,m_ref_ratio);            
            
            interp_coef = m_fineCC[dir].get(fine_iv*iv_off,0) -
                          m_crseCC[dir].get(coarse_iv*iv_off,0);
              
            int coarse_comp = a_src_comp;
            int fine_comp   = a_dest_comp;
                                          
            coeff = volc(coarse_iv,0)/(volf(fine_iv,0)*refScale);
                                    
            for (; coarse_comp < a_src_comp + a_num_comp; ++coarse_comp, ++fine_comp)
              {
                Real slope = coeff * interp_coef * slope_fab(coarse_iv,coarse_comp);                  
                //Real slope =  interp_coef * slope_fab(coarse_iv,coarse_comp);                  
                fine_data_fab(fine_iv,fine_comp) += slope;
              }
                         
              
            //fine_data_fab(fine_iv,0) = 1.0;
            //fine_data_fab(fine_iv,1) =  1.0*cos( (fine_iv[1]+0.5)*dphif);
            //fine_data_fab(fine_iv,2) = -1.0*sin( (fine_iv[1]+0.5)*dphif);
            //fine_data_fab(fine_iv,3) = 0.0;
            //fine_data_fab(fine_iv,4) = 0.154604375362396;
            //fine_data_fab(fine_iv,5) = 0.0;
            //fine_data_fab(fine_iv,6) = 0.0;
            //fine_data_fab(fine_iv,7) = 0.0;
          }
          
    #ifndef NDEBUG                
          FORT_VIEWBOXDATA(
        CHF_FRA(fine_data_fab));
    #endif    
      }
  }
  
  DataIterator dit = a_fine_data.boxLayout().dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {        
    FArrayBox& fine_data_fab = a_fine_data[dit()];           
    m_csh->PLFP_InterpVarsToCons(fine_data_fab, fine_data_fab.box(), iFineLevel);
  }
}

void
PiecewiseLinearFillPatchMHDAM::printIntVectSets() const
{
  DataIterator lit = m_fine_interp.boxLayout().dataIterator();
  int ibox = 0;
  for (lit.begin(); lit.ok(); ++lit)
    {
#if CHOMBO_VERSION_MAJOR < 4      
      cout << "grid " << lit().intCode() << ": " << endl;
#else
      cout << "grid " << ibox << ": " << endl;
#endif      
      cout << "fine ivs" << endl;
      cout << m_fine_interp[lit()] << endl;

      for (int dir = 0; dir < SpaceDim; ++dir)
        {
          cout << "coarse centered ivs [" << dir << "]: " << endl;
          cout << m_coarse_centered_interp[dir][lit()] << endl;
          cout << "coarse lo ivs [" << dir << "]: " << endl;
          cout << m_coarse_lo_interp[dir][lit()] << endl;
          cout << "coarse hi ivs [" << dir << "]: " << endl;
          cout << m_coarse_hi_interp[dir][lit()] << endl;
        }
      ibox++;
    }
}
#include "NamespaceFooter.H"
