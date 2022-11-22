#include "CSHCartesian.H"

#include <BoxIterator.H> 

#include "CSHCartesianF_F.H"
#include "DebugF_F.H"
#include "MHDAMDefs.H"
#include "EquationSystem.H"

                            
CSHCartesian::CSHCartesian()
 : CoordinateSystemHandler()
{
  m_dxdydz = NULL;
}

CSHCartesian::~CSHCartesian()
{
  if (m_dxdydz!=NULL) delete[] m_dxdydz;

}

void CSHCartesian::define(
            eCoordinateSystem  a_CS,
            int                a_max_level,
            const Vector<int>& a_ref_ratios,
            const Box&         a_prob_domain,
            RealVect           a_domainLength,
            const EquationSystem & a_eqSys)
   
{
  CoordinateSystemHandler::define(a_CS,a_max_level,a_ref_ratios,a_prob_domain, a_eqSys);
  
  m_max_level = a_max_level;
  m_dxdydz = new RealVect[m_max_level+1];
  
  int level,dim;
  
  level = 0;
  
  Real longestDim = MAX(a_domainLength[0],a_domainLength[1]);
  if (SpaceDim == 3) longestDim = MAX(longestDim,a_domainLength[2]);
  
  Real dx = longestDim / a_prob_domain.longside();
  m_dxdydz[level] = RealVect(D_DECL(dx, dx, dx));

  // Written for future use
  // for (dim = 0; dim <= SpaceDim; dim++) m_dxdydz[level][dim] = m_domainLength[dim] / m_prob_domain.size(dim);
  
  for (level = 1; level <= m_max_level; ++level)
  for (dim = 0; dim < SpaceDim; dim++) 
    m_dxdydz[level][dim] = m_dxdydz[level-1][dim] / m_ref_ratios[level-1];    
    
  m_loCurv = RealVect(D_DECL(0.0,0.0,0.0));
  m_hiCurv = RealVect(D_DECL(m_dxdydz[level][0]*m_prob_domain.size(0),
                             m_dxdydz[level][1]*m_prob_domain.size(1),
                             m_dxdydz[level][2]*m_prob_domain.size(2)));
}

void CSHCartesian::setGridSpacing(const FArrayBox & a_du,
                                              int   a_dir)
{
  MayDay::Error("CSHCartesian::setGridSpacing must not be used");  
}

bool CSHCartesian::constStep(int a_dir)
{
  return true;
}    

Real CSHCartesian::dx(int a_dir, int a_level)
{
  CH_assert(constStep(a_dir));
  
  return m_dxdydz[a_level][a_dir];
}

// Returns grid spacing. Works for all types of directions.
void CSHCartesian::dx(FArrayBox & a_dx,
                  const Box & a_box,
                  int         a_dir,
                  int         a_level)
{
  CH_assert(a_level <= m_max_level);
  
  a_dx.setVal(m_dxdydz[a_level][a_dir],a_box,0);
}

// Returns cell volumes  
void CSHCartesian::getCellVolumes(FArrayBox & a_vol,
                              const Box & a_box,
                              int         a_level)
{  
  CH_assert(a_level <= m_max_level);
  
  a_vol.setVal(D_TERM(m_dxdydz[a_level][0],*m_dxdydz[a_level][1],*m_dxdydz[a_level][2]) ,a_box,0);
}
  
// Get cell centers in a_dir direction  
void CSHCartesian::getCellCenters(FArrayBox & a_centers,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level)
{
  CH_assert(a_level <= m_max_level);
  
  Real dx = m_dxdydz[a_level][a_dir];
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    a_centers.set(iv, 0, dx*(iv[a_dir]+0.5));    
  }
  
#ifndef NDEBUG            
      FORT_VIEWBOXDATACONST(
      CHF_CONST_FRA(a_centers));
#endif

}

void CSHCartesian::getCellCenter( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level)
{
  D_TERM(
    a_coords[0] = m_dxdydz[a_level][0]*(a_iv[0]+0.5);,
    a_coords[1] = m_dxdydz[a_level][1]*(a_iv[1]+0.5);,
    a_coords[2] = m_dxdydz[a_level][2]*(a_iv[2]+0.5));
  
}

void CSHCartesian::getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level)
{
  CH_assert(a_level <= m_max_level);
  
  Real dx = m_dxdydz[a_level][a_dir];
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    a_coords.set(iv, 0, dx*iv[a_dir]);    
  }
  
#ifndef NDEBUG            
      FORT_VIEWBOXDATACONST(
      CHF_CONST_FRA(a_coords));
#endif
}

void CSHCartesian::getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level)
{
  CH_assert(a_level <= m_max_level);
    
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    D_TERM(
      a_coords.set(iv, 0, m_dxdydz[a_level][0]*iv[0]);,
      a_coords.set(iv, 1, m_dxdydz[a_level][1]*iv[1]);,
      a_coords.set(iv, 2, m_dxdydz[a_level][2]*iv[2]));    
  }
  
#ifndef NDEBUG            
      FORT_VIEWBOXDATACONST(
      CHF_CONST_FRA(a_coords));
#endif
}


// Get cartesian coordiantes of nodes
void CSHCartesian::getNodeCoordsCartesian( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level)
{
  CH_assert(a_level <= m_max_level);
  CH_assert(a_coords.box().contains(a_box));
    
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    D_TERM(
      a_coords.set(iv, 0, m_dxdydz[a_level][0]*iv[0]);,
      a_coords.set(iv, 1, m_dxdydz[a_level][1]*iv[1]);,
      a_coords.set(iv, 2, m_dxdydz[a_level][2]*iv[2]));    
  }
}


// Get cartesian coordiantes of a node
void CSHCartesian::getNodeCoordsCartesian( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level)
{  
  D_TERM(
    a_coords[0] = m_dxdydz[a_level][0]*a_iv[0];,
    a_coords[1] = m_dxdydz[a_level][1]*a_iv[1];,
    a_coords[2] = m_dxdydz[a_level][2]*a_iv[2]);            
}                              

// Get the closest index for a given curvilinear coordinate
void CSHCartesian::getClosestIndex(int  & a_i,
                           const Real & a_coord,
                           int          a_dir,                           
                           int          a_level)
{
  if ((a_dir>=0) && (a_dir<SpaceDim))
    a_i = floor(a_coord/m_dxdydz[a_level][a_dir]);
  else 
    MayDay::Error("CSHCartesian::getClosestIndex Wrong dimension");  
}                           

void CSHCartesian::getClosestCell(IntVect  & a_iv,
                            const RealVect & a_coords,
                            int              a_level)
{
  D_TERM(
    a_iv[0] = floor(a_coords[0]/m_dxdydz[a_level][0]);,
    a_iv[1] = floor(a_coords[1]/m_dxdydz[a_level][1]);,
    a_iv[2] = floor(a_coords[2]/m_dxdydz[a_level][2]));      
}

// Transform curvilinear coordiantes to cartesian
void CSHCartesian::transCurvCoordsToCartesian( RealVect & a_cs,
                              const RealVect & a_crv)
{
  a_cs = a_crv;
}

// Transform cartesian coordiantes to curvilinear
void CSHCartesian::transCartesianCoordsToCurv( RealVect & a_crv,
                              const RealVect & a_cs)
{
  a_crv = a_cs;
}

void CSHCartesian::getAreas(FArrayBox & a_areas,
                        const Box & a_box,
                        int         a_dir,
                        int         a_level)
{
  a_areas.setVal(1.0);
  //MayDay::Error("CSHCartesian::getAreas must not be used");  
}

// Transform fluxes calculated by Riemann problem class for using them in UpdateState
void CSHCartesian::transFluxesForUpdateState(                        
                        FArrayBox & a_F,                  
                  const Box       & a_box,
                  int               a_dir,
                  int               a_level)
{
  if (m_CoordinateSystem == CS_Cylindrical)
  {
    Real dy = m_dxdydz[a_level][1];
    FORT_TRANSFLUXES_CYL_UPDSTATE(
         CHF_FRA(a_F),         
         CHF_CONST_REAL(dy),     
         CHF_CONST_INT(a_dir),     
         CHF_BOX(a_box));
  }  
}

void CSHCartesian::multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
}

void CSHCartesian::multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const FArrayBox & a_areas,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
}

void CSHCartesian::scalingFactor(
                        FArrayBox & a_scale, 
                  const Box       & a_box,                        
                        int         a_level)
{
  CH_assert(a_level <= m_max_level);
  
 if (m_CoordinateSystem == CS_Cylindrical)
  {
    Real dy = m_dxdydz[a_level][1];
    FORT_SCALINGFACTOR_CYL(
         CHF_FRA1(a_scale,0),         
         CHF_CONST_REAL(dy),              
         CHF_BOX(a_box));
  } else 
  {
    Real dx = m_dxdydz[a_level][0];
    a_scale.setVal(1.0/dx);
  }
}

// square root of orthogonal components (g_ii) of the metric tensor
void CSHCartesian::sqrtMetricTensorCoeff(
                        FArrayBox & a_g,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
  a_g.setVal(1.0);
}