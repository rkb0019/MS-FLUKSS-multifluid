#include "FORT_PROTO.H"
#include "CH_HDF5.H"
#include <BoxIterator.H>

#include "CSHSpherical.H"
#include "MHDAMDefs.H"


#include "CSHSphericalF_F.H"


extern "C"
{
#define FORT_INITSPHERICAL FORTRAN_NAME( INITSPHERICAL ,initspherical )
void FORT_INITSPHERICAL
(
  int*  const nlevels, 
  Real* const lo, 
  Real* const hi,   
  int*  const sizes,
  int*  const ref_ratios,
  int*  const rspacing, 
  Real* const expCoeff
);

#define FORT_SETGRIDSPACING FORTRAN_NAME( SETGRIDSPACING ,setgridspacing )
void FORT_SETGRIDSPACING
(
  const Real *const a_dr, 
  const int  *const a_dir
);


#define FORT_FINALIZESPHERICAL FORTRAN_NAME( FINALIZEPOLAR ,finalizespherical )
void FORT_FINALIZESPHERICAL();
}


                            
CSHSpherical::CSHSpherical()
{
  m_gridData = NULL;
  m_Rspacing = R_Const;  
  m_expAlpha = 1.0;
  D_TERM(
     m_constStep[0]=true;,
     m_constStep[1]=true;,
     m_constStep[2]=true;);
}

CSHSpherical::~CSHSpherical()
{  
  FORT_FINALIZESPHERICAL(); 
  if (m_gridData!=NULL) delete[] m_gridData;
}

void CSHSpherical::define(
            eCoordinateSystem   a_CS,
            int                 a_max_level,
            const Vector<int>&  a_ref_ratios,
            const Box&          a_prob_domain,
            RealVect            a_domainLength,
            const Vector<Real>& a_domainBox,
            const EquationSystem & a_eqSys
            )
   
{
  CoordinateSystemHandler::define(a_CS,a_max_level,a_ref_ratios,a_prob_domain, a_eqSys);
  m_gridData = new SphericalLevel[m_max_level+1];
  
  int level,dim,i;
  int sizes[CH_SPACEDIM]={D_DECL(m_prob_domain.size(0),m_prob_domain.size(1),m_prob_domain.size(2))};
  
  level = 0;
  for (dim = 0; dim < SpaceDim; dim++) m_gridData[level].m_dudvdw_const[dim] = a_domainLength[dim] / m_prob_domain.size(dim);
  
  for (level = 1; level <= m_max_level; ++level)
  for (dim   = 0; dim<SpaceDim; dim++) 
    m_gridData[level].m_dudvdw_const[dim] = m_gridData[level-1].m_dudvdw_const[dim] / m_ref_ratios[level-1];  
    
  int rsize   = m_prob_domain.size(0);
  int phisize = m_prob_domain.size(1);
  
  CH_assert(m_max_level+1<=m_ref_ratios.size());
  int*  refRatios2 =  new int[m_max_level+1];
  for (i=0; i<=m_max_level; i++) refRatios2[i] = m_ref_ratios[i];
    
  
  m_loCurv = RealVect(D_DECL(a_domainBox[0],a_domainBox[1],a_domainBox[2]));
  m_hiCurv = RealVect(D_DECL(a_domainBox[SpaceDim],a_domainBox[SpaceDim+1],a_domainBox[SpaceDim+2]));
  
  int Rspacing = (int)m_Rspacing;  
  
  if (m_Rspacing!=R_Const) m_constStep[0] = false;
  
  FORT_INITSPHERICAL(&m_max_level, m_loCurv.dataPtr(), m_hiCurv.dataPtr(), sizes, refRatios2, &Rspacing, &m_expAlpha);
  
  delete[] refRatios2;
  

  
  // Allocating m_gridData
  
  Box boxD = m_prob_domain;
  
  for (int level = 0; level <= m_max_level; ++level)
  {
    if (level > 0)
      boxD.refine(m_ref_ratios[level-1]);
      
    m_gridData[level].m_prob_domain = boxD;
      
  }
         
  for (int dir = 0; dir < SpaceDim; ++dir)
  {
    IntVect iv_off (IntVect::Zero);
    iv_off[dir] = 1;
        
    boxD.define(m_prob_domain.smallEnd()*iv_off, 
                m_prob_domain.bigEnd()  *iv_off);
              
    Box b;    
    for (int level = 0; level <= m_max_level; ++level)
    {
      if (level > 0)
        boxD.refine(m_ref_ratios[level-1]);
        
      boxD.define(boxD.smallEnd()*iv_off, 
                  boxD.bigEnd()  *iv_off);
      
      Box b = boxD;         
      b.grow(dir,4);      // four guest cells
      m_gridData[level].m_dudvdw[dir].define(b,1);
      m_gridData[level].m_uvwc[dir].define(b,1);    
    
      b = boxD;
      b.surroundingNodes(dir);
      b.grow(dir,4);
      m_gridData[level].m_uvwn[dir].define(b,1);
      
      FORT_SETLEVELGRIDSPACING(
        CHF_FRA1(m_gridData[level].m_dudvdw[dir],0),
        CHF_FRA1(m_gridData[level].m_uvwc[dir],0),
        CHF_FRA1(m_gridData[level].m_uvwn[dir],0),
        CHF_CONST_INT(dir),
        CHF_CONST_INT(level));

    }
  }
    
    
}

// Input parameters
void CSHSpherical::input( ParmParse & parser, int verbosity )
{
  CoordinateSystemHandler::input( parser, verbosity );
  
  int Rspacing = 0;
  parser.query("Rspacing",Rspacing);
  m_Rspacing = (eRspacing)(Rspacing);
  
  if (m_Rspacing == R_Exp)
    parser.query("Rexpalpha",m_expAlpha);  
  
}

void CSHSpherical::setGridSpacing(const FArrayBox & a_du, int a_dir)
{
  CH_assert(a_du.box().size((a_dir+1)%SpaceDim) == 1);
  if (SpaceDim == 3) CH_assert(a_du.box().size((a_dir+2)%SpaceDim) == 1);
  CH_assert(a_du.box().smallEnd(a_dir) == m_prob_domain.smallEnd(a_dir));
  CH_assert(a_du.box().bigEnd(a_dir)   == m_prob_domain.bigEnd(a_dir));  
  
  
  FORT_SETGRIDSPACING(a_du.dataPtr(),&a_dir);      
  
  for (int level = 0; level <= m_max_level; ++level)
  {
    FORT_SETLEVELGRIDSPACING(
      CHF_FRA1(m_gridData[level].m_dudvdw[a_dir],0),
      CHF_FRA1(m_gridData[level].m_uvwc[a_dir],0),
      CHF_FRA1(m_gridData[level].m_uvwn[a_dir],0),
      CHF_CONST_INT(a_dir),
      CHF_CONST_INT(level));
  }    
  
  m_constStep[a_dir]=false;
        
  
}



  
// Returns true if grid spacing is constant in this direction  
bool CSHSpherical::constStep(int a_dir)
{
  //if (a_dir == 0) return false; else return m_constStep[a_dir];
  
  return m_constStep[a_dir];
}


// Returns grid spacing if it is constant in this direction  
Real CSHSpherical::dx(int a_dir, int a_level) 
{
  CH_assert(m_constStep[a_dir]==true);
  
  return m_gridData[a_level].m_dudvdw_const[a_dir];
  
}



// Returns grid spacing. Works for all types of directions.

void CSHSpherical::dx(FArrayBox & a_dx,
                const Box & a_box,
                int         a_dir,
                int         a_level)
{
  a_dx.copy(m_gridData[a_level].m_dudvdw[a_dir],a_box);           
}

// Returns cell volumes  
void CSHSpherical::getCellVolumes(FArrayBox & a_vol,
                            const Box & a_box,
                            int         a_level)
{
  FORT_VOLUMESPHERICAL(
       CHF_FRA1(a_vol,0),     
       CHF_BOX(a_box),     
       CHF_CONST_INT(a_level));

}                            

// Get cell centers in a_dir direction  
void CSHSpherical::getCellCenters(FArrayBox & a_centers,
                            const Box & a_box,
                            int         a_dir,
                            int         a_level)
{
  CH_assert(m_gridData[a_level].m_uvwc[a_dir].box().contains(a_box));
  a_centers.copy(m_gridData[a_level].m_uvwc[a_dir],a_box);             
}

void CSHSpherical::getCellCenter( RealVect  & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level)
{
  IntVect iv_off;    
  D_TERM(    
    iv_off = IntVect::Zero; iv_off[0] = 1; a_coords[0] = m_gridData[a_level].m_uvwc[0].get(a_iv*iv_off,0);,
    iv_off = IntVect::Zero; iv_off[1] = 1; a_coords[1] = m_gridData[a_level].m_uvwc[1].get(a_iv*iv_off,0);,
    iv_off = IntVect::Zero; iv_off[2] = 1; a_coords[2] = m_gridData[a_level].m_uvwc[2].get(a_iv*iv_off,0););

}


void CSHSpherical::getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level)
{
  a_coords.copy(m_gridData[a_level].m_uvwn[a_dir],a_box);             
}

void CSHSpherical::getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,
                              int         a_level)
{
  BoxIterator bit( a_box );  
  IntVect iv; Real D_DECL(r,phi,theta); 
 
  for( bit.begin(); bit.ok(); ++bit )
  {
    iv = bit();
    D_TERM(
      r     = m_gridData[a_level].m_uvwn[0].get(iv*BASISV(0),0);,
      phi   = m_gridData[a_level].m_uvwn[1].get(iv*BASISV(1),0);,
      theta = m_gridData[a_level].m_uvwn[2].get(iv*BASISV(2),0));
      
    D_TERM(
      a_coords.set(iv,0,r);,
      a_coords.set(iv,1,phi);,
      a_coords.set(iv,2,theta););
  }   
}


// Get cartesian coordiantes of nodes
void CSHSpherical::getNodeCoordsCartesian(FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level)
{
  FORT_GETNODECOORDSSPHERICAL_CS_FAB(
    CHF_FRA(a_coords),
    CHF_BOX(a_box),    
    CHF_CONST_INT(a_level));
    
}

// Get cartesian coordiantes of a node
void CSHSpherical::getNodeCoordsCartesian( RealVect & a_coords,
                                      const IntVect & a_iv,                              
                                      int             a_level)
{
  FORT_GETNODECOORDSSPHERICAL_CS(
        CHF_REALVECT(a_coords),
        CHF_CONST_INTVECT(a_iv),
        CHF_CONST_INT(a_level));
}

void CSHSpherical::getClosestCell(IntVect  & a_iv,
                            const RealVect & a_coords,
                            int              a_level)
{  
  IntVect iv1,iv2; Real a1,a2; int dir;
  for (dir = 0; dir < SpaceDim; dir++)
  {
    a_iv[dir]=-1;
    
    IntVect iv_off (IntVect::Zero);
    iv_off[dir] = 1;
    
    BoxIterator bit( m_gridData[a_level].m_uvwc[dir].box() );    
    for( bit.begin(); bit.ok(); ++bit )
    {
      iv1 = bit();
      iv2 = iv1 + iv_off*IntVect::Unit;
      
      a1  = m_gridData[a_level].m_uvwn[dir].get(iv1,0);
      a2  = m_gridData[a_level].m_uvwn[dir].get(iv2,0);
      if ((a1<=a_coords[dir])&&(a_coords[dir]<=a2))
      {
        a_iv[dir]=iv1[dir];
        break;
      }
    }
  }  
}

// Get the closest index for a given curvilinear coordinate
void CSHSpherical::getClosestIndex(int  & a_i,
                           const Real & a_coord,
                           int          a_dir,                           
                           int          a_level)
{
  IntVect iv1,iv2; Real a1,a2;
  
  a_i=-1;
  
  IntVect iv_off (IntVect::Zero);
  iv_off[a_dir] = 1;
  
  BoxIterator bit( m_gridData[a_level].m_uvwc[a_dir].box() );    
  for( bit.begin(); bit.ok(); ++bit )
  {
    iv1 = bit();
    iv2 = iv1 + iv_off*IntVect::Unit;
    
    a1  = m_gridData[a_level].m_uvwn[a_dir].get(iv1,0);
    a2  = m_gridData[a_level].m_uvwn[a_dir].get(iv2,0);
    if ((a1<=a_coord)&&(a_coord<=a2))
    {
      a_i=iv1[a_dir];
      break;
    }
  }
    
}

// Transform curvilinear coordiantes to cartesian
void CSHSpherical::transCurvCoordsToCartesian( RealVect & a_cs,
                                         const RealVect & a_crv)
{
#if CH_SPACEDIM == 2
  a_cs[0]=a_crv[0]*cos(a_crv[1]);
  a_cs[1]=a_crv[0]*sin(a_crv[1]);  
#endif

#if CH_SPACEDIM == 3 
  a_cs[0]=a_crv[0]*sin(a_crv[2])*cos(a_crv[1]);
  a_cs[1]=a_crv[0]*sin(a_crv[2])*sin(a_crv[1]);
  a_cs[2]=a_crv[0]*cos(a_crv[2]);
#endif
}

// Transform cartesian coordiantes to curvilinear
void CSHSpherical::transCartesianCoordsToCurv( RealVect & a_crv,
                              const RealVect & a_cs)
{
  a_crv[0] = sqrt(D_TERM(a_cs[0]*a_cs[0],+a_cs[1]*a_cs[1],+a_cs[2]*a_cs[2]));
  
  Real r   = MAX(1e-8,sqrt(a_cs[0]*a_cs[0]+a_cs[1]*a_cs[1]));  
  Real phi = acos(a_cs[0]/r);
  if (a_cs[1]<0.0) phi = d_2PI-phi;
  
  a_crv[1] = phi;
  
#if CH_SPACEDIM == 3 
  r        = MAX(1e-8,sqrt(a_cs[0]*a_cs[0]+a_cs[1]*a_cs[1]+a_cs[2]*a_cs[2]));  
  a_crv[2] = acos(a_cs[2]/r);  
#endif

}


// Get cells areas
void CSHSpherical::getAreas(FArrayBox & a_areas,
                      const Box & a_box,
                      int         a_dir,
                      int         a_level)
{
  if (CH_SPACEDIM == 2)
  {
    FORT_AREASPOLAR(
      CHF_FRA1(a_areas,0),
      CHF_BOX(a_box),
      CHF_CONST_INT(a_dir),
      CHF_CONST_INT(a_level));
  }
  if (CH_SPACEDIM == 3)
  {
    FORT_AREASSPHERICAL(
      CHF_FRA1(a_areas,0),
      CHF_BOX(a_box),
      CHF_CONST_INT(a_dir),
      CHF_CONST_INT(a_level));
  }
  
}

// Modifies conservative variables in a_interp to  values
// that will be used PiecewiseLinearFillPatchMHDAM for filling ghost cells
void CSHSpherical::PLFP_ConsToInterpVars(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level)
{          
  FORT_COMPUTE_VELOCITIES(
        CHF_FRA(a_interp),   
        CHF_CONST_I1D(m_velVectors,m_nVelVectors),     
        CHF_BOX(a_box));  
  
  transCartesianVectToCurv(a_interp,m_vectors,m_nVectors,a_box,a_level);         
}


    
// Reverse operation to PLFP_ConsToInterpVars
void CSHSpherical::PLFP_InterpVarsToCons(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level)
{
    FORT_COMPUTE_MOMENTUM(
        CHF_FRA(a_interp),   
        CHF_CONST_I1D(m_velVectors,m_nVelVectors),     
        CHF_BOX(a_box));  
    
  transCurvVectToCartesian(a_interp,m_vectors,m_nVectors,a_box,a_level);  
}

// Transform Cartesian components of a vector into curvilinear
void CSHSpherical::transCartesianVectToCurv(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level)
{    
  IndexType bType =  a_U.box().ixType();  
  CH_assert( D_TERM(
       (bType.ixType(0) == IndexType::CELL),
    && (bType.ixType(1) == IndexType::CELL),
    && (bType.ixType(2) == IndexType::CELL) ) );
    
      
  if (SpaceDim == 2)
  {
    FORT_TRANSCARTESIANTOPOLAR(
      CHF_FRA(a_U), 
      CHF_CONST_I1D(a_vectors,a_numvectors),
      CHF_BOX(a_box),
      CHF_CONST_INT(a_level) );                  
  }
  if (SpaceDim == 3)
  {
    Box & pd  = m_gridData[a_level].m_prob_domain;
    
    /*Box boxC = a_box;
    if (boxC.smallEnd(2)<pd.smallEnd(2))
    {
      Box boxB = a_box;
      boxB.setBig(2,pd.smallEnd(2)-1);
      boxC.setSmall(2,pd.smallEnd(2));
      
      FORT_TRANSCARTESIANTOSPHERICALZ(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_BOX(boxB),
        CHF_CONST_INT(a_level) );                        
    }
    
    if (boxC.bigEnd(2)>pd.bigEnd(2))
    {
      Box boxB = a_box;
      boxB.setSmall(2,pd.bigEnd(2)+1);
      boxC.setBig(2,pd.bigEnd(2));
      
      FORT_TRANSCARTESIANTOSPHERICALZ(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_BOX(boxB),
        CHF_CONST_INT(a_level) );                        
    }*/
            
    FORT_TRANSCARTESIANTOSPHERICAL(
      CHF_FRA(a_U), 
      CHF_CONST_I1D(a_vectors,a_numvectors),
      CHF_BOX(a_box),
      CHF_CONST_INT(a_level) );                      
            
  }
  
  
}
                  
  // Transform curvilinear components of a vector into Cartesian
void CSHSpherical::transCurvVectToCartesian(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level)
{  
    
  IndexType bType =  a_U.box().ixType();  
  if ( D_TERM(
       (bType.ixType(0) == IndexType::CELL),
    && (bType.ixType(1) == IndexType::CELL),
    && (bType.ixType(2) == IndexType::CELL)) )
  {    
    int iType = 0;
    if (SpaceDim == 2)
    {
      FORT_TRANSPOLARTOCARTESIAN(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_CONST_INT(iType),
        CHF_BOX(a_box),
        CHF_CONST_INT(a_level) );                  
    }
    if (SpaceDim == 3)
    {
      FORT_TRANSSPHERICALTOCARTESIAN(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_CONST_INT(iType),
        CHF_BOX(a_box),
        CHF_CONST_INT(a_level) );                  
    }
    return;
  }
  
  int iFace = -1;
  if (D_TERM(
       (bType.ixType(0) == IndexType::NODE),
    && (bType.ixType(1) == IndexType::CELL),
    && (bType.ixType(2) == IndexType::CELL)) )
    iFace = 0;
  else if (D_TERM(
       (bType.ixType(0) == IndexType::CELL),
    && (bType.ixType(1) == IndexType::NODE),
    && (bType.ixType(2) == IndexType::CELL)) )
    iFace = 1;
  else if ( (SpaceDim == 3) && D_TERM(
       (bType.ixType(0) == IndexType::CELL),
    && (bType.ixType(1) == IndexType::CELL),
    && (bType.ixType(2) == IndexType::NODE)) )
    iFace = 2;
    
  if (iFace == -1) MayDay::Error("CSHSpherical::transCurvVectToCartesian doesn't support required type of box");
  
  if (SpaceDim == 2)
  {    
    FORT_TRANSPOLARTOCARTESIAN(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_CONST_INT(iFace),
        CHF_BOX(a_box),
        CHF_CONST_INT(a_level) );   
  }  
  if (SpaceDim == 3)
  {   
      FORT_TRANSSPHERICALTOCARTESIAN(
        CHF_FRA(a_U), 
        CHF_CONST_I1D(a_vectors,a_numvectors),
        CHF_CONST_INT(iFace),
        CHF_BOX(a_box),
        CHF_CONST_INT(a_level) );    
  }
  
}

// Transform fluxes calculated by Riemann problem class for using them in UpdateState
void CSHSpherical::transFluxesForUpdateState(                        
                        FArrayBox & a_F,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{    
  transCurvVectToCartesian(a_F,m_vectors,m_nVectors, a_box, a_level);    
}

void CSHSpherical::scalingFactor(
                        FArrayBox & a_scale, 
                  const Box       & a_box,                        
                        int         a_level)
{
  getCellVolumes(a_scale, a_box, a_level);
  a_scale.invert(1.0);  
}

bool CSHSpherical::scaleFineFluxes()
{
  return false;
}


// Calculates fluxes for using them in FluxRegistry.incrementCoarse
void CSHSpherical::modifyFluxesFRincrementCoarse(
                        FArrayBox & fluxC, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level)
{
  fluxC.copy(flux);
  transFluxesForUpdateState(fluxC,fluxC.box(),idir,a_level);  
}                  

  // Calculates fluxes for using them in FluxRegistry.incrementFine
void CSHSpherical::modifyFluxesFRincrementFine(
                        FArrayBox & fluxF, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level)
{
  fluxF.copy(flux);  
  int refRatio = m_ref_ratios[a_level-1];   
  
  if (SpaceDim == 3)
  {
    FORT_FR_INCREMENTFINE_SPHERICAL(
          CHF_FRA(fluxF), 
          CHF_CONST_I1D(m_vectors,m_nVectors),
          CHF_CONST_INT(refRatio),
          CHF_BOX(fluxF.box()),
          CHF_CONST_INT(idir),
          CHF_CONST_INT(a_level) );
  }
  if (SpaceDim == 2)
  {
    FORT_FR_INCREMENTFINE_POLAR(
          CHF_FRA(fluxF), 
          CHF_CONST_I1D(m_vectors,m_nVectors),
          CHF_CONST_INT(refRatio),
          CHF_BOX(fluxF.box()),
          CHF_CONST_INT(idir),
          CHF_CONST_INT(a_level) );
  }
}

void CSHSpherical::sqrtMetricTensorCoeff(
                        FArrayBox & a_g,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
  if (a_dir == 0)
  {
    a_g.setVal(1.0);
    return;
  }
  if ( ((SpaceDim == 2) && (a_dir == 1)) ||
       ((SpaceDim == 3) && (a_dir == 2)) )
  {
    // sqrt(g) = r
    FORT_METRICTENSOR_R(
          CHF_FRA1(a_g,0),           
          CHF_BOX(a_box),          
          CHF_CONST_INT(a_level));    
  } else
  {    
    CH_assert((SpaceDim == 3) && (a_dir == 1));    
    
    // sqrt(g) = r*sin(theta)
    FORT_METRICTENSOR11_SPHERICAL(
          CHF_FRA1(a_g,0),           
          CHF_BOX(a_box),          
          CHF_CONST_INT(a_level));    
  }
  
}
                  
// Boundary conditions through the z-axis. 
void CSHSpherical::zaxisBC(
                       FArrayBox  & a_W,                 
                 const Box        & a_box,
                       int          a_sign)
{  
  for (int i=0; i<m_nVectors; i++)
  {
    a_W.mult(-1.0, a_box, m_vectors[i]+1, 2);
  }
  
}

#ifdef CH_USE_HDF5
// Reads coordinate specific data  
void CSHSpherical::readGeomInfo(hid_t a_fileID, bool a_dirs[CH_SPACEDIM])
{    
  //a_handle.setGroup("/");
  //a_handle.setGroup("geometry");

#ifdef H516
  hid_t groupID = H5Gopen(a_fileID, "/geometry");
#else
  hid_t groupID = H5Gopen2(a_fileID, "/geometry", H5P_DEFAULT);
#endif

  if ((groupID<0) && (procID() == 0))  
  {
    MayDay::Warning("CSHSpherical::readGeomInfo Geometry information is missing in a checkpoint file");
    return;
  }
  
  std::string dataspace_names[3];
  dataspace_names[0] = "dr";
  dataspace_names[1] = "dphi";
  dataspace_names[2] = "dtheta";

  for (int dir = 0;dir < SpaceDim; ++dir)
  if (a_dirs[dir]==true)
  { 
    IntVect iv_off(IntVect::Zero);
    iv_off[dir]=1;  
    hsize_t dimsf = 0;
    
    //H5E_auto_t efunc; void* edata;

#ifdef H516
    H5E_auto_t efunc; void* edata;
    H5Eget_auto(&efunc, &edata);      
    H5Eset_auto(NULL, NULL);
    hid_t dataset = H5Dopen(groupID, dataspace_names[dir].c_str());  
    H5Eset_auto(efunc, edata);
#else
    H5E_auto2_t efunc; void* edata;
    H5Eget_auto2(H5E_DEFAULT,&efunc, &edata); //added 2 here
    H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
    hid_t dataset = H5Gopen2(groupID, dataspace_names[dir].c_str(),H5P_DEFAULT);  
    H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
    
    if (dataset < 0) continue;
    
    //std::string err(std::string(dataspace_names[dir]));
    
    hid_t dataspace = H5Dget_space(dataset);
    if (H5Sget_simple_extent_ndims(dataspace)!=1)
    {
      std::string err(dataspace_names[dir]+std::string(" array should be 1D array in checkpoint file"));
      MayDay::Error(err.c_str());   
    }
    H5Sget_simple_extent_dims(dataspace, &dimsf, NULL);
    if (dimsf!=m_prob_domain.size(dir))
    {
      std::string err("CSHSpherical::readGeomInfo. Size of the "+std::string(dataspace_names[dir])+" array in geometry section does not agree with the size defined in inputs file");
      MayDay::Error(err.c_str());   
    }
      
    FArrayBox buf(Box(m_prob_domain.smallEnd()*iv_off, 
                      m_prob_domain.bigEnd()  *iv_off),1);    
                  
    herr_t status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.dataPtr());
      
    H5Sclose(dataspace);
    H5Dclose(dataset); 

    setGridSpacing(buf, dir);    
  }
  //a_handle.setGroup("/");
  H5Gclose(groupID);
}
  
// Writes coordinate specific data
void CSHSpherical::writeGeomInfo(hid_t a_fileID)
{
  //a_handle.setGroup("/");
  //a_handle.setGroup("geometry");
  
#ifdef H516
  hid_t groupID = H5Gcreate(a_fileID, "/geometry", 0);
#else
  hid_t groupID = H5Gcreate2(a_fileID, "/geometry", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif  

  //hid_t aid  = H5Screate(H5S_SCALAR);                                 
  //hid_t attr = H5Acreate(a_handle.groupID(), "SpaceDim", H5T_NATIVE_INT, aid, H5P_DEFAULT); 
  //H5Awrite(attr, H5T_NATIVE_INT, &SpaceDim);   
  //H5Aclose(attr);                                                     
  //H5Sclose(aid);                                                     
  
  Real phys_domain_buf[2*CH_SPACEDIM] = {
    D_DECL(m_loCurv[0],m_loCurv[1],m_loCurv[2]),
    D_DECL(m_hiCurv[0],m_hiCurv[1],m_hiCurv[2]) };      
  hsize_t dimsf = 2*CH_SPACEDIM;
  hid_t dataspace = H5Screate_simple(1, &dimsf, NULL);  

#ifdef H516        
  hid_t attr = H5Acreate(groupID, "phys_domain", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
  hid_t attr = H5Acreate2(groupID, "phys_domain", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#endif
 
  H5Awrite(attr, H5T_NATIVE_DOUBLE, phys_domain_buf);   
  H5Aclose(attr);                                                     
  H5Sclose(dataspace);                                                     
  
  int step_const[CH_SPACEDIM] = {  D_DECL(
    m_constStep[0] ? 1 : 0,
    m_constStep[1] ? 1 : 0,
    m_constStep[2] ? 1 : 0) };    
  dimsf = CH_SPACEDIM;
  dataspace = H5Screate_simple(1, &dimsf, NULL);   

#ifdef H516       
  attr = H5Acreate(groupID, "step_const", H5T_NATIVE_INT, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(groupID, "step_const", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); 
#endif

  H5Awrite(attr, H5T_NATIVE_INT, step_const);   
  H5Aclose(attr);                                                     
  H5Sclose(dataspace);                     
  
  char const * dataspace_names[3] = {"dr","dphi","dtheta"};

  for (int dir = 0;dir < SpaceDim; ++dir)
  if (m_constStep[dir] == false)
  { 
    IntVect iv_off(IntVect::Zero);
    iv_off[dir]=1;  
    Box b(m_prob_domain.smallEnd()*iv_off, 
          m_prob_domain.bigEnd()  *iv_off);  
              
    FArrayBox spacing(b,1);                     
    dx(spacing, b, dir, 0);
    
    //Real * dr_buffer = new Real[drBox.size(0)];    
    //int i;
    //for (i=0;i<drBox.size(0);i++)
    //{
    //  dr_buffer[i]=dr.get(IntVect(D_DECL(i,0,0)),0);
    //}
    
    //a_handle.setGroup("/");
    //a_handle.setGroup("geometry");
    
    dimsf = b.size(dir);
    
    dataspace       = H5Screate_simple(1, &dimsf, NULL); 

#ifdef H516  
    hid_t dataset   = H5Dcreate(groupID, dataspace_names[dir], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);  
#else
    hid_t dataset   = H5Dcreate2(groupID, dataspace_names[dir], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
#endif

    herr_t status   = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, spacing.dataPtr());
      
    H5Sclose(dataspace);
    H5Dclose(dataset);    
        
    //delete[] dr_buffer;
  }
  //a_handle.setGroup("/");
  
  H5Gclose(groupID);

}
#endif

