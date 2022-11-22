#include "CSHandler.H"
#include "CSHandlerF_F.H"
#include "EquationSystem.H"
  

CoordinateSystemHandler::CoordinateSystemHandler()
{
  m_vectors = NULL;
  m_velVectors = NULL;
}

CoordinateSystemHandler::~CoordinateSystemHandler()
{
  if (m_vectors!=NULL) delete[] m_vectors;
  if (m_velVectors!=NULL) delete[] m_velVectors;
}

void CoordinateSystemHandler::define(
            eCoordinateSystem   a_CS,
            int                 a_max_level,
            const Vector<int> & a_ref_ratios,
            const Box         & a_prob_domain,
            const EquationSystem & a_eqSys)
{
  m_CoordinateSystem          = a_CS;
  m_max_level   = a_max_level;
  m_ref_ratios  = a_ref_ratios;
  m_prob_domain = a_prob_domain;  
  
  Vector<int> vecVars;
  a_eqSys.vectorVars(vecVars);
  
  int i;
  
  m_nVectors = vecVars.size();
  m_vectors  = new int[m_nVectors];
  for (i=0; i < m_nVectors; ++i) m_vectors[i] = vecVars[i];    
  
  a_eqSys.velocityVars(vecVars);
  m_nVelVectors = vecVars.size();
  m_velVectors  = new int[m_nVelVectors];
  for (i=0; i < m_nVelVectors; ++i) m_velVectors[i] = vecVars[i];  
   
}

// Input parameters  
void CoordinateSystemHandler::input( ParmParse & parser, int verbosity )
{
  m_verbosity = verbosity;
}

const RealVect &  CoordinateSystemHandler::smallEndCurv () const
{
  return m_loCurv;
}
  
const RealVect &  CoordinateSystemHandler::bigEndCurv () const
{
  return m_hiCurv;
}


CoordinateSystemHandler::eCoordinateSystem CoordinateSystemHandler::coordinateSystem( void ) const
{
  return m_CoordinateSystem;
}

// Modifies conservative variables in a_interp to  values
  // that will be used PiecewiseLinearFillPatchMHDAM for filling ghost cells
void CoordinateSystemHandler::PLFP_ConsToInterpVars(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level)
{
}

void CoordinateSystemHandler::getNodeCoords( NodeFArrayBox & a_coords,                              
                              int         a_level)
{
   // node-centered FAB
  FArrayBox& thisFAB = a_coords.getFab();  
  getNodeCoords(thisFAB, thisFAB.box(), a_level);  
}                              

void CoordinateSystemHandler::getNodeCoordsCartesian( NodeFArrayBox & a_coords,                              
                              int         a_level)
{
  // node-centered FAB
  FArrayBox& thisFAB = a_coords.getFab();  
  getNodeCoordsCartesian(thisFAB, thisFAB.box(), a_level);  
}                              

 
  // Reverse operation to PLFP_ConsToInterpVars
void CoordinateSystemHandler::PLFP_InterpVarsToCons(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level)
{
}
                           
  // Transform Cartesian components of a vector into curvilinear
void CoordinateSystemHandler::transCartesianVectToCurv(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level)
{
}

void CoordinateSystemHandler::transCartesianVectToCurv(
                        FArrayBox & a_U,                   
                  const Box       & a_box,
                  int               a_level)
{
  transCartesianVectToCurv(a_U, m_vectors, m_nVectors, a_box, a_level);
}                  

// Transform curvilinear components of a vector into Cartesian
void CoordinateSystemHandler::transCurvVectToCartesian(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level)
{
}

// Transform curvilinear components of a vector into Cartesian
void CoordinateSystemHandler::transCurvVectToCartesian(
                        FArrayBox & a_U,                   
                  const Box       & a_box,
                  int               a_level)
{
  transCurvVectToCartesian(a_U, m_vectors, m_nVectors, a_box, a_level);
}

  
void CoordinateSystemHandler::multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
  FArrayBox areas(a_box,1);
  getAreas(areas, a_box, a_dir, a_level);
  
  int iBgn  = a_comps.begin();
  int iEnd  = a_comps.end();
  
  FORT_MULTIPLYFLUXESBYAREA(
    CHF_FRA(a_F),
    CHF_CONST_INT(iBgn),
    CHF_CONST_INT(iEnd),
    CHF_CONST_FRA1(areas,0),
    CHF_BOX(a_box));
  
}

void CoordinateSystemHandler::multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const FArrayBox & a_areas,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
  int iBgn  = a_comps.begin();
  int iEnd  = a_comps.end();
  
   FORT_MULTIPLYFLUXESBYAREA(
    CHF_FRA(a_F),
    CHF_CONST_INT(iBgn),
    CHF_CONST_INT(iEnd),
    CHF_CONST_FRA1(a_areas,0),
    CHF_BOX(a_box));
}

void CoordinateSystemHandler::transFluxesForUpdateState(                        
                        FArrayBox & a_F,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level)
{
}

bool CoordinateSystemHandler::scaleFineFluxes()
{
  return true;
}

void CoordinateSystemHandler::modifyFluxesFRincrementCoarse(
                        FArrayBox & fluxC, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level)
{
  fluxC.copy(flux);
}                        
                        
void CoordinateSystemHandler::modifyFluxesFRincrementFine(
                        FArrayBox & fluxF, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level)
{
  fluxF.copy(flux);
}

#ifdef CH_USE_HDF5
// Reads coordinate specific data  
void CoordinateSystemHandler::readGeomInfo(hid_t a_fileID)
{
  bool dirs[CH_SPACEDIM] = {D_DECL(true,true,true)};
  readGeomInfo(a_fileID, dirs);
}

void CoordinateSystemHandler::readGeomInfo(hid_t a_fileID, bool a_dirs[CH_SPACEDIM])
{
}

// Writes coordinate specific data
void CoordinateSystemHandler::writeGeomInfo(hid_t a_fileID)
{
}
#endif