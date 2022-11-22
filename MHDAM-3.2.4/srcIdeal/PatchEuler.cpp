#include <limits>

#include "BoxIterator.H"
#include "IntVectSet.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"

#include "PatchMHDAMF_F.H"
#include "PatchEuler.H"
#include "PatchEulerF_F.H"
#include "SourceCalculator.H"
#include "CSHandler.H"
#include "LevelSetF_F.H"

#include "LGintegrator.H"

PatchEuler::PatchEuler():PatchMHDAM()
{
  
}

// Factory method - this object is its own factory.  It returns a pointer
// to new PatchMHDAM object with its initial and boundary condtions, slope
// parameters, and artificial viscosity information defined.
PatchMHDAM* PatchEuler::new_patchMHDAM() const
{
                                                          // Make the new object
  PatchMHDAM* retval = static_cast<PatchMHDAM*>(new PatchEuler());

  copyTo( retval );
                                                        // Return the new object
  return retval;
}

//                                      Copy internal data to another PatchMHDAM
void PatchEuler::copyTo( PatchMHDAM * pPatch ) const
{
  PatchMHDAM::copyTo( pPatch );
}

//                                       Compute dt and returns a cell with the minimum dt. 
Real PatchEuler::computeDt( const FArrayBox& a_U,
                            const Box&     a_box,
                            IntVect&       a_minDtCell)
{
  CH_assert(isDefined());
  CH_assert(a_U.contains(a_box));
  
  if (m_bLSonly == true)
  {
    return computeDtLevelSet( a_U, a_box, a_minDtCell);    
  }
  
  FArrayBox dtFab(a_box,1);
  dtFab.setVal(0.0);
    
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))  
  {
    Real dx = m_csh->dx(0,m_level);


    int startRho = URHO;

    FORT_MAXWAVESPEED_E( CHF_FRA1(dtFab,0),
                       CHF_CONST_FRA(a_U),
                       CHF_CONST_INT(startRho),
                       CHF_BOX(a_box));        
                         
    dtFab.invert(dx);
  }
  
  if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Spherical) )
  {        
    dtFab.setVal(numeric_limits<Real>::max());
    
    FArrayBox USph(a_box,a_U.nComp());
    USph.copy(a_U);
        
    m_csh->transCartesianVectToCurv(USph,a_box,m_level);    
    
    int startRho = URHO;
        
    FORT_MINDT_SPHERICAL_E(CHF_FRA1(dtFab,0),
          CHF_CONST_FRA(USph),   
          CHF_CONST_INT(startRho),
          CHF_CONST_INT(m_level),
          CHF_BOX(a_box));
                    
  }

  
  Real dt = m_PhPr->computeDt(a_U, dtFab, a_box, a_minDtCell);
  
  return dt;  
}


void PatchEuler::tagCells(const FArrayBox&  a_U,
                          const Box&        a_box,
                                IntVectSet& a_tags)
{
  int dir;

  const Box& b = a_box;
  IntVectSet& localTags = a_tags;
  
  BaseFab<int> lockedCells(a_box,1);
  m_PhPr->lockedCellsRegrid(lockedCells, a_U, a_box);

  // Compute relative gradient of density
  if ((RefineParams[REF_RHO].m_PerformRef) &&
      (RefineParams[REF_RHO].m_maxlevel >= m_level))
  {
    FArrayBox gradDensityFab(b,SpaceDim);

    for (dir = 0; dir < SpaceDim; ++dir)
    {
      const Box bCenter = b & grow(m_domain,-BASISV(dir));
      const Box bLo     = b & adjCellLo(bCenter,dir);
      const int hasLo = ! bLo.isEmpty();
      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi = ! bHi.isEmpty();

      FORT_GETRELGRAD(CHF_FRA1(gradDensityFab,dir),
                    CHF_CONST_FRA1(a_U,URHO),
                    CHF_CONST_INT(dir),
                    CHF_BOX(bLo),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(bHi),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(bCenter));
    }

    FArrayBox gradDensityMagFab(b,1);
    FORT_MAGNITUDE(CHF_FRA1(gradDensityMagFab,0),
                 CHF_CONST_FRA(gradDensityFab),
                 CHF_BOX(b));

    // Create tags based on undivided gradient of density
    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();      

      if ((lockedCells(iv) == 0) && (gradDensityMagFab(iv) >= RefineParams[REF_RHO].m_Threshold))
      {
        localTags |= iv;
      }
      
    }
  }

  m_PhPr->tagCells(a_U, a_box, a_tags);
}

/// This method returns interval that includes indices of all variables used in adaptation criteria.
void PatchEuler::tagCellVarsInterval(Interval& a_interval)
{
  if (RefineParams[REF_RHO].m_PerformRef)
  {
    a_interval.define(URHO,URHO);
  }
}


void PatchEuler::fluxesHancock(       FArrayBox & a_FMinus,
                                      FArrayBox & a_FPlus,
                                const FArrayBox & a_WMinus,
                                const FArrayBox & a_WPlus,
                                const int &       a_dir,
                                const Box &       a_box)
{  
    int startRho  = URHO;
    int startRhoW = WRHO;
    
    FORT_FLUXESHANCOCK_E(
        CHF_FRA(a_FMinus), 
        CHF_CONST_FRA(a_WMinus),          
        CHF_CONST_INT(a_dir),
        CHF_CONST_INT(startRho),
        CHF_CONST_INT(startRhoW),
        CHF_BOX(a_box));
        
    FORT_FLUXESHANCOCK_E(
        CHF_FRA(a_FPlus), 
        CHF_CONST_FRA(a_WPlus),          
        CHF_CONST_INT(a_dir),
        CHF_CONST_INT(startRho),
        CHF_CONST_INT(startRhoW),
        CHF_BOX(a_box));
    
}                                


void PatchEuler::fluxesRP(       FArrayBox & a_F,
                           const FArrayBox & a_WLeft,
                           const FArrayBox & a_WRight,
                           const int &       a_dir,
                           const Box &       a_box)
{
  m_RS->fluxes( a_F, a_WLeft, a_WRight, a_dir, WRHO, a_box );
}

void PatchEuler::addExplicitSources(       FArrayBox    & a_U,
                                     const FArrayBox    & a_W,
                                     const FArrayBox    & a_SOut,
                                           FArrayBox    & a_S,
                                     const FluxBox      & a_Bn,
                                           FArrayBox    & a_divB,
                                           BaseFab<int> & a_REG,
                                     const Real         & a_dt,   
                                     const FArrayBox    & a_scale,                                    
                                     const Box          & a_box)
{
  CH_assert(a_W.box().contains(a_box));

  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  int startRho  = URHO;
  int startRhoW = WRHO;

  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);  
    
    FORT_SOURCEAXISYMMETRIC_E( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(startRho),
                               CHF_CONST_INT(startRhoW),
                               CHF_BOX(a_box) );
                               
    int iBGN = URHO;
    int iEND = UENG;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                       CHF_CONST_FRA(a_S),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(a_box) );
    
  }
  if (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym)
  {                          
    FORT_SOURCEAXISYMMETRIC_POLAR_E( CHF_FRA(a_S),
                           CHF_CONST_FRA(a_W),
                           CHF_CONST_REAL(a_dt), 
                           CHF_CONST_INT(startRho),
                           CHF_CONST_INT(startRhoW),  
                           CHF_CONST_INT(m_level),
                           CHF_BOX(a_box) );   
                           
    int iBGN = URHO;
    int iEND = UENG;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                       CHF_CONST_FRA(a_S),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(a_box) );
  }
    
  
  m_PhPr->explicitSource( a_U, a_S, a_W, a_REG, a_dt, a_box );

  SourceCalculator * pSoCal = m_PhPr->getSourceCalculator();
  if( pSoCal != NULL )
  {    
    pSoCal->addExternalSources( a_U, a_SOut, a_W, a_dt, m_csh->dx(0,m_level), a_box );
  }
}

void PatchEuler::addExplicitSources(       FArrayBox    & a_U,
                                     const FArrayBox    & a_W,
                                     const FArrayBox    & a_SOut,
                                           FArrayBox    & a_S,
                                     const FluxBox      & a_Bn,
                                           FArrayBox    & a_divB,
                                     const FluxBox      & a_Un,
                                           BaseFab<int> & a_REG,
                                     const Real         & a_dt,   
                                     const FArrayBox    & a_scale,                                    
                                     const Box          & a_box)
{
  CH_assert(a_W.box().contains(a_box));

  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  int startRho  = URHO;
  int startRhoW = WRHO;

  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);  
    
    FORT_SOURCEAXISYMMETRIC_E( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(startRho),
                               CHF_CONST_INT(startRhoW),
                               CHF_BOX(a_box) );
                               
    int iBGN = URHO;
    int iEND = UENG;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                       CHF_CONST_FRA(a_S),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(a_box) );
    
  }
  if (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym)
  {                          
    FORT_SOURCEAXISYMMETRIC_POLAR_E( CHF_FRA(a_S),
                           CHF_CONST_FRA(a_W),
                           CHF_CONST_REAL(a_dt), 
                           CHF_CONST_INT(startRho),
                           CHF_CONST_INT(startRhoW),  
                           CHF_CONST_INT(m_level),
                           CHF_BOX(a_box) );   
                           
    int iBGN = URHO;
    int iEND = UENG;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                       CHF_CONST_FRA(a_S),
                       CHF_CONST_INT(iBGN),
                       CHF_CONST_INT(iEND),
                       CHF_BOX(a_box) );
  }
    
  
  m_PhPr->explicitSource( a_U, a_S, a_W, a_REG, a_dt, a_box );

  SourceCalculator * pSoCal = m_PhPr->getSourceCalculator();
  if( pSoCal != NULL )
  {    
    pSoCal->addExternalSources( a_U, a_SOut, a_W, a_dt, m_csh->dx(0,m_level), a_box );
  }
}

void PatchEuler::postprocessing(       FArrayBox & a_U,
                                 const FArrayBox & a_Uold,
                                 const Real      & a_dt,
                                 const Box       & a_box)
{
  int startRho = URHO;

  FORT_POSTPROCESSING_E( CHF_FRA(a_U),
                         CHF_CONST_FRA(a_Uold),
                         CHF_CONST_REAL(a_dt),
                         CHF_CONST_INT(startRho),
                         CHF_BOX(a_box) );
                         
  if (m_eqSys->numTrackingSurfaces()>0)
  {
    int ils = m_eqSys->lsIndexCons(0);
    int nls = m_eqSys->numTrackingSurfaces();    
    
    FORT_POSTPROCESSING_LS(CHF_FRA(a_U),       
                CHF_CONST_FRA(a_Uold),                                         
                CHF_CONST_INT(ils),
                CHF_CONST_INT(nls),
                CHF_BOX(a_box));        
  }
}

void PatchEuler::preprocessing ( const FArrayBox    & a_W,
                                       FArrayBox    & a_S,
                                       BaseFab<int> & a_R,
                                 const Box          & a_box)
{
}

