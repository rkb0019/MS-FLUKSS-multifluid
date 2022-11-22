#include <limits>

#include <BoxIterator.H>
#include <IntVectSet.H>

#include "PatchMHDMF.H"
#include "PatchMHDMFF_F.H"
#include "PatchIdealMHDF_F.H"
#include "PatchEulerF_F.H"
#include "PatchMHDAMF_F.H"
#include "SourceCalculator.H"
#include "MultiFluidProblem.H"
#include "DivCleaningF_F.H"
#include "CSHandler.H"
#include "EquationSystem.H"
#include "EqSysMHDMF.H"
#include "LevelSetF_F.H"
#include "PickupIonsF_F.H"
#include "LoHiCenter.H"
#include "LoHiSide.H"

#include "LGintegrator.H"

// Temporarily for source terms
#include "DednerF_F.H"


PatchMHDMF::PatchMHDMF(int a_nFluids):PatchMHDAM()
{
  CH_assert(a_nFluids>=1);
  m_nFluids = a_nFluids;
  m_iDivBMethod        = 1;
}

// Input parameters
void PatchMHDMF::input( ParmParse & parser, int verbosity )
{  
}

// Factory method - this object is its own factory.  It returns a pointer
// to new PatchMHDAM object with its initial and boundary condtions, Riemann solvers,
// reconstruction parametrs, coordinate system and time approximation information defined.
PatchMHDAM* PatchMHDMF::new_patchMHDAM() const
{
                                                          // Make the new object
  PatchMHDAM* retval = static_cast<PatchMHDAM*>(new PatchMHDMF(m_nFluids  ));

  copyTo( retval );
                                                        // Return the new object
  return retval;
}

//                                      Copy internal data to another PatchMHDAM
void PatchMHDMF::copyTo( PatchMHDAM * pPatch ) const
{
  PatchMHDAM::copyTo( pPatch );

  PatchMHDMF* pP = dynamic_cast<PatchMHDMF*>(pPatch);
  if (pP == NULL) MayDay::Error("PatchMHDMF::copyTo. Wrong argument");

  pP->m_nFluids = m_nFluids;  
}


//                                       Compute dt and returns a cell with the minimum dt.
Real PatchMHDMF::computeDt( const FArrayBox& a_U,
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

  int iFluid, startRho;

  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);


  if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) ||
       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
  {
    dtFab.setVal(0.0);

    Real dx = m_csh->dx(0,m_level);

    FORT_MAXWAVESPEED(CHF_FRA1(dtFab,0),
                      CHF_CONST_FRA(a_U),
                      CHF_BOX(a_box));

    for (iFluid = 1; iFluid < m_nFluids; ++iFluid)
    {
      startRho = eqSys->densityIndexCons(iFluid);

      FORT_MAXWAVESPEED_E(CHF_FRA1(dtFab,0),
                           CHF_CONST_FRA(a_U),
                           CHF_CONST_INT(startRho),
                           CHF_BOX(a_box));
    }


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

    FORT_MINDT_SPHERICAL(CHF_FRA1(dtFab,0),
                CHF_CONST_FRA(USph),
                CHF_CONST_INT(m_level),
                CHF_BOX(a_box));

    for (iFluid = 1; iFluid < m_nFluids; ++iFluid)
    {
      startRho = eqSys->densityIndexCons(iFluid);

      FORT_MINDT_SPHERICAL_E(CHF_FRA1(dtFab,0),
                CHF_CONST_FRA(USph),
                CHF_CONST_INT(startRho),
                CHF_CONST_INT(m_level),
                CHF_BOX(a_box));
    }

    m_eqSys->computeDt( USph, dtFab, a_box, m_level, a_minDtCell );
  }

  Real dt = m_PhPr->computeDt(a_U, dtFab, a_box, a_minDtCell);

  return dt;
}


void PatchMHDMF::tagCells(const FArrayBox&  a_U,
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
      const int hasLo   = ! bLo.isEmpty();
      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi   = ! bHi.isEmpty();

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

  if ((RefineParams[REF_B].m_PerformRef) &&
      (RefineParams[REF_B].m_maxlevel >= m_level))
  {
    FArrayBox B(a_U.box(),1);
    B.setVal(0.0);

    int iBGN = UBX;
    int iEND = UBZ;

    FORT_GETVECTMAGNITUDE(
              CHF_FRA1(B,0),
              CHF_CONST_FRA(a_U),
              CHF_CONST_INT(iBGN),
              CHF_CONST_INT(iEND),
              CHF_BOX(a_U.box()));

    FArrayBox gradBFab(b,SpaceDim);

    for (dir = 0; dir < SpaceDim; ++dir)
    {
      const Box bCenter = b & grow(m_domain,-BASISV(dir));
      const Box bLo     = b & adjCellLo(bCenter,dir);
      const int hasLo   = ! bLo.isEmpty();
      const Box bHi     = b & adjCellHi(bCenter,dir);
      const int hasHi = ! bHi.isEmpty();

      FORT_GETRELGRAD(CHF_FRA1(gradBFab,dir),
                    CHF_CONST_FRA1(B,0),
                    CHF_CONST_INT(dir),
                    CHF_BOX(bLo),
                    CHF_CONST_INT(hasLo),
                    CHF_BOX(bHi),
                    CHF_CONST_INT(hasHi),
                    CHF_BOX(bCenter));
    }

    FArrayBox gradBMagFab(b,1);
    FORT_MAGNITUDE(CHF_FRA1(gradBMagFab,0),
                 CHF_CONST_FRA(gradBFab),
                 CHF_BOX(b));

    // Tag where gradient exceeds threshold
    BoxIterator bit(b);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();

      if ((lockedCells(iv) == 0) && (gradBMagFab(iv) >= RefineParams[REF_B].m_Threshold))
      {
        localTags |= iv;
      }
    }

  }
  m_PhPr->tagCells(a_U, a_box, a_tags);

  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {
    IntVectSet newTags; IntVect iv_new; int j;
    IVSIterator ivsit (a_tags);

    int jmin = m_domain.domainBox().smallEnd()[1];
    int jmax = m_domain.domainBox().bigEnd()[2];

    for (ivsit.begin(); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      if (( iv[2] == 0) || (iv[2] == 1))
      {
        IntVect iv_new(D_DECL(iv[0],0,0));
        for (j=jmin; j<=jmax; ++j)
        {
          iv_new[1] = j;
          newTags |= iv_new;
        }
      }
    }
    a_tags |= newTags;
  }


}

/// This method returns interval that includes indices of all variables used in adaptation criteria.
void PatchMHDMF::tagCellVarsInterval(Interval& a_interval)
{
  if ((RefineParams[REF_RHO].m_PerformRef) && (RefineParams[REF_B].m_PerformRef))
  {
    a_interval.define(URHO,UBZ);
  }

  if ((RefineParams[REF_RHO].m_PerformRef) && (!RefineParams[REF_B].m_PerformRef))
  {
    a_interval.define(URHO,URHO);
  }

  if ((!RefineParams[REF_RHO].m_PerformRef) && (RefineParams[REF_B].m_PerformRef))
  {
    a_interval.define(UBX,UBZ);
  }
}

void PatchMHDMF::fluxesHancock(       FArrayBox & a_FMinus,
                                      FArrayBox & a_FPlus,
                                const FArrayBox & a_WMinus,
                                const FArrayBox & a_WPlus,
                                const int &       a_dir,
                                const Box &       a_box )
{
  FORT_FLUXESHANCOCK(
        CHF_FRA(a_FMinus),
        CHF_CONST_FRA(a_WMinus),
        CHF_CONST_INT(a_dir),
        CHF_BOX(a_box));

  FORT_FLUXESHANCOCK(
        CHF_FRA(a_FPlus),
        CHF_CONST_FRA(a_WPlus),
        CHF_CONST_INT(a_dir),
        CHF_BOX(a_box));

  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);

  TurbulenceModel* pTM = eqSys->getTurbulenceModel();

  if( pTM != NULL )
  {
    pTM->primToFlux( a_FMinus, a_WMinus, a_dir, a_box );
    pTM->primToFlux( a_FPlus,  a_WPlus,  a_dir, a_box );
  }

  PickupIons* pPI = eqSys->getPickupIons();

  if( pPI != NULL )
  {

    BaseFab<int> region;
    FArrayBox S;    
    preprocessing( a_WMinus, S, region, a_box );
    pPI->primToFlux( a_FMinus, a_WMinus, region, a_dir, a_box );
//old below, does not use region
//    pPI->primToFlux( a_FMinus, a_WMinus, a_dir, a_box );
    preprocessing( a_WPlus, S, region, a_box );
    pPI->primToFlux( a_FPlus,  a_WPlus, region, a_dir, a_box );
//old below, does not use region
//    pPI->primToFlux( a_FPlus,  a_WPlus,  a_dir, a_box );
  }

  for (int iFluid = 1; iFluid < m_nFluids; ++iFluid)
  {
    int startRho  = eqSys->densityIndexCons(iFluid);
    int startRhoW = eqSys->densityIndexPrim(iFluid);

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

//FF: this was missing, creating artifact in the LS solution with the Hancock time approx
//FF: FLUXESHANCOCK_LS added to LevelSetF.ChF
  if (m_eqSys->numTrackingSurfaces()>0)
  {
    for (int i = 0; i < m_eqSys->numTrackingSurfaces(); ++i)
    {
      int ils = m_eqSys->lsIndexCons(i);

        FORT_FLUXESHANCOCK_LS(
        CHF_FRA(a_FMinus),
        CHF_CONST_FRA(a_WMinus),
        CHF_CONST_INT(a_dir),
        CHF_CONST_INT(ils),
        CHF_BOX(a_box));

        FORT_FLUXESHANCOCK_LS(
        CHF_FRA(a_FPlus),
        CHF_CONST_FRA(a_WPlus),
        CHF_CONST_INT(a_dir),
        CHF_CONST_INT(ils),
        CHF_BOX(a_box));


    } 


  }



}


void PatchMHDMF::fluxesRP(       FArrayBox & a_F,
                           const FArrayBox & a_WLeft,
                           const FArrayBox & a_WRight,
                           const int &       a_dir,
                           const Box &       a_box )
{
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);

  m_RS->fluxes( a_F, a_WLeft, a_WRight, a_dir,  WRHO, a_box );

  TurbulenceModel* pTM = eqSys->getTurbulenceModel();

  if( pTM != NULL )
  {
    pTM->upwindFluxes( a_F, a_WLeft, a_WRight, a_dir, WRHO, a_box );
  }

  PickupIons* pPI = eqSys->getPickupIons();

  if( pPI != NULL )
  {
    BaseFab<int> Region;
    FArrayBox S;    
    preprocessing( a_WLeft, S, Region, a_box );
    pPI->upwindFluxes( a_F, a_WLeft, a_WRight, Region, a_dir, WRHO, a_box );
//old below, uses simple upwinding and does not use region
//this is bad for PUIs at the TS and HP, ok away from shocks
//    pPI->upwindFluxes( a_F, a_WLeft, a_WRight, a_dir, WRHO, a_box );
  }

  for (int iFluid = 1; iFluid < m_nFluids; ++iFluid)
    m_RSGD->fluxes( a_F, a_WLeft, a_WRight, a_dir, eqSys->densityIndexPrim(iFluid), a_box );
}

void PatchMHDMF::computeDivB(   FArrayBox & a_divB,
                          const FluxBox   & a_Bn,
                          const FArrayBox & a_W,
                          const Box       & a_box)
{
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();

  a_divB.setVal(0.0);
  for( int idir = 0; idir < SpaceDim; idir++ )
  {
//    const FArrayBox & Bn = a_Bn[idir];
             // Compute each component of the divergence of the Magnetic Field
    if (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical)
    {
      // This code should be revised in the future to add support of non-unifrom meshes
      // For now it will be compiled (but no gurantee to work!!!)
      Real dx = m_csh->dx(idir,m_level);
      FORT_COMPUTEDIVB_CYL( CHF_FRA1(a_divB,0),
                            CHF_CONST_FRA1(a_Bn[idir],0),
                            CHF_CONST_REAL(dx),
                            CHF_CONST_INT(idir),
                            CHF_BOX(a_box));
    }
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
    {
      FORT_COMPUTEDIVB( CHF_FRA1(a_divB,0),
                        CHF_CONST_FRA1(a_Bn[idir],0),
                        CHF_CONST_INT(idir),
                        CHF_BOX(a_box));
    }

    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) )
    {

      FORT_COMPUTEDIVB_POLAR( CHF_FRA1(a_divB,0),
                  CHF_CONST_FRA1(a_Bn[idir],0),
                  CHF_CONST_INT(m_level),
                  CHF_CONST_INT(idir),
                  CHF_BOX(a_box));
    }
    if ( CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {
      FArrayBox areas(a_Bn[idir].box(),1);
      m_csh->getAreas(areas, areas.box(), idir, m_level);
      FORT_COMPUTEDIVB_SPHERICAL( CHF_FRA1(a_divB,0),
                  CHF_CONST_FRA1(a_Bn[idir],0),
                  CHF_CONST_FRA1(areas,0),
                  CHF_CONST_INT(idir),
                  CHF_BOX(a_box));
    }
  }

  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);
    FORT_AXISYMMETRICDBZBYDZ( CHF_FRA1(a_divB,0),
                              CHF_CONST_FRA(a_W),
                              CHF_CONST_REAL(dx),
                              CHF_BOX(a_box) );
  }

}



void PatchMHDMF::computeBn(       FArrayBox & a_Bn,
                            const FArrayBox & a_W,
                            const FArrayBox & a_WMinus,
                            const FArrayBox & a_WPlus ,
                            const int       & a_method,
                            const int       & a_dir,
                            const Box       & a_box )
{
  if( a_method == 0 ) return;

  if( m_eqSys->numStates() <= WBX ) return;

  CH_assert( a_Bn.box().contains(a_box) );

  if( m_verbosity >= 3 )
  {
    pout() << "computeBn:in" << endl;
  }

  if( a_method == 1 )
  {
         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
    FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);

                                        // Computes a_Bn for all edges
                                        // that are not on the physical boundary
    CH_assert( shiftWLeft.box().contains(a_box) );
    CH_assert( shiftWRight.box().contains(a_box) );

    FORT_COMPUTEBN( CHF_FRA1(a_Bn,0),
                    CHF_CONST_FRA(shiftWLeft),
                    CHF_CONST_FRA(shiftWRight),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(a_box));

                                         // Now fill Bn on the physical boundary
    Box BnBox = a_Bn.box();
                                   // In periodic case, this doesn't do anything
    if( !m_domain.isPeriodic(a_dir) )
    {
      SideIterator si;
      for( si.begin(); si.ok(); ++si )
      {
        Side::LoHiSide side = si();

        Box tmp = BnBox;
                            // Determine which side and thus shifting directions
        int sign;
        if( side == Side::Lo )
        {
          sign = -1;
        }
        else
        {
          sign = 1;
        }

        tmp.shiftHalf(a_dir,sign);

                                 // Is there a domain boundary next to this grid
        if( !m_domain.contains(tmp) )
        {
          tmp &= m_domain;

                          // Find the strip of cells next to the domain boundary
          Box boundaryBox;
          if( side == Side::Lo )
          {
            boundaryBox = bdryLo(tmp,a_dir);
          }
          else
          {
            boundaryBox = bdryHi(tmp,a_dir);
          }

          CH_assert(      boundaryBox.type(a_dir) == BnBox.type(a_dir));
          CH_assert( shiftWLeft.box().type(a_dir) == BnBox.type(a_dir));
          CH_assert(shiftWRight.box().type(a_dir) == BnBox.type(a_dir));

          if( side == Side::Lo )
          {
            a_Bn.copy( shiftWRight, boundaryBox, WBX+a_dir, boundaryBox, 0, 1 );
          }
          if( side == Side::Hi )
          {
            a_Bn.copy( shiftWLeft,  boundaryBox, WBX+a_dir, boundaryBox, 0, 1 );
          }
        }
      }
    }
                            // Shift the left and right primitive variable boxes
                            // back to their original position
    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
  } else {
    FORT_COMPUTEBNCD( CHF_FRA1(a_Bn,0),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
  }

  if( m_verbosity >= 3 )
  {
    pout() << "computeBn:out" << endl;
  }
}

void PatchMHDMF::computeDivU(   FArrayBox & a_divU,
                          const FluxBox   & a_Un,
                          const FArrayBox & a_W,
                          const Box       & a_box)
{
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();

  a_divU.setVal(0.0);
  for( int idir = 0; idir < SpaceDim; idir++ )
  {
//    const FArrayBox & Un = a_Un[idir];
             // Compute each component of the divergence of the velocity
//    if (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical)
//    {
//      // This code should be revised in the future to add support of non-unifrom meshes
//      // For now it will be compiled (but no gurantee to work!!!)
//      Real dx = m_csh->dx(idir,m_level);
//      FORT_COMPUTEDIVU_CYL( CHF_FRA1(a_divU,0),
//                            CHF_CONST_FRA1(a_Un[idir],0),
//                            CHF_CONST_REAL(dx),
//                            CHF_CONST_INT(idir),
//                            CHF_BOX(a_box));
//    }
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
    {
      FORT_COMPUTEDIVU( CHF_FRA1(a_divU,0),
                        CHF_CONST_FRA1(a_Un[idir],0),
                        CHF_CONST_INT(idir),
                        CHF_BOX(a_box));
    }

//    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) ||
//         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) )
//    {
//
//      FORT_COMPUTEDIVU_POLAR( CHF_FRA1(a_divU,0),
//                  CHF_CONST_FRA1(a_Un[idir],0),
//                  CHF_CONST_INT(m_level),
//                  CHF_CONST_INT(idir),
//                  CHF_BOX(a_box));
//    }
    if ( CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {
      FArrayBox areas(a_Un[idir].box(),1);
      m_csh->getAreas(areas, areas.box(), idir, m_level);
      FORT_COMPUTEDIVU_SPHERICAL( CHF_FRA1(a_divU,0),
                  CHF_CONST_FRA1(a_Un[idir],0),
                  CHF_CONST_FRA1(areas,0),
                  CHF_CONST_INT(idir),
                  CHF_BOX(a_box));
    }
  }

/*  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);
    FORT_AXISYMMETRICDBZBYDZ( CHF_FRA1(a_divU,0),
                              CHF_CONST_FRA(a_W),
                              CHF_CONST_REAL(dx),
                              CHF_BOX(a_box) );
  }*/

    Box divUbox;
    divUbox = a_divU.box();
    FArrayBox scale(divUbox, 1);
    m_csh->scalingFactor(scale, divUbox, m_level);
    a_divU.mult(scale);

}

//needed to compute divU
void PatchMHDMF::computeUn(       FArrayBox & a_Un,
                            const FArrayBox & a_W,
                            const FArrayBox & a_WMinus,
                            const FArrayBox & a_WPlus ,
                            const int       & a_method,
                            const int       & a_dir,
                            const Box       & a_box )
{
  CH_assert( a_Un.box().contains(a_box) );

  if( m_verbosity >= 3 )
  {
    pout() << "computeUn:in" << endl;
  }

  if( a_method == 1 )
  {
         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
    FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);

                                        // Computes a_Un for all edges
                                        // that are not on the physical boundary
    CH_assert( shiftWLeft.box().contains(a_box) );
    CH_assert( shiftWRight.box().contains(a_box) );

    FORT_COMPUTEUN( CHF_FRA1(a_Un,0),
                    CHF_CONST_FRA(shiftWLeft),
                    CHF_CONST_FRA(shiftWRight),
                    CHF_CONST_INT(a_dir),
                    CHF_BOX(a_box));

                                         // Now fill Un on the physical boundary
    Box UnBox = a_Un.box();
                                   // In periodic case, this doesn't do anything
    if( !m_domain.isPeriodic(a_dir) )
    {
      SideIterator si;
      for( si.begin(); si.ok(); ++si )
      {
        Side::LoHiSide side = si();

        Box tmp = UnBox;
                            // Determine which side and thus shifting directions
        int sign;
        if( side == Side::Lo )
        {
          sign = -1;
        }
        else
        {
          sign = 1;
        }

        tmp.shiftHalf(a_dir,sign);

                                 // Is there a domain boundary next to this grid
        if( !m_domain.contains(tmp) )
        {
          tmp &= m_domain;

                          // Find the strip of cells next to the domain boundary
          Box boundaryBox;
          if( side == Side::Lo )
          {
            boundaryBox = bdryLo(tmp,a_dir);
          }
          else
          {
            boundaryBox = bdryHi(tmp,a_dir);
          }

          CH_assert(      boundaryBox.type(a_dir) == UnBox.type(a_dir));
          CH_assert( shiftWLeft.box().type(a_dir) == UnBox.type(a_dir));
          CH_assert(shiftWRight.box().type(a_dir) == UnBox.type(a_dir));

          if( side == Side::Lo )
          {
            a_Un.copy( shiftWRight, boundaryBox, WVELX+a_dir, boundaryBox, 0, 1 );
          }
          if( side == Side::Hi )
          {
            a_Un.copy( shiftWLeft,  boundaryBox, WVELX+a_dir, boundaryBox, 0, 1 );
          }
        }
      }
    }
                            // Shift the left and right primitive variable boxes
                            // back to their original position
    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
  } else {
    FORT_COMPUTEUNCD( CHF_FRA1(a_Un,0),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(a_dir),
                      CHF_BOX(a_box));
  }

  if( m_verbosity >= 3 )
  {
    pout() << "computeUn:out" << endl;
  }
}
void PatchMHDMF::computeDivB(        FArrayBox& a_divB,                              
                              const  FArrayBox& a_U,
                              const  Box&       a_box)
{
  // Get the number of various variables
  int numFlux  = m_eqSys->numFluxes();
  int numPrim  = m_eqSys->numPrimitives();
  int numSlope = m_eqSys->numSlopes();
  int numCons  = a_U.nComp();
  
  Box UBox     = a_U.box();
  Box WBox     = UBox;
  Box WBoxB    = UBox;
      WBoxB   &= m_domain;
    
                                     // Primitive variables    
  FArrayBox W( WBox, numPrim );                                     
      
  Box valueBox = a_box;  
  valueBox &= m_domain;
  
  Box slopeBox = valueBox;
  slopeBox.grow( 1 );
  
  // Calculate the primitive variables from the conserved variables
  m_eqSys->stateToPrim( W, a_U, WBoxB );
  m_csh->transCartesianVectToCurv(W, WBoxB, m_level);

  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    m_PhPr->fillGhostCells( W, a_U, dir, m_currentTime );
  }
    
    
  FArrayBox WMinus[SpaceDim];  
  FArrayBox WPlus [SpaceDim];
  FluxBox   Bn(valueBox,1);  
    

////////////////////////////////////////////////////////////////  Reconstruction

  for( int dir = 0; dir < SpaceDim; dir++ )
  {
    WMinus[dir].clear();
    WPlus [dir].clear();

                    // Size the intermediate, extrapolated primitive variables
    WMinus[dir].resize( slopeBox, numPrim  );
    WPlus [dir].resize( slopeBox, numPrim  );

                                       // Compute primitive variables on faces
    m_Reconstruction->faceValues( W, WMinus[dir], WPlus[dir], m_level, dir, slopeBox, m_eqSys, m_csh );
         
    m_Reconstruction->checkPositivity( WMinus[dir], WPlus[dir], W, slopeBox );
  }


  // Boxes for face centered state - used for MUSCL approach
  Box faceBox[SpaceDim];  
  
  for( int dir = 0; dir < SpaceDim; ++dir )
  {
     //  This box is face centered in direction "dir", is one bigger than the
     // input box in all directions except "dir" and stays one cell away from
     // the domain boundary in "dir"
    faceBox[dir] = valueBox;
    faceBox[dir].grow(dir,1);//faceBox[dir].grow(1);
    faceBox[dir] &= m_domain;
    faceBox[dir].grow(dir,-1);
    faceBox[dir].surroundingNodes(dir);   
  }
  
  for( int dir = 0; dir < SpaceDim; dir++ )
  {       
    computeBn( Bn[dir], W, WMinus[dir], WPlus[dir], m_iDivBMethod, dir, faceBox[dir]);
  }
  
  computeDivB(a_divB, Bn, W, a_box);
            
}                              

void PatchMHDMF::correctBn(       FArrayBox & a_Bn,
                            const FArrayBox & a_phi,
                            const int       & a_dir,
                            const Box       & a_box )
{

  if( m_eqSys->numStates() <= WBX ) return;

  Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);

  CH_assert( a_Bn.box().contains(a_box) );

  if( m_verbosity >= 3 )
  {
    pout() << "correctBn:in" << endl;
  }

  FORT_CORRECTBN( CHF_FRA1(a_Bn,0),
                  CHF_CONST_FRA1(a_phi,0),
                  CHF_CONST_REAL(dx),
                  CHF_CONST_INT(a_dir),
                  CHF_BOX(a_box));

  if( m_verbosity >= 3 )
  {
    pout() << "computeBn:out" << endl;
  }
}

//old version, unchanged
void PatchMHDMF::addExplicitSources(       FArrayBox    & a_U,
                                     const FArrayBox    & a_W,
                                     const FArrayBox    & a_SOut,
                                           FArrayBox    & a_S,
                                     const FluxBox      & a_Bn,
                                           FArrayBox    & a_divB,
                                           BaseFab<int> & a_REG,
                                     const Real         & a_dt,
                                     const FArrayBox    & a_scale,
                                     const Box          & a_box )
{
  CH_assert(a_W.box().contains(a_box));
  

  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);

  if( m_verbosity >= 3 )
  {
    pout() << "addExplicitSources:in" << endl;
    pout() << "UBox   :"; a_U.box().p();
    pout() << "WBox   :"; a_W.box().p();
    pout() << "SOutBox:"; a_SOut.box().p();
    pout() << "SBox   :"; a_S.box().p();
    pout() << "divBBox:"; a_divB.box().p();
  }

  int iBGN = UMOMX;
  int iEND = UBZ;
  
  FArrayBox divB( a_box, 1 );
  computeDivB( divB, a_Bn, a_W, a_box);

  if( m_b8waveUse == true ) // Powell 8-wave divergence cleaning
  {
                                         // The divergence of the Magnetic Field
    a_S.setVal(0.0);
    
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
    {
      Real dx = m_csh->dx(0,m_level);
      double dtbydx = a_dt/dx;

      FORT_SOURCE8WAVES( CHF_FRA(a_S),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_FRA1(divB,0),
                         CHF_CONST_REAL(dtbydx),
                         CHF_BOX(a_box) );
    }
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) )
    {
      FORT_SOURCE8WAVES_POLAR( CHF_FRA(a_S),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_FRA1(divB,0),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_INT(m_level),
                 CHF_BOX(a_box) );
    }
    if ( CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {
      //divB*=a_scale;
      FORT_SOURCE8WAVES_SPHERICAL( CHF_FRA(a_S),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_FRA1(divB,0),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_FRA1(a_scale,0),
                 CHF_CONST_INT(m_level),
                 CHF_BOX(a_box));
    }

    iBGN   = UMOMX;
    iEND   = UBZ;

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

    if (a_divB.nComp() > 0) a_divB.copy(divB);
  }
  
     
                                                     // Geometrical source terms
  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);
    iBGN   = URHO;
    if (m_nFluids == 1)
    {
      iEND = UBZ;
      FORT_SOURCEAXISYMMETRIC( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_BOX(a_box) );
                               
      if (eqSys->isUsedCorrectionPotential() == 1) iEND = UBZ+1;
      
    } else
    {
      FORT_SOURCEAXISYMMETRIC_MF( CHF_FRA(a_S),
                                  CHF_CONST_FRA(a_W),
                                  CHF_CONST_INT(iRhoN),
                                  CHF_CONST_INT(m_nFluids),
                                  CHF_CONST_REAL(a_dt),
                                  CHF_CONST_REAL(dx),
                                  CHF_BOX(a_box) );
      iEND = eqSys->densityIndexCons(m_nFluids-1)+UENG;
    }
        
    if (eqSys->isUsedCorrectionPotential() == 1)
    {
      FORT_SOURCEAXISYMMETRICDEDNER( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_BOX(a_box) );            
    }

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

  }

  if (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical)
  {
    Real dx = m_csh->dx(0,m_level);
    if (m_nFluids==1)
    {
      FORT_SOURCEAXISYMMETRIC_CYL( CHF_FRA(a_S),
                                   CHF_CONST_FRA(a_W),
                                   CHF_CONST_REAL(a_dt),
                                   CHF_CONST_REAL(dx),
                                   CHF_BOX(a_box) );
    } else
    {
      FORT_SOURCEAXISYMMETRIC_CYL_MF( CHF_FRA(a_S),
                                      CHF_CONST_FRA(a_W),
                                      CHF_CONST_INT(iRhoN),
                                      CHF_CONST_INT(m_nFluids),
                                      CHF_CONST_REAL(a_dt),
                                      CHF_CONST_REAL(dx),
                                      CHF_BOX(a_box) );
    }

    iBGN   = UMOMY;
    iEND   = (m_nFluids < 2 ? UMOMY :  UMOMY + eqSys->densityIndexCons(m_nFluids-1)); // UMOMY1+(m_nFluids-2)*(URHO2-URHO1)

  }
  if (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym)
  {
    FORT_SOURCEAXISYMMETRIC_POLAR( CHF_FRA(a_S),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_INT(m_level),
                             CHF_BOX(a_box) );

    iBGN = URHO;
    iEND = UBX;

    for (int iFluid = 1; iFluid < m_nFluids; ++iFluid)
    {
      int startRho = eqSys->densityIndexCons(iFluid);

      FORT_SOURCEAXISYMMETRIC_POLAR_E( CHF_FRA(a_S),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_INT(startRho),
                             CHF_CONST_INT(startRho),
                             CHF_CONST_INT(m_level),
                             CHF_BOX(a_box) );
    }
    if (m_nFluids > 1) iEND = eqSys->densityIndexCons(m_nFluids-1)+UENG;

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

  }
  
  if (eqSys->isUsedCorrectionPotential() == 1)
  {
    Real dx = m_csh->dx(0,m_level);
    int iRho = eqSys->densityIndexCons(0);
    /*FORT_DEDNERSOURCETERMS(
         CHF_FRA(a_S),
         CHF_CONST_FRA(a_W),
         CHF_CONST_FRA1(divB,0),
         CHF_CONST_INT(iRho),
         CHF_CONST_REAL(dx),     
         CHF_CONST_REAL(a_dt),     
         CHF_BOX(a_box));
    int iBgn = UENG;
    int iEnd = UBZ;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                 CHF_CONST_FRA(a_S),
                 CHF_CONST_INT(iBgn),
                 CHF_CONST_INT(iEnd),
                 CHF_BOX(a_box) );*/
                 
    FORT_PSIPARABOLIC(
         CHF_FRA(a_U),         
         CHF_CONST_INT(iRho),
         CHF_CONST_REAL(dx),     
         CHF_CONST_REAL(a_dt),     
         CHF_BOX(a_box));
  }

  eqSys->explicitSource( a_U, a_S, a_W, a_dt, m_level, a_box );
  m_PhPr->explicitSource( a_U, a_S, a_W, a_REG, a_dt, a_box );

  SourceCalculator * pSoCal = m_PhPr->getSourceCalculator();
  if( pSoCal != NULL )
  {
    Real dx = (m_csh->constStep(0)==true ? m_csh->dx(0,m_level) : 1.0);
    pSoCal->addExternalSources( a_U, a_SOut, a_W, a_dt, dx, a_box );
  }

  if( m_verbosity >= 3 )
  {
    pout() << "addExplicitSources:out" << endl;
  }
}

//new version, passes region and calculates pdivu rather than ugradp for PUIs
void PatchMHDMF::addExplicitSources(       FArrayBox    & a_U,
                                     const FArrayBox    & a_W,
                                     const FArrayBox    & a_SOut,
                                           FArrayBox    & a_S,
                                     const FluxBox      & a_Bn,
                                           FArrayBox    & a_divB,
                                     const FluxBox      & a_Un,
                                           BaseFab<int> & a_REG,
                                     const Real         & a_dt,
                                     const FArrayBox    & a_scale,
                                     const Box          & a_box )
{
  CH_assert(a_W.box().contains(a_box));
  

  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);

  if( m_verbosity >= 3 )
  {
    pout() << "NEWaddExplicitSources:in" << endl;
    pout() << "UBox   :"; a_U.box().p();
    pout() << "WBox   :"; a_W.box().p();
    pout() << "SOutBox:"; a_SOut.box().p();
    pout() << "SBox   :"; a_S.box().p();
    pout() << "divBBox:"; a_divB.box().p();
  }

  int iBGN = UMOMX;
  int iEND = UBZ;
  
  FArrayBox divB( a_box, 1 );
  computeDivB( divB, a_Bn, a_W, a_box);

  if( m_b8waveUse == true ) // Powell 8-wave divergence cleaning
  {
                                         // The divergence of the Magnetic Field
    a_S.setVal(0.0);
    
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
    {
      Real dx = m_csh->dx(0,m_level);
      double dtbydx = a_dt/dx;

      FORT_SOURCE8WAVES( CHF_FRA(a_S),
                         CHF_CONST_FRA(a_W),
                         CHF_CONST_FRA1(divB,0),
                         CHF_CONST_REAL(dtbydx),
                         CHF_BOX(a_box) );
    }
    if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) )
    {
      FORT_SOURCE8WAVES_POLAR( CHF_FRA(a_S),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_FRA1(divB,0),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_INT(m_level),
                 CHF_BOX(a_box) );
    }
    if ( CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {
      //divB*=a_scale;
      FORT_SOURCE8WAVES_SPHERICAL( CHF_FRA(a_S),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_FRA1(divB,0),
                 CHF_CONST_REAL(a_dt),
                 CHF_CONST_FRA1(a_scale,0),
                 CHF_CONST_INT(m_level),
                 CHF_BOX(a_box));
    }

    iBGN   = UMOMX;
    iEND   = UBZ;

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

    if (a_divB.nComp() > 0) a_divB.copy(divB);
  }
  
     
                                                     // Geometrical source terms
  if (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric)
  {
    Real dx = m_csh->dx(0,m_level);
    iBGN   = URHO;
    if (m_nFluids == 1)
    {
      iEND = UBZ;
      FORT_SOURCEAXISYMMETRIC( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_BOX(a_box) );
                               
      if (eqSys->isUsedCorrectionPotential() == 1) iEND = UBZ+1;
      
    } else
    {
      FORT_SOURCEAXISYMMETRIC_MF( CHF_FRA(a_S),
                                  CHF_CONST_FRA(a_W),
                                  CHF_CONST_INT(iRhoN),
                                  CHF_CONST_INT(m_nFluids),
                                  CHF_CONST_REAL(a_dt),
                                  CHF_CONST_REAL(dx),
                                  CHF_BOX(a_box) );
      iEND = eqSys->densityIndexCons(m_nFluids-1)+UENG;
    }
        
    if (eqSys->isUsedCorrectionPotential() == 1)
    {
      FORT_SOURCEAXISYMMETRICDEDNER( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_REAL(dx),
                               CHF_BOX(a_box) );            
    }

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

  }

  if (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical)
  {
    Real dx = m_csh->dx(0,m_level);
    if (m_nFluids==1)
    {
      FORT_SOURCEAXISYMMETRIC_CYL( CHF_FRA(a_S),
                                   CHF_CONST_FRA(a_W),
                                   CHF_CONST_REAL(a_dt),
                                   CHF_CONST_REAL(dx),
                                   CHF_BOX(a_box) );
    } else
    {
      FORT_SOURCEAXISYMMETRIC_CYL_MF( CHF_FRA(a_S),
                                      CHF_CONST_FRA(a_W),
                                      CHF_CONST_INT(iRhoN),
                                      CHF_CONST_INT(m_nFluids),
                                      CHF_CONST_REAL(a_dt),
                                      CHF_CONST_REAL(dx),
                                      CHF_BOX(a_box) );
    }

    iBGN   = UMOMY;
    iEND   = (m_nFluids < 2 ? UMOMY :  UMOMY + eqSys->densityIndexCons(m_nFluids-1)); // UMOMY1+(m_nFluids-2)*(URHO2-URHO1)

  }
  if (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym)
  {
    FORT_SOURCEAXISYMMETRIC_POLAR( CHF_FRA(a_S),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_INT(m_level),
                             CHF_BOX(a_box) );

    iBGN = URHO;
    iEND = UBX;

    for (int iFluid = 1; iFluid < m_nFluids; ++iFluid)
    {
      int startRho = eqSys->densityIndexCons(iFluid);

      FORT_SOURCEAXISYMMETRIC_POLAR_E( CHF_FRA(a_S),
                             CHF_CONST_FRA(a_W),
                             CHF_CONST_REAL(a_dt),
                             CHF_CONST_INT(startRho),
                             CHF_CONST_INT(startRho),
                             CHF_CONST_INT(m_level),
                             CHF_BOX(a_box) );
    }
    if (m_nFluids > 1) iEND = eqSys->densityIndexCons(m_nFluids-1)+UENG;

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );

  }
  
  if (eqSys->isUsedCorrectionPotential() == 1)
  {
    Real dx = m_csh->dx(0,m_level);
    int iRho = eqSys->densityIndexCons(0);
    /*FORT_DEDNERSOURCETERMS(
         CHF_FRA(a_S),
         CHF_CONST_FRA(a_W),
         CHF_CONST_FRA1(divB,0),
         CHF_CONST_INT(iRho),
         CHF_CONST_REAL(dx),     
         CHF_CONST_REAL(a_dt),     
         CHF_BOX(a_box));
    int iBgn = UENG;
    int iEnd = UBZ;
    FORT_ADDSOURCES( CHF_FRA(a_U),
                 CHF_CONST_FRA(a_S),
                 CHF_CONST_INT(iBgn),
                 CHF_CONST_INT(iEnd),
                 CHF_BOX(a_box) );*/
                 
    FORT_PSIPARABOLIC(
         CHF_FRA(a_U),         
         CHF_CONST_INT(iRho),
         CHF_CONST_REAL(dx),     
         CHF_CONST_REAL(a_dt),     
         CHF_BOX(a_box));
  }

//take out for ugradp - needed for pdivu source term for PUIs
  FArrayBox divU( a_box, 1 );
  computeDivU( divU, a_Un, a_W, a_box);
//end take out for ugradp

  eqSys->explicitSource( a_U, a_S, a_W, a_REG, divU, a_dt, m_level, a_box );
//old equation system source terms below, does not have ugradp for PUIs
//also does not pass region
//  eqSys->explicitSource( a_U, a_S, a_W, a_dt, m_level, a_box );
  m_PhPr->explicitSource( a_U, a_S, a_W, a_REG, a_dt, a_box );

  SourceCalculator * pSoCal = m_PhPr->getSourceCalculator();
  if( pSoCal != NULL )
  {
    Real dx = (m_csh->constStep(0)==true ? m_csh->dx(0,m_level) : 1.0);
    pSoCal->addExternalSources( a_U, a_SOut, a_W, a_dt, dx, a_box );
  }

  if( m_verbosity >= 3 )
  {
    pout() << "addExplicitSources:out" << endl;
  }
}


void PatchMHDMF::postprocessing(       FArrayBox & a_U,
                                      const FArrayBox & a_Uold,
                                      const Real      & a_dt,
                                      const Box       & a_box)
{
  if (m_nFluids==1)
  {
    FORT_POSTPROCESSING( CHF_FRA(a_U),
                          CHF_CONST_FRA(a_Uold),
                          CHF_CONST_REAL(a_dt),
                          CHF_BOX(a_box) );
  } else
  {
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    int iRhoN  = eqSys->densityIndexCons(1);
    FORT_POSTPROCESSING_MF( CHF_FRA(a_U),
                          CHF_CONST_FRA(a_Uold),
                          CHF_CONST_INT(iRhoN),
                          CHF_CONST_INT(m_nFluids),
                          CHF_CONST_REAL(a_dt),
                          CHF_BOX(a_box) );

  }
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

  PickupIons* pPI = m_eqSys->getPickupIons();

  if( pPI != NULL )
  {
    if( pPI->modelID() == PickupIons::PI_TWO_EQNS )
    {
      int iRhoPI = pPI->primInterval().begin();
      FORT_POSTPROCESSING_PI( CHF_FRA(a_U),
                              CHF_CONST_INT(iRhoPI),
                              CHF_CONST_INT(iRhoPI),
                              CHF_BOX(a_box));
    }
  }
}

                     // Things to do after advancing this level by one time step
void PatchMHDMF::postTimeStep( void )
{
}

void PatchMHDMF::preprocessing ( const FArrayBox    & a_W,
                                            FArrayBox    & a_S,
                                            BaseFab<int> & a_R,
                                      const Box          & a_box)
{
  if (m_nFluids >= 2)
  {
    MultiFluidProblem * pMFP = (MultiFluidProblem *)(m_PhPr);

    pMFP->defineRegions( a_W, a_S, a_R, a_box);
  }
}

       // Calculate magnetic field in cell centers using electric field on edges
void PatchMHDMF::recalculateMagneticField(       FArrayBox & a_U,
                                           const FArrayBox & a_Uold,
                                           const EdgeBox   & a_E,
                                           const Real      & a_dt,
                                           const Box       & a_box )
{
  int iBx = UBX;

  Real dx = m_csh->dx(0,m_level);

  const FArrayBox& EX = a_E[0];
  const FArrayBox& EY = a_E[1];
  const FArrayBox& EZ = a_E[2];

  FORT_RECALCULATEB( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_Uold),
                     CHF_CONST_FRA1(EX,0),
                     CHF_CONST_FRA1(EY,0),
                     CHF_CONST_FRA1(EZ,0),
                     CHF_CONST_INT(iBx),
                     CHF_CONST_REAL(dx),
                     CHF_CONST_REAL(a_dt),
                     CHF_BOX(a_box) );
}

                                            // Calculate electric field on edges
void PatchMHDMF :: calculateElectricField( const FArrayBox & a_W,
                                           const FArrayBox & a_Wold,
                                           const FArrayBox   a_F[CH_SPACEDIM],
                                                 EdgeBox   & a_E,
                                           const Real      & a_dt,
                                           const Box       & a_box )
{
  switch( m_DivergenceCleaning ){
    case DC_CT_BS:
      {
        int iBx = UBX;

        const FArrayBox& FX = a_F[0];
        const FArrayBox& FY = a_F[1];
        const FArrayBox& FZ = a_F[SpaceDim-1];

        FArrayBox& EX = a_E[0];
        FArrayBox& EY = a_E[1];
        FArrayBox& EZ = a_E[2];

        FORT_ELECTRICFIELDBS( CHF_CONST_FRA(FX),
                              CHF_CONST_FRA(FY),
                              CHF_CONST_FRA(FZ),
                              CHF_FRA1(EX,0),
                              CHF_FRA1(EY,0),
                              CHF_FRA1(EZ,0),
                              CHF_CONST_INT(iBx),
                              CHF_BOX(a_box) );
      } break;
    case DC_CT_GS0:
      {
        int iUx = WVELX;
        int iBx = UBX;

        const FArrayBox& FX = a_F[0];
        const FArrayBox& FY = a_F[1];
        const FArrayBox& FZ = a_F[SpaceDim-1];

        FArrayBox& EX = a_E[0];
        FArrayBox& EY = a_E[1];
        FArrayBox& EZ = a_E[2];

        FORT_ELECTRICFIELDGS0( CHF_CONST_FRA(a_W),
                               CHF_CONST_FRA(FX),
                               CHF_CONST_FRA(FY),
                               CHF_CONST_FRA(FZ),
                               CHF_FRA1(EX,0),
                               CHF_FRA1(EY,0),
                               CHF_FRA1(EZ,0),
                               CHF_CONST_INT(iUx),
                               CHF_CONST_INT(iBx),
                               CHF_BOX(a_box) );
      } break;
    case DC_CT_GS1:
      {
        int iUx = WVELX;
        int iBx = UBX;

        const FArrayBox& FX = a_F[0];
        const FArrayBox& FY = a_F[1];
        const FArrayBox& FZ = a_F[SpaceDim-1];

        FArrayBox& EX = a_E[0];
        FArrayBox& EY = a_E[1];
        FArrayBox& EZ = a_E[2];

        FORT_ELECTRICFIELDGS1( CHF_CONST_FRA(a_W),
                               CHF_CONST_FRA(FX),
                               CHF_CONST_FRA(FY),
                               CHF_CONST_FRA(FZ),
                               CHF_FRA1(EX,0),
                               CHF_FRA1(EY,0),
                               CHF_FRA1(EZ,0),
                               CHF_CONST_INT(iUx),
                               CHF_CONST_INT(iBx),
                               CHF_BOX(a_box) );
      } break;
  }
}
