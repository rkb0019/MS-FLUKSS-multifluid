#include "BoxIterator.H"
#include "IntVectSet.H"

#include "HelioTILTProblem.H"
#include "HelioTILTF_F.H"
#include "HeliosphericF_F.H"

#include "MHDAMDefs.H"

#include "LGintegrator.H"
#include "DebugF_F.H"
#include "EqSysMHDMF.H"
#include "RiemannSolver.H"
#include "TMBreechEtAl2008.H"
#include "HeliosphericTurbF_F.H"
#include "PITwoEquations.H"
#include "HeliosphericPIF_F.H"
#include "HeliosphericPlasmaBCF_F.H"

// Null constructor
HelioTILTProblem::HelioTILTProblem()
 : HeliosphericProblem()
{  
}

HelioTILTProblem::~HelioTILTProblem()
{  
}


// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* HelioTILTProblem::new_PhysProblem()
{
  HelioTILTProblem* retval = new HelioTILTProblem();
  
  retval->copy_PhysProblem(this);  

  return static_cast<PhysProblem*>(retval);

}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void HelioTILTProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  const HelioTILTProblem* PP = dynamic_cast<const HelioTILTProblem*>(a_PP);
  
  if (PP == NULL) MayDay::Error("HelioTILTProblem::copy_PhysProblem. Wrong argument");
  
  HeliosphericProblem::copy_PhysProblem(a_PP);
      
    
 
}

void HelioTILTProblem::defineMesh(const ProblemDomain & a_prob_domain,
                                     const Vector<Real>  & a_domainBox)
{
  HeliosphericProblem::defineMesh(a_prob_domain, a_domainBox);
 
      
}


                                                    // Set up initial conditions
void HelioTILTProblem::initialize( LevelData<FArrayBox>& a_U )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? dx = m_csh->dx(0,m_level) : 1.0);

  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);  
  int iHCS = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);

                                          // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
                                                     // Storage for current grid
    FArrayBox& U = a_U[dit()];

                                                          // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;

                                        // Set up initial condition in this grid
    if( m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian )
    {
      // FORT_HELIOTILTINIT( CHF_FRA(U),
      //                 CHF_CONST_REAL(dx),
      //                 CHF_CONST_INT(iRhoN),
      //                 CHF_CONST_INT(fluids),
      //                 CHF_CONST_INT(iHCS),
      //                 CHF_BOX(uBox) );

        FORT_HELIOINIT( CHF_FRA(U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_BOX(uBox) );
                      
      if (m_subproblem == HPBC_SOLARCYCLE)
      {
        Real t = 0.0;
        FORT_HELIOREINIT_CYCLE( CHF_FRA(U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(t),
                      CHF_CONST_REAL(m_initR),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(uBox) );
      }

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOINIT_TM( CHF_FRA(U),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(iRhoZ2),
                               CHF_BOX(uBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->consInterval().begin();
            FORT_HELIOINIT_PI( CHF_FRA(U),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(iRhoPI),
                               CHF_BOX(uBox) );
          }
        }
      }
    }
    
    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)  
    {
      // FORT_HELIOTILTINITSPHERICAL( CHF_FRA(U),                    
      //               CHF_CONST_INT(iRhoN),
      //               CHF_CONST_INT(fluids),
      //               CHF_CONST_INT(iHCS),
      //               CHF_CONST_INT(m_level),
      //               CHF_BOX(uBox) );

        FORT_HELIOINITSPHERICAL( CHF_FRA(U),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iHCSb),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );

//why is this here?
/*        FORT_HELIOINITPLASMASPHERICAL( CHF_FRA(U),
                      CHF_CONST_INT(m_subproblem),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iHCSb),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );*/

      if (m_subproblem == HPBC_SOLARCYCLE)
      {
        Real t = 0.0;
        FORT_HELIOREINIT_CYCLE_SPHERICAL( CHF_FRA(U),
                      CHF_CONST_REAL(t),
                      CHF_CONST_REAL(m_initR),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );
      }

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOINITSPHERICAL_TM( CHF_FRA(U),
                                        CHF_CONST_INT(iRhoZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_BOX(uBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->consInterval().begin();
            FORT_HELIOINITSPHERICAL_PI( CHF_FRA(U),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_BOX(uBox) );
          }
        }
      }
    }
  }
}


                                                             // Fill ghost cells
void HelioTILTProblem::fillGhostCells(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;

  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  
  int fluids = nFluids();
  
  int iHCS = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    Real dx = m_csh->dx(0,m_level);
    
    int indW = WBox.bigEnd( a_dir );
    int indD =  tmp.bigEnd( a_dir );

    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                           // Fill the ghost cells
      FORT_HELIOTILTGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_CONST_INT(iHCS),                    
                    CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2  = pTurbMod->primInterval().begin();
            FORT_HELIOGS_TM( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOGS_PI( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }
    }

    indW     = WBox.smallEnd( a_dir );
    indD     =  tmp.smallEnd( a_dir );

    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indW, indD - indW );

                                                           // Fill the ghost cells
      FORT_HELIOTILTGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_CONST_INT(iHCS),                    
                    CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2  = pTurbMod->primInterval().begin();
            FORT_HELIOGS_TM( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOGS_PI( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }
    }
  }

  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {
    int indW,indD,nGS, jsize;    
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain;  
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    nGS  = abs(indW-indD);

    jsize = m_domain.domainBox().size(1);

    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
      
      if( a_dir == 0 )
      {
        // FORT_HELIOTILTGSSPHERICAL(
        //     CHF_FRA(a_W),
        //     CHF_CONST_INT(sign),
        //     CHF_CONST_INT(a_dir),
        //     CHF_CONST_INT(iRhoN),
        //     CHF_CONST_INT(fluids),
        //     CHF_CONST_INT(iHCS),
        //     CHF_CONST_INT(m_level),
        //     CHF_CONST_REAL(a_time),
        //     CHF_BOX(boundaryBox));
          FORT_HELIOGSSPHERICAL(
            CHF_FRA(a_W),
            CHF_CONST_FRA(a_U),
            CHF_CONST_INT(sign),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(jsize),
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),
            CHF_CONST_INT(m_ls_indices.m_iHCS),
            CHF_CONST_INT(m_ls_indices.m_iHCSb),
            CHF_CONST_INT(m_ls_indices.m_iRegTr),
            CHF_CONST_INT(m_level),
            CHF_CONST_REAL(a_time),
            CHF_BOX(boundaryBox));

        // FORT_HELIOGSNEUTRALS(
        //   CHF_FRA(a_W),           
        //   CHF_CONST_FRA(a_U),           
        //   CHF_CONST_INT(sign),          
        //   CHF_CONST_INT(iRhoN),
        //   CHF_CONST_INT(fluids),          
        //   CHF_CONST_INT(m_level),
        //   CHF_CONST_REAL(a_time),
        //   CHF_BOX(boundaryBox));            

        if( m_eqSys != NULL )
        {
          TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
          if( pTurbMod != NULL )
          {
            if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
            {
              int iZ2  = pTurbMod->primInterval().begin();
              FORT_HELIOGSSPHERICAL_TM( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_HELIOGSSPHERICAL_PI( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }
        }
      }

      if (a_dir == 2)
      {
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
      }
    }

    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
  
      if (a_dir == 0)
      {
        // FORT_HELIOTILTGSSPHERICAL(
        //   CHF_FRA(a_W),
        //   CHF_CONST_INT(sign),
        //   CHF_CONST_INT(a_dir),
        //   CHF_CONST_INT(iRhoN),
        //   CHF_CONST_INT(fluids),
        //   CHF_CONST_INT(iHCS),
        //   CHF_CONST_INT(m_level),
        //   CHF_CONST_REAL(a_time),
        //   CHF_BOX(boundaryBox));
          FORT_HELIOGSSPHERICAL(
            CHF_FRA(a_W),
            CHF_CONST_FRA(a_U),
            CHF_CONST_INT(sign),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(jsize),
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),
            CHF_CONST_INT(m_ls_indices.m_iHCS),
            CHF_CONST_INT(m_ls_indices.m_iHCSb),
            CHF_CONST_INT(m_ls_indices.m_iRegTr),
            CHF_CONST_INT(m_level),
            CHF_CONST_REAL(a_time),
            CHF_BOX(boundaryBox));

//        FORT_HELIOGSNEUTRALS(
//          CHF_FRA(a_W),           
//          CHF_CONST_FRA(a_U),           
//          CHF_CONST_INT(sign),          
//          CHF_CONST_INT(iRhoN),
//          CHF_CONST_INT(fluids),          
//          CHF_CONST_INT(m_level),
//          CHF_CONST_REAL(a_time),
//          CHF_BOX(boundaryBox));                      

        if( m_eqSys != NULL )
        {
          TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
          if( pTurbMod != NULL )
          {
            if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
            {
              int iZ2  = pTurbMod->primInterval().begin();
              FORT_HELIOGSSPHERICAL_TM( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_HELIOGSSPHERICAL_PI( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }
        }
      }
      if (a_dir == 2)
      {
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
      }
    }
  }
}

                                                          // Set boundary fluxes
void HelioTILTProblem::fluxBC(       FArrayBox&      a_F,
                                        FArrayBox&      a_Bn,
                                  const FArrayBox&      a_WMinus,
                                  const FArrayBox&      a_WPlus,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side,
                                  const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
  {
    HeliosphericProblem::fluxBC( a_F, a_Bn, a_WMinus, a_WPlus, a_dir, a_side,  a_time);
    return;
  }

  int sign;
  Box FBox = a_F.box();
  Box tmp = FBox;
  
  Real dx = m_csh->dx(0,m_level);
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);  

                            // Determine which side and thus shifting directions
  if( a_side == Side::Lo )
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }

  tmp.shiftHalf(a_dir,sign);

                                 // Is there a domain boundary next to this grid
  if( !m_domain.contains( tmp ) )
  {
    tmp &= m_domain;

    Box boundaryBox;

                          // Find the strip of cells next to the domain boundary
    if( a_side == Side::Lo )
    {
      boundaryBox = bdryLo(tmp,a_dir);
    }
    else
    {
      boundaryBox = bdryHi(tmp,a_dir);
    }
    
    m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );                        

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
    FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);

    int fluids = nFluids();
                                                      // Set the boundary fluxes
    // FORT_HELIOTILTBC( CHF_FRA(a_F),
    //               CHF_FRA1(a_Bn,0),
    //               CHF_CONST_FRA(shiftWLeft),
    //               CHF_CONST_FRA(shiftWRight),
    //               CHF_CONST_INT(sign),
    //               CHF_CONST_REAL(dx),
    //               CHF_CONST_INT(a_dir),
    //               CHF_CONST_INT(iRhoN),
    //               CHF_CONST_INT(fluids),
    //               CHF_BOX(boundaryBox) );
        FORT_HELIOBC( CHF_FRA(a_F),
                  CHF_FRA1(a_Bn,0),
                  CHF_CONST_FRA(shiftWLeft),
                  CHF_CONST_FRA(shiftWRight),
                  CHF_CONST_INT(sign),
                  CHF_CONST_REAL(dx),
                  CHF_CONST_INT(a_dir),
                  CHF_CONST_INT(iRhoN),
                  CHF_CONST_INT(fluids),
                  CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2    = pTurbMod->primInterval().begin();
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOBC_TM( CHF_FRA(a_F),
                             CHF_CONST_FRA(shiftWLeft),
                             CHF_CONST_FRA(shiftWRight),
                             CHF_CONST_INT(sign),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_CONST_INT(iRhoZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOBC_PI( CHF_FRA(a_F),
                             CHF_CONST_FRA(shiftWLeft),
                             CHF_CONST_FRA(shiftWRight),
                             CHF_CONST_INT(sign),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }

    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
  }
}


//                              Creates tagged cells for dynamic mesh refinement
void HelioTILTProblem::tagCells( const FArrayBox&  a_U,
                              const Box&        a_box,
                              IntVectSet& a_tags )
{  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)                           
  {
    const Box& b = a_box;

    BoxIterator bit( b );
    
    Real dx = m_csh->dx(0,m_level);

    for( bit.begin(); bit.ok(); ++bit )
    {
      const IntVect& iv = bit();

      Real D_DECL(x,y,z);

      D_TERM(
        x = (iv[0] + 0.5)*dx - m_sunXYZ[0];,
        y = (iv[1] + 0.5)*dx - m_sunXYZ[1];,
        z = (iv[2] + 0.5)*dx - m_sunXYZ[2];
      );

      Real dist2 = D_TERM( x*x, + y*y, + z*z );

      //if (( 0.9*m_initR*m_initR < dist2 ) && ( dist2 < 1.1*m_initR*m_initR ))
      if (( dist2 > m_initR*m_initR  ) && (x > 0.0))
      {
        a_tags |= iv;
      }
    }
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)  
  {
    BoxIterator bit( a_box ); RealVect cv_coord,cs_coord;

    for( bit.begin(); bit.ok(); ++bit )
    {
      const IntVect& iv = bit();
      m_csh->getCellCenter(cv_coord,iv,m_level);
      m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);

      if (m_level == 0)
      {
        if ( (cv_coord[0]>40.0) && (cv_coord[0]<50.0) ) a_tags |= iv;
      }       
/*
      Real theta1 = (90.0-(m_sunTILT + 3.0))*d_PI_180;
      Real theta2 = (90.0+(m_sunTILT + 3.0))*d_PI_180;
      
      Real phi    = 15.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;
      
      if (m_level == 0)
      {
        if ( (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && (cv_coord[0]<160.0) &&
            ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
        
      }    
      if (m_level == 1)
      {              
        phi    = 45.0;
        phi1   = (90.0  - phi)*d_PI_180;
        phi2   = (270.0 + phi)*d_PI_180;
        
        if ((cv_coord[0]>m_initR) && (cv_coord[0]<160.0) &&
            (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
            ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
        
      }
      if (m_level == 2)
      {
        if ( (cv_coord[0]>60.0) && (cv_coord[0]<150.0) &&      
             (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
             ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
        
      }*/
      
 /*     if (m_level == 0)
      {
        if ((cv_coord[0]>40.0) &&        
            (cv_coord[2]>d_PI_4) && (cv_coord[2]<d_3PI_4) && 
            (cs_coord[0]>0.0) ) a_tags |= iv;
        
      }    
      if (m_level == 1)
      {
        if ((cv_coord[0]>65.0) &&        
            (cv_coord[2]>d_PI_4) && (cv_coord[2]<d_3PI_4) && 
            (cs_coord[0]>0.0) ) a_tags |= iv;
        
      }
      if (m_level == 2)
      {
        if ((cv_coord[0]>75.0) &&        
            ( ((cv_coord[1]>0.0) && (cv_coord[1]<d_PI_4)) || (cv_coord[1]>7.0/4.0*d_PI)  ) &&
            (cv_coord[2]>d_PI_4) && (cv_coord[2]<d_3PI_4) && 
            (cs_coord[0]>0.0) ) a_tags |= iv;
        
      } */
    }
  }
  
}

//                             Check geometrical limitations for grid adaptation
void HelioTILTProblem::lockedCellsRegrid( BaseFab<int> & a_flag,
                                       const FArrayBox&  a_U,
                                       const Box&     a_box)
{
  a_flag.setVal(0);
        
}                            

//                            Return boundary condition flags for all boundaries
void HelioTILTProblem::getBCFlags( eBoundaryConditions leftBC,
                                      eBoundaryConditions rightBC,
                                      eBoundaryConditions bottomBC,
                                      eBoundaryConditions topBC,
                                      eBoundaryConditions frontBC,
                                      eBoundaryConditions behindBC )
{
  leftBC   = BC_Continuous;
  rightBC  = BC_Fixed;
  bottomBC = BC_Continuous;
  topBC    = BC_Continuous;
  frontBC  = BC_Continuous;
  behindBC = BC_Continuous;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions HelioTILTProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    return ( a_sd == Side::Lo ) ? BC_Continuous : BC_Fixed;
  } else {
    return BC_Continuous;
  }
}

