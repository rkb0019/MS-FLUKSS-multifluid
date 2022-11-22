#include "LayoutIterator.H"
#include <PiecewiseLinearFillPatch.H> 
#include "FORT_PROTO.H"

#include "KineticSources3D.H"
#include "KineticSources3DF_F.H"
#include "AMRLevelIdealMHD.H"
#include "PatchIdealMHDF_F.H"
#include "LGintegrator.H"
#include "KS2DIntergrator.H"
#include "EosCommon.H"
#include "DebugF_F.H"
#include "TecplotIO.H"
#include "EqSysMHDMF.H"
#include "SWLISMF_F.H"


  
extern "C"
{
// Prototype for Fortran procedure mc_init 

#define FORT_MC3D_INIT_CALL_1 FORTRAN_NAME( MC3D_INIT_CALL_1 ,mc3d_init_call_1 )
void FORT_MC3D_INIT_CALL_1
(  
  int*  const au_nlevels     
);


#define FORT_MC3D_INIT_CALL_2 FORTRAN_NAME( MC3D_INIT_CALL_2 ,mc3d_init_call_2 )
void FORT_MC3D_INIT_CALL_2
(
  Real* const au_LISM_nH,
  Real* const au_LISM_vH,
  Real* const au_LISM_TH,   
  Real* const au_total_neutrals, 
  int*  const au_restart,
  char* const au_load_file      
);

#define FORT_SETUP_LEVEL_3D FORTRAN_NAME( SETUP_LEVEL_3D ,setup_level_3d )
void FORT_SETUP_LEVEL_3D
(
  int*  const au_level,
  int*  const au_ref_ratio,
  Real* const au_Lo,
  Real* const au_Hi,
  Real* const au_dx,
  int*  const au_size,  
  Real* const au_source_m, 
  Real* const au_source_px, 
  Real* const au_source_py, 
  Real* const au_source_pz, 
  Real* const au_source_mvsq,
  Real* const au_plasma,
  Real* const au_nchex, 
  Real* const au_H_dens, 
  Real* const au_H_ux, 
  Real* const au_H_uy, 
  Real* const au_H_uz, 
  Real* const au_H_temp
);



// Prototype for Fortran procedure mc_neutrals
#define FORT_MC_NEUTRALS_3D FORTRAN_NAME( RUN_MC_NEUTRALS_3D, run_mc_neutrals_3d)
void FORT_MC_NEUTRALS_3D(
    Real* const run_time    
);


#define FORT_MC_OUTPUT_RAW_3D FORTRAN_NAME( RUN_MC_OUTPUT_RAW_3D ,run_mc_output_raw_3d )
void FORT_MC_OUTPUT_RAW_3D
(    
  const char* const chkfile     
);

#define FORT_FINALIZE_KINETIC3D_F FORTRAN_NAME( FINALIZE_KINETIC3D_F ,finalize_kinetic3d_f )
void FORT_FINALIZE_KINETIC3D_F
(    
);


}



                                                                   // Costructor
KineticSources3D::KineticSources3D()
{
  m_totalNeutrals = 0.0;
    
  m_runtimeYears = -1.0;        
  m_timeFactor = 1.0;                 
  
  m_timeIntervalYears = -1.0;   
  m_timeIntervalMHDAM = -1.0;   
  m_nextTimeMHDAM = 0.0;   
  
  m_stepInterval = -1;   
  
  m_firstCall = true;
  
  m_chk_prefix = "n";
  
  m_cur_step = -1;
}
                                                                   // Destructor
KineticSources3D::~KineticSources3D()
{
  for (int ilev=0;ilev<m_ngrids;ilev++)
  {
    LevelData<FArrayBox>* LD = m_plasmaDataLD[ilev];
    if (LD!=NULL) delete LD;
    
    FArrayBox* FAB = m_plasmaDataFAB[ilev];    
    if (procID() > 0) 
    if (FAB!=NULL) delete FAB;  
    
    FAB = m_sourceTerms[ilev];  
    if (FAB!=NULL) delete FAB;    
    
    FAB = m_neutralData[ilev];  
    if (FAB!=NULL) delete FAB;            
  }
  
  FORT_FINALIZE_KINETIC3D_F();
}
                              // Factory method - this object is its own factory
SourceCalculator * KineticSources3D :: new_SourceCalculator( void )
{
  SourceCalculator* retval = new KineticSources3D();

  return retval;
}

                                                             // Input parameters
void KineticSources3D :: input( ParmParse & parser, int verbosity )
{     
  int i; 
  
  parser.query( "lismN",   m_lismN );
  parser.query( "lismV",   m_lismV );
  parser.query( "lismT",   m_lismT );
  parser.query( "netN",    m_netN  );  
  
  D_TERM(parser.query( "XC",      m_sunXYZ[0]    );,
         parser.query( "YC",      m_sunXYZ[1]    );,
         parser.query( "ZC",      m_sunXYZ[2]    ););  
         
  parser.query( "lismUX",  m_lismUX  );
  parser.query( "lismUY",  m_lismUY  );
  parser.query( "lismUZ",  m_lismUZ  );
         
  parser.query("chk_prefix",m_chk_prefix);    
  
  
  int maxLevel = 0;
  parser.query("max_level",maxLevel);
  
  Real domainLength = 1.0;
  parser.query("domain_length",domainLength);
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim);
  for (int i = 0; i < SpaceDim; ++i) numCells[i]=0;
  parser.queryarr("num_cells",numCells,0,SpaceDim);
  
  // For now this is Jacob's requirement
  //if (numCells[0] != 2*numCells[1]) 
  //  MayDay::Error("KineticSources3D::input, num_cells[0] must be equal 2*num_cells[1]");
  
  Real dx = domainLength/numCells[0];
  for (i = 0; i < SpaceDim; ++i) m_sunXYZ[i] = ((int)floor(m_sunXYZ[i]/dx+0.5))*dx;
  
  
  D_TERM(
    m_sunIJK[0] = (int)floor(m_sunXYZ[0]/dx+0.5);,
    m_sunIJK[1] = (int)floor(m_sunXYZ[1]/dx+0.5);,
    m_sunIJK[2] = (int)floor(m_sunXYZ[2]/dx+0.5););  

  D_TERM(  
      CH_assert(fabs(m_sunIJK[0]*dx-m_sunXYZ[0])<1e-6);,
      CH_assert(fabs(m_sunIJK[1]*dx-m_sunXYZ[1])<1e-6);,
      CH_assert(fabs(m_sunIJK[2]*dx-m_sunXYZ[2])<1e-6););
  
  
  char buf[20];std::vector<Real> tmpVect;    
  parser.query( "kinetic_ngrids",      m_ngrids );
  if (maxLevel + 1 < m_ngrids)
    MayDay::Error("KineticSources3D::input, number of additional levels should be greater than number of kinetic levels");
  
  m_ngrids++;// Total number of grids, base level added
  
  m_gridsBoundaries.resize(m_ngrids);  RealVect rv;
  for (i=1;i<m_ngrids;i++)
  {
    sprintf(buf, "kinetic_grid%ibb", i); // We read data specified relatively the Sun position
    parser.queryarr( buf, tmpVect, 0, 2*CH_SPACEDIM);    
    D_TERM(
      rv[0] = tmpVect[0];,rv[1] = tmpVect[1];,rv[2] = tmpVect[2];);
    rv+=m_sunXYZ;    
    m_gridsBoundaries[i].m_smallEnd = rv;    
    
    if (i>1) // Base level is specified later, in initExternalSC
    if (D_TERM( (m_gridsBoundaries[i].m_smallEnd[0] <= m_gridsBoundaries[i-1].m_smallEnd[0])  ,
                || (m_gridsBoundaries[i].m_smallEnd[1] <= m_gridsBoundaries[i-1].m_smallEnd[1])  ,
                || (m_gridsBoundaries[i].m_smallEnd[2] <= m_gridsBoundaries[i-1].m_smallEnd[2]) ))
      MayDay::Error("Check kinetic boundaries. Higher level box must lie inside of the lowere level box completely. Check lower boundaries");
      
    D_TERM(
      rv[0] = tmpVect[0+CH_SPACEDIM];,rv[1] = tmpVect[1+CH_SPACEDIM];,rv[2] = tmpVect[2+CH_SPACEDIM];);
    rv+=m_sunXYZ;    
    m_gridsBoundaries[i].m_bigEnd = rv;    
    
    if (i>1) // Base level is specified later, in initExternalSC
    if (D_TERM( (m_gridsBoundaries[i].m_bigEnd[0] >= m_gridsBoundaries[i-1].m_bigEnd[0])  ,
                || (m_gridsBoundaries[i].m_bigEnd[1] >= m_gridsBoundaries[i-1].m_bigEnd[1])  ,
                || (m_gridsBoundaries[i].m_bigEnd[2] >= m_gridsBoundaries[i-1].m_bigEnd[2]) ))
      MayDay::Error("Check kinetic boundaries. Higher level box must lie inside of the lowere level box completely. Check upper boundaries");
  }  
          
  m_runtimeYears = -1.0;
  parser.query("kinetic_runtime",  m_runtimeYears);    
  
  m_initialRuntimeYears = -1.0;
  parser.query("kinetic_initial_runtime",  m_initialRuntimeYears);    
        
  m_timeIntervalYears = -1.0;     
  parser.query("kinetic_time_interval",  m_timeIntervalYears);    
  m_timeFactor        = (eos_AU/m_lismV)/(60.0*60.0*24.0*365.0);   
  m_timeIntervalMHDAM = m_timeIntervalYears/m_timeFactor;     // Transform run time to our units.  
  
  m_stepInterval = -1;
  parser.query("kinetic_step_interval",  m_stepInterval);        
    
  m_totalNeutrals = 20000;
  parser.query( "total_neutrals",  m_totalNeutrals);
  
  if (parser.contains("kinetic_restart_file"))  parser.query("kinetic_restart_file", m_restartFile);    
  
  m_photoionize = 0;
  parser.query( "kinetic_photoionize",  m_photoionize);
  if (m_photoionize!=1) m_photoionize = 0;
  
  m_verbosity = verbosity;
  
                                                             // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << endl << "Source term calculator parameters:"      << endl;        
    pout() << "    kinetic_runtime       = " << m_runtimeYears  << endl;    
    pout() << "    total_neutrals        = " << m_totalNeutrals  << endl;
    pout() << "    kinetic_restart_file  = " << m_restartFile  << endl;
    pout() << "    kinetic_time_interval = " << m_timeIntervalYears  << endl;
    pout() << "    kinetic_step_interval = " << m_stepInterval  << endl;    
    pout() << "    kinetic_time_interval = " << m_timeIntervalMHDAM  << " (MHDAM)" << endl;
    pout() << "    kinetic_photoionize   = " << m_photoionize  << endl;    
  }
  
  m_nextTimeMHDAM = 0.0;
  m_firstCall     = true;
}


// Initialization of external source calculator.
void KineticSources3D :: initExternalSC(Vector<AMRLevel*> a_levels)
{  
  int ilev; AMRLevelIdealMHD* curLevel;  
  
  if (m_ngrids > a_levels.size())
  {
    m_ngrids = a_levels.size();
    if (procID() == 0)
      MayDay::Warning("KineticSources3D: number of kinetic levels ('kinetic_ngrids' parameter) is larger than number of AMR levels");
  }
  
  
  Real dx;
  
  if( m_verbosity >= 3 )
  {
    pout() << "    KineticSources3D::initExternalSC "  << endl;    
  }
  
  Real lismV[3]={m_lismV*m_lismUX,m_lismV*m_lismUY,m_lismV*m_lismUZ};
  
  Real total_neutrals = m_totalNeutrals;
  
  int  restart = (m_restartFile.size() > 0 ? 1 : 0);
  char loadfile[60]={0}; loadfile[0] = 0;
  if (restart == 1) m_restartFile.copy(loadfile,m_restartFile.size());
  
  int verbosity = m_verbosity;
  verbosity = 3;  
  
  FORT_MC3D_INIT_CALL_1(&m_ngrids);
  
  
  Box b; FArrayBox * pFAB;  
  
  int ref_ratio,tmp_size[CH_SPACEDIM];
  
  RealVect rLo,rHi;
  
  m_plasmaDataLD. resize(m_ngrids,NULL);
  m_plasmaDataFAB.resize(m_ngrids,NULL);  
  m_sourceTerms  .resize(m_ngrids,NULL);
  m_neutralData  .resize(m_ngrids,NULL);  
    
  
  ilev = 0;
  curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[ilev]);    
  dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);
  b  = curLevel->problemDomain().domainBox();  
  m_kineticBoxes.push_back(b);
  
  Vector<Box> vmainBox;
  Vector<int> vmainProc;
            
  vmainBox.push_back(b);
  vmainProc.push_back(0);
  DisjointBoxLayout dbMainBox(vmainBox,vmainProc);
  m_plasmaDataLD[ilev] = new LevelData<FArrayBox>(dbMainBox,numVarsForSC());
  
  if (procID() == 0)
  {    
    FArrayBox & tmpFAB = (*m_plasmaDataLD[ilev])[m_plasmaDataLD[ilev]->dataIterator()()];
    pFAB = &tmpFAB;
  } else
  {
    pFAB = new FArrayBox(b,numVarsForSC());
  }
  
  m_plasmaDataFAB[ilev] = pFAB;
        
  pFAB = new FArrayBox(b,SCOMP);
  m_sourceTerms[ilev] = pFAB;
  
  m_gridsBoundaries[ilev].m_smallEnd = RealVect(b.smallEnd())*dx;
  m_gridsBoundaries[ilev].m_bigEnd   = RealVect(b.bigEnd()+IntVect::Unit)  *dx;
  ref_ratio = 1;
  D_TERM(tmp_size[0]=b.size(0);,tmp_size[1]=b.size(1);,tmp_size[2]=b.size(2););
  
  rLo =  m_gridsBoundaries[ilev].m_smallEnd;
  rLo -= m_sunXYZ;
  rHi =  m_gridsBoundaries[ilev].m_bigEnd;
  rHi -= m_sunXYZ;
  
  // Setup m_neutralData, that should have extra cell in each dimension
  Box b_neutrals = b;
  b_neutrals.grow(1);
  pFAB = new FArrayBox(b_neutrals,NCOMP);
  m_neutralData[ilev] = pFAB;
    
    
  FORT_SETUP_LEVEL_3D(&ilev, &ref_ratio, 
    rLo.dataPtr(),
    rHi.dataPtr(),
    &dx, tmp_size,
    m_sourceTerms[ilev]->dataPtr(URHO),
    m_sourceTerms[ilev]->dataPtr(UMOMX),
    m_sourceTerms[ilev]->dataPtr(UMOMY),
    m_sourceTerms[ilev]->dataPtr(UMOMZ),
    m_sourceTerms[ilev]->dataPtr(UENG),
    
    m_plasmaDataFAB[ilev]->dataPtr(),
    
    m_neutralData[ilev]->dataPtr(NCHEX-NCHEX),
    m_neutralData[ilev]->dataPtr(NRHO -NCHEX),
    m_neutralData[ilev]->dataPtr(NVELX-NCHEX),
    m_neutralData[ilev]->dataPtr(NVELY-NCHEX),
    m_neutralData[ilev]->dataPtr(NVELZ-NCHEX),
    m_neutralData[ilev]->dataPtr(NTEMP-NCHEX)   
    );
  
  
  for (ilev=1;ilev<m_ngrids;ilev++)
  {
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[ilev]);  
    dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);
    
    IntVect iLo(D_DECL (    
      (int)floor(m_gridsBoundaries[ilev].m_smallEnd[0]/dx+0.5),
      (int)floor(m_gridsBoundaries[ilev].m_smallEnd[1]/dx+0.5),
      (int)floor(m_gridsBoundaries[ilev].m_smallEnd[2]/dx+0.5) ));      
      
    iLo.coarsen(a_levels[ilev-1]->refRatio());
    iLo.scale  (a_levels[ilev-1]->refRatio());
            
    IntVect iHi(D_DECL (    
      (int)floor(m_gridsBoundaries[ilev].m_bigEnd[0]/dx+0.5),
      (int)floor(m_gridsBoundaries[ilev].m_bigEnd[1]/dx+0.5),
      (int)floor(m_gridsBoundaries[ilev].m_bigEnd[2]/dx+0.5) ));      
      
    iHi.coarsen(a_levels[ilev-1]->refRatio());
    iHi.scale  (a_levels[ilev-1]->refRatio());
    
    m_gridsBoundaries[ilev].m_smallEnd = RealVect(iLo)*dx;
    m_gridsBoundaries[ilev].m_bigEnd   = RealVect(iHi+IntVect::Unit)*dx;
    
    b.define(iLo,iHi);
    D_TERM(tmp_size[0]=b.size(0);,tmp_size[1]=b.size(1);,tmp_size[2]=b.size(2););
    
    
    m_kineticBoxes.push_back(b);
  
    vmainBox.clear();    
            
    vmainBox.push_back(b);  
    DisjointBoxLayout dbMainBox_a(vmainBox,vmainProc);
    m_plasmaDataLD[ilev] = new LevelData<FArrayBox>(dbMainBox_a,numVarsForSC());
  
    if (procID() == 0)
    {
      FArrayBox & tmpFAB = (*m_plasmaDataLD[ilev])[m_plasmaDataLD[ilev]->dataIterator()()];
      pFAB = &tmpFAB;
    } else
    {
      pFAB = new FArrayBox(b,numVarsForSC());
    }
  
    m_plasmaDataFAB[ilev] = pFAB;
    
    // Fill plasma with LISM. 
    // We need this becuase if we restart with modified kinetic levels some plasma data might not be initialized before regridding happens.            
    pFAB->setVal(m_lismN*1e+6,KRHO);
    pFAB->setVal(m_lismV*m_lismUX*1e-2,KVELX);
    pFAB->setVal(m_lismV*m_lismUY*1e-2,KVELY);
    pFAB->setVal(m_lismV*m_lismUZ*1e-2,KVELZ);    
    pFAB->setVal(m_lismT,KTEMP);
    pFAB->setVal(0.1    ,KREG);
    
        
    pFAB = new FArrayBox(b,SCOMP);
    m_sourceTerms[ilev] = pFAB;
    
    ref_ratio = a_levels[ilev-1]->refRatio();
    
    rLo =  m_gridsBoundaries[ilev].m_smallEnd;
    rLo -= m_sunXYZ;
    rHi =  m_gridsBoundaries[ilev].m_bigEnd;
    rHi -= m_sunXYZ;
    
    // Setup m_neutralData, that should have extra cell in each dimension
    b_neutrals = b;
    b_neutrals.grow(1);
    pFAB = new FArrayBox(b_neutrals,NCOMP);
    m_neutralData[ilev] = pFAB;
      
      
    FORT_SETUP_LEVEL_3D(&ilev, &ref_ratio, 
      rLo.dataPtr(),
      rHi.dataPtr(),
      &dx, tmp_size,
      m_sourceTerms[ilev]->dataPtr(URHO),
      m_sourceTerms[ilev]->dataPtr(UMOMX),
      m_sourceTerms[ilev]->dataPtr(UMOMY),
      m_sourceTerms[ilev]->dataPtr(UMOMZ),
      m_sourceTerms[ilev]->dataPtr(UENG),
      
      m_plasmaDataFAB[ilev]->dataPtr(),
      
      m_neutralData[ilev]->dataPtr(NCHEX-NCHEX),
      m_neutralData[ilev]->dataPtr(NRHO -NCHEX),
      m_neutralData[ilev]->dataPtr(NVELX-NCHEX),
      m_neutralData[ilev]->dataPtr(NVELY-NCHEX),
      m_neutralData[ilev]->dataPtr(NVELZ-NCHEX),
      m_neutralData[ilev]->dataPtr(NTEMP-NCHEX)   
      );  
      
   }  
   
   FORT_MC3D_INIT_CALL_2(&m_netN, lismV, &m_lismT, &total_neutrals, &restart,  loadfile);
    
  
}


//                                       Check time for source terms calculation
bool KineticSources3D :: checkTime( Real dTime, int iStep )
{
  //return true;
  
  if (iStep == m_cur_step)
    MayDay::Error("KineticSources3D::checkTime must never be called twice for the same step");
  
  bool bRetCode  = false;
  m_cur_step = iStep;
  
      
  if( m_timeIntervalMHDAM > 0.0 )
  {
    if( dTime >= m_nextTimeMHDAM )
    {
      m_nextTimeMHDAM += m_timeIntervalMHDAM;
      bRetCode     = true;
    }
    if( dTime >= m_nextTimeMHDAM )
    {
      m_nextTimeMHDAM = dTime + m_timeIntervalMHDAM;
      bRetCode     = true;
    }
  }

  if( m_stepInterval > 0 )
  {
    if( iStep - (iStep/m_stepInterval)*m_stepInterval == 0 )
    {
       bRetCode    = true;
    }
 }

  return bRetCode;
}

// number of variables that are needed from  main code for external source calculator
int KineticSources3D::numVarsForSC()
{
  return KCOMP;
}

//                                             Number of calculated source terms
//                                                    and neutral parameters
int KineticSources3D :: numSourceTerms( void )
{
  return SCOMP+ // source terms
         NCOMP; // neutral data
}


void KineticSources3D::CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep)
{ 
  if ((abs(m_verbosity)>=1) && (procID() == 0))
  {
    pout() << "KineticSources3D::CalculateSources: step/time for new source terms" << endl;
  }
  
  timeval time1,time2,time3,time4;
  
  gettimeofday(&time3,NULL);
  
  EquationSystem * eqSys = static_cast<AMRLevelIdealMHD*>(a_levels[0])->getpatchMHDAM()->getPhysProblem()->equationSystem();
  
  int ilev; AMRLevelIdealMHD* curLevel;
  
  for (ilev=0;ilev<m_ngrids;ilev++)
  {
    Vector<Box> vkineticBoxes;
    Vector<int> vkineticProcs;
    
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[ilev]);  
    Real dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);
    
    // Create LevelData that reflects are kinetic region and original distributed mesh
        
    const Box mainBox = m_kineticBoxes[ilev];
                
    const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
    
    Box bpart;
    LayoutIterator lit = levelDomain.layoutIterator();
    
    for (lit.begin();lit.ok();++lit)              
    {            
      const Box& b = levelDomain[lit()];                
      bpart = b & mainBox;
      if (bpart.isEmpty()) continue;
      
      vkineticBoxes.push_back(bpart);
      vkineticProcs.push_back(levelDomain.procID(lit.i()));                  
    }
        
    DisjointBoxLayout    dbKinetic(vkineticBoxes,vkineticProcs);    
    LevelData<FArrayBox> ldKinetic(dbKinetic,numVarsForSC());
    
          
    DataIterator dit_d = levelDomain.dataIterator();
    DataIterator dit_k = dbKinetic.dataIterator();
    for (dit_d.begin(); dit_d.ok(); ++dit_d)  
    {        
      const Box& b = levelDomain[dit_d()];
      if (!mainBox.intersects(b)) continue;      
            
      const FArrayBox& UFab = curLevel->getStateNew()[dit_d()];	            
      FArrayBox& WKin = ldKinetic[dit_k()];	            
      
      CH_assert(b.contains(WKin.box()));
      
      bpart = b & mainBox;
      CH_assert(bpart == WKin.box());
      
      FORT_PREPAREDATAFORKINETICSC3D(
        CHF_CONST_FRA(UFab),
        CHF_FRA(WKin),
        CHF_BOX(bpart));
                    
      FArrayBox Wtmp(bpart,eqSys->numPrimitives());
      BaseFab<int>  Reg(bpart,1);
      
      eqSys->stateToPrim(Wtmp,UFab,bpart);         
                
      FORT_DEFINE_REGIONS_2F( CHF_CONST_FRA(Wtmp),
                            CHF_FIA1(Reg,0),
                            CHF_CONST_REAL(dx),
                            CHF_BOX(bpart) );
                            
                            
      BoxIterator bit( bpart ); int r;
      for( bit.begin(); bit.ok(); ++bit )
      {
        const IntVect& iv = bit();
        Reg.getVal(&r,iv,0,1);
        switch (r){
          case 1: WKin.set(iv, KREG, 0.1);break;
          case 2: WKin.set(iv, KREG, 2.1);break;
          case 3: WKin.set(iv, KREG, 3.1);break;
          default:WKin.set(iv, KREG, 0.1); break;
        }
        if (fabs(iv[0]*dx-m_sunXYZ[0])>1000.0) WKin.set(iv, KREG, 0.1); // It is always region 1 in distant tail. We need this to recombine neutrals
      }
      
      ++dit_k;           
    }
    
    ldKinetic.copyTo(ldKinetic.interval(), *(m_plasmaDataLD[ilev]), ldKinetic.interval());
        
    #ifdef CH_MPI        
    int BufferSize;  
    BufferSize = mainBox.numPts()*numVarsForSC();
    MPI_Bcast(m_plasmaDataFAB[ilev]->dataPtr(0),BufferSize,MPI_CH_REAL,0 /* rank of broadcast root */,Chombo_MPI::comm);  
    #endif  
    
        
    if ((0)&&(ilev>0) && (procID() == 0))
    {  
      char file_name[50];  
              
      sprintf(file_name,"lev%i.dat",ilev);
      FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"Z\" \"rho\" \"ux\" \"uy\" \"uz\" \"T\" \"Reg\"");    
      WriteFArrayBoxToTecplotFile(tfile, *(m_plasmaDataFAB[ilev]), m_plasmaDataFAB[ilev]->box(), Interval(0,m_plasmaDataFAB[ilev]->nComp()-1), dx);
      CloseTecplotFile(tfile);
    }
    
  }


   
  Real runtimeYears = m_runtimeYears;    
  int  restart = (m_restartFile.size() > 0 ? 1 : 0);
  if ((restart == 0)&&(m_firstCall)) runtimeYears = m_initialRuntimeYears;
  
  
  if ((abs(m_verbosity)>=2) && (procID() == 0))
  {
    pout() << "FORT_MC_NEUTRALS_3D" << endl;
  }
  
  
  
  gettimeofday(&time1,NULL);
  FORT_MC_NEUTRALS_3D(&runtimeYears);  
  gettimeofday(&time2,NULL);
  
  if (procID() == 0)
  {
    double stepTime = (double)(time2.tv_sec-time1.tv_sec) + 1e-6*((double)(time2.tv_usec-time1.tv_usec));
    char buf[20];
    sprintf(buf,"%.4g",stepTime);
    pout() << "MC_NEUTRALS_3D execution time: " << buf << " sec" << endl;
  }
  
  
  
  m_firstCall = false;
  
  // Output results, test only
  for (ilev=0;ilev<m_ngrids;ilev++)  
  if ((0) && (ilev>0) && (procID() == 0))
  {    
    Real dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);
    
    char file_name[50];                
    sprintf(file_name,"storig_lev%i_proc%i.dat",ilev,procID());
    FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"SC\"\n VARIABLES=\"X\" \"Y\" \"Z\" \"SRHO\" \"SMOMX\" \"SMOMY\" \"SMOMZ\" \"SENG\"");    
    WriteFArrayBoxToTecplotFile(tfile, *(m_sourceTerms[ilev]), m_sourceTerms[ilev]->box(), Interval(0,m_sourceTerms[ilev]->nComp()-1), dx);
    CloseTecplotFile(tfile);
        
    sprintf(file_name,"neutrals_lev%i_proc%i.dat",ilev,procID());
    tfile = OpenTecplotFile(file_name,"TITLE = \"SC\"\n VARIABLES=\"X\" \"Y\" \"Z\" \"NCHEX\" \"NRHO\" \"NVELX\" \"NVELY\" \"NVELZ\" \"NTEMP\"");    
    WriteFArrayBoxToTecplotFile(tfile, *(m_neutralData[ilev]), m_neutralData[ilev]->box(), Interval(0,m_neutralData[ilev]->nComp()-1), dx);
    CloseTecplotFile(tfile);
  }
  
  
  
  // Copy to data to Level 0
  curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[0]);  
  LevelData<FArrayBox>& LSCData0 = curLevel->getSCData();   
  const DisjointBoxLayout& levelDomain0 = curLevel->getStateNew().disjointBoxLayout();  
  Box Sbox = m_kineticBoxes[0];
  DataIterator dit = levelDomain0.dataIterator();    
  for (dit.begin(); dit.ok(); ++dit)  
  {        
    const Box& b = levelDomain0[dit()]; 
    // if (!Sbox.intersects(b)) continue; We do not need for level 0    
    FArrayBox& SFab = LSCData0[dit()];               
    
    SFab.copy(*m_sourceTerms[0],b,0,b,0,    SCOMP); //  a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp,  a_numcomp    
    SFab.copy(*m_neutralData[0],b,0,b,SCOMP,NCOMP); 
    
    FORT_MODIFYSOURCEAFTERKINETICSC3D(
        CHF_FRA(SFab),
        CHF_BOX(b));            

  }   
   
  // Fill data on finer levels  
  for(ilev = 1; ilev < a_levels.size(); ++ilev)
  if (a_levels[ilev]->boxes().size() > 0) 
  {                                         
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[ilev]);    
          
    AMRLevelIdealMHD* amrMHDCoarserPtr  = curLevel->getCoarserLevel(); 
    LevelData<FArrayBox>& LSCData       = curLevel->getSCData();
    LevelData<FArrayBox>& LSCDataCoarse = amrMHDCoarserPtr->getSCData();
    
    const DisjointBoxLayout& grids = LSCData.disjointBoxLayout();    
  
    // Interpolate from coarser ilev 
    if (amrMHDCoarserPtr!=NULL)
    {
      int nRefCrse = amrMHDCoarserPtr->refRatio();
      
      FineInterp fineInterpST;  
      fineInterpST.define(grids,
                          LSCData.nComp(),
                          nRefCrse,
                          curLevel->problemDomain());
                          
      fineInterpST.interpToFine(LSCData,LSCDataCoarse);
    }
    
    if (ilev>=m_ngrids) continue;
    
    // Substitute interpolated values with exact values    
                 
    const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
    Sbox = m_kineticBoxes[ilev];
    dit = levelDomain.dataIterator();    
    for (dit.begin(); dit.ok(); ++dit)  
    {        
      Box b = levelDomain[dit()]; 
      b &= Sbox;
      if (b.isEmpty()) continue;
      FArrayBox& SFab = LSCData[dit()];               
      
      SFab.copy(*m_sourceTerms[ilev],b,0,b,0    ,SCOMP); //  a_src, a_srcbox, a_srccomp, a_destbox, a_destcomp,  a_numcomp    
      SFab.copy(*m_neutralData[ilev],b,0,b,SCOMP,NCOMP); 
  
      FORT_MODIFYSOURCEAFTERKINETICSC3D(
        CHF_FRA(SFab),
        CHF_BOX(b));            
      
    }
      
    // Filling ghost cells (becuase we do not do this in main code)
    if (LSCData.ghostVect() != IntVect::Zero)
    {
      int numGhost = LSCData.ghostVect()[0];
      if (amrMHDCoarserPtr!=NULL)
      {      
        const DisjointBoxLayout& levelDomain = LSCData.disjointBoxLayout();
        
        PiecewiseLinearFillPatch pwl(levelDomain,
                                 LSCDataCoarse.disjointBoxLayout(),
                                 LSCData.nComp(),
                                 amrMHDCoarserPtr->problemDomain(),
                                 amrMHDCoarserPtr->refRatio(),
                                 numGhost);
  
        pwl.fillInterp(LSCData,
                   LSCDataCoarse,
                   LSCDataCoarse,
                   1.0,
                   0,
                   0,
                   LSCData.nComp());
      }
      
      LSCData.exchange(Interval(0,LSCData.nComp()-1));
    }

  //  char filename[60];filename[0]=0;
  //  sprintf(filename,"mhdk_st%i.dat",m_level);
  //  FILE* tfile = OpenTecplotFile(filename,"TITLE = \"source terms\"\n VARIABLES=\"X\" \"Y\" \"SMOM\" \"SENG\"");    
  //    
  //  const DisjointBoxLayout& levelDomain = m_SCData.disjointBoxLayout();  
  //  DataIterator dit = m_SCData.dataIterator();
  //  for (dit.begin(); dit.ok(); ++dit)  
  //  {            
  //    const Box& b = levelDomain[dit()]; 
  //    FArrayBox& SFab = m_SCData[dit()];                   
  //    WriteFArrayBoxToTecplotFile(tfile, SFab, SFab.box(), Interval(0,1), m_dx);            
  //  }            
  //  CloseTecplotFile(tfile);
    
  }
   
   
  gettimeofday(&time4,NULL);
  
  if (procID() == 0)
  {
    double stepTime = (double)(time4.tv_sec-time3.tv_sec) + 1e-6*((double)(time4.tv_usec-time3.tv_usec));
    char buf[20];
    sprintf(buf,"%.4g",stepTime);
    pout() << "KineticSources3D::CalculateSources time: " << buf << " sec" << endl;
  }
    
  /*if (m_lastKineticFiles.size() > 3)
  {    
    std::string fName(m_lastKineticFiles.front());
    //if (procID() == 0) remove(fName.c_str());
    m_lastKineticFiles.pop();    
  }     
  char iter_str[100];iter_str[0] = 0;
  sprintf(iter_str, "%s%06i.hdf5",  m_chk_prefix.c_str(), a_curStep );
  //writeCheckpointFile(iter_str);
  m_lastKineticFiles.push(std::string(iter_str));*/
   
}

//                                                     Add external source terms
void KineticSources3D :: addExternalSources(       FArrayBox & a_U,
                                           const FArrayBox & a_S,
                                           const FArrayBox & a_W,
                                           const Real      & a_dt,
                                           const Real      & a_dx,
                                           const Box       & a_box)
{
  FORT_ADDKINETICSOURCE3D( CHF_FRA(a_U),
                         CHF_CONST_FRA(a_S),
                         CHF_CONST_REAL(a_dt),
                         CHF_BOX(a_box) );
}

void KineticSources3D::writeCheckpointFile(const char* chkFileName) const
{
  //return;
  FORT_MC_OUTPUT_RAW_3D(chkFileName);
}

// Number additional variables for writing to plot file  
int KineticSources3D::numPlotVars()
{
  return numSourceTerms();
}
    
 // Names of the additional variables for writing to plot file 
Vector<std::string> KineticSources3D::plotNames()
{
  Vector<std::string> retval;
  
  retval.push_back("srho");
  retval.push_back("x-smom");
  retval.push_back("y-smom");
  retval.push_back("z-smom");
  retval.push_back("seng");
  retval.push_back("nchex");
  retval.push_back("nrho");
  retval.push_back("x-nvel");
  retval.push_back("y-nvel");
  retval.push_back("z-nvel");
  retval.push_back("ntemp");  
    
  return retval;
}
