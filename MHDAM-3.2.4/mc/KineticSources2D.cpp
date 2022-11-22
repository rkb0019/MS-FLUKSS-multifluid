#include "LayoutIterator.H"
#include <PiecewiseLinearFillPatch.H> 
#include "FORT_PROTO.H"

#include "KineticSources2D.H"
#include "KineticSources2DF_F.H"
#include "AMRLevelIdealMHD.H"
#include "PatchIdealMHDF_F.H"
#include "LGintegrator.H"
#include "KS2DIntergrator.H"
#include "EosCommon.H"
#include "DebugF_F.H"
#include "TecplotIO.H"


  
extern "C"
{

// Prototype for Fortran procedure mc_init 
#define FORT_MC_INIT_2D FORTRAN_NAME( RUN_MC_INIT_2D ,run_mc_init_2d )
void FORT_MC_INIT_2D
(
    Real* const p_xmax_in, 
    Real* const p_zmin_in, 
    Real* const p_zmax_in, 
    int*  const ngrids_in,
    int*  const grids_nx_in, 
    int*  const grids_nz_in, 
    Real* const grids_boundariesx_in, 
    Real* const grids_boundariesz_in, 
    int*  const restart_in, 
    char* const loadfile, 
    Real*  total_neutrals, 
    Real* const lism_nh_in, 
    Real* const lism_vh_in, 
    Real* const lism_th_in,
    int* const verbosity_in,
    int* const photoionize_in
);

#define FORT_RUN_PASS_PLASMA_GRID FORTRAN_NAME( RUN_PASS_PLASMA_GRID ,run_pass_plasma_grid )
void FORT_RUN_PASS_PLASMA_GRID
(
  int*  const kk,
  int*  const nx, 
  int*  const ny, 
  Real* const grid, 
  int*  const region
);



// Prototype for Fortran procedure mc_neutrals
#define FORT_MC_NEUTRALS_2D FORTRAN_NAME( RUN_MC_NEUTRALS_2D, run_mc_neutrals_2d)
void FORT_MC_NEUTRALS_2D(
    Real* const run_time, 
    Real* const min_nchex
);

#define FORT_RUN_RETURN_NEUTRAL_GRID FORTRAN_NAME(RUN_RETURN_NEUTRAL_GRID, run_return_neutral_grid)
void FORT_RUN_RETURN_NEUTRAL_GRID(
  int*  const kk,
  int*  const nx, 
  int*  const ny, 
  int*  const region,
  Real* const grid);
  
#define FORT_RUN_RETURN_SOURCE_GRID FORTRAN_NAME(RUN_RETURN_SOURCE_GRID, run_return_source_grid)
void FORT_RUN_RETURN_SOURCE_GRID(
  int*  const kk,
  int*  const nx,   
  int*  const ny,   
  Real* const grid
);


#define FORT_MC_OUTPUT_RAW_2D FORTRAN_NAME( RUN_MC_OUTPUT_RAW_2D ,run_mc_output_raw_2d )
void FORT_MC_OUTPUT_RAW_2D
(    
  const char* const chkfile     
);


}



                                                                   // Costructor
KineticSources2D::KineticSources2D()
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
KineticSources2D::~KineticSources2D()
{
}
                              // Factory method - this object is its own factory
SourceCalculator * KineticSources2D :: new_SourceCalculator( void )
{
  SourceCalculator* retval = new KineticSources2D();

  return retval;
}

                                                             // Input parameters
void KineticSources2D :: input( ParmParse & parser, int verbosity )
{     
  int i; 
  
  parser.query( "lismN",   m_lismN );
  parser.query( "lismV",   m_lismV );
  parser.query( "lismT",   m_lismT );
  parser.query( "netN",    m_netN  );  
  parser.query( "XC",      m_sunXC );
  parser.query( "YC",      m_sunYC );
  parser.query( "ZC",      m_sunZC );
  
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
  //  MayDay::Error("KineticSources2D::input, num_cells[0] must be equal 2*num_cells[1]");
  
  Real dx = domainLength/numCells[0];
  
  m_sunIC = (int)floor(m_sunXC/dx+0.5);
  m_sunJC = (int)floor(m_sunYC/dx+0.5);
  m_sunKC = (int)floor(m_sunZC/dx+0.5);    
  m_sunXC = m_sunIC*dx;
  m_sunYC = m_sunJC*dx;
  m_sunZC = m_sunKC*dx;
  
  
  char buf[20];std::vector<Real> tmpVect;    
  parser.query( "kinetic_ngrids",      m_ngrids );
  if (maxLevel + 1 < m_ngrids)
    MayDay::Error("KineticSources2D::input, number of additional levels should be greater than number of kinetic levels");
  
  m_ngrids++;// Total number of grids, base level added
  
  m_gridsBoundaries.resize(m_ngrids);  
  for (i=1;i<m_ngrids;i++)
  {
    sprintf(buf, "kinetic_grid%ibb", i);
    parser.queryarr( buf, tmpVect, 0, CH_SPACEDIM);    
    RealVect rv(D_DECL(tmpVect[0], tmpVect[1], tmpVect[2]));
    m_gridsBoundaries[i] = rv;    
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
void KineticSources2D :: initExternalSC(Vector<AMRLevel*> a_levels)
{  
  int i,ilev; AMRLevelIdealMHD* curLevel;  
  
  CH_assert(m_ngrids <= a_levels.size());
  
  IntVect SunPos(D_DECL(m_sunIC,m_sunJC,m_sunKC));
  
  Real dx;
  
  if( m_verbosity >= 3 )
  {
    pout() << "    KineticSources2D::initExternalSC "  << endl;    
  }
      
  m_gridsDim.resize(m_ngrids);
  ilev = 0;
  curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[ilev]);    
  dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);
  m_gridsDim[ilev] = curLevel->problemDomain().domainBox().size()-SunPos;
  m_gridsBoundaries[ilev][0] = (m_gridsDim[ilev][0])*dx;
  m_gridsBoundaries[ilev][1] = (m_gridsDim[ilev][1])*dx;
  
  // Now this is Jacob's requirement
  //assert(curLevel->problemDomain().domainBox().size(0) == 2*curLevel->problemDomain().domainBox().size(1));
  
  for (i=1;i<m_ngrids;i++)
  {
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);  
    dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,i);
    IntVect gridDim (D_DECL (    
      (int)floor(m_gridsBoundaries[i][0]/dx+0.5),
      (int)floor(m_gridsBoundaries[i][1]/dx+0.5),
      (int)floor(m_gridsBoundaries[i][2]/dx+0.5) ));      
    
    // Boundaries of finer and coarse levels must have common edge/plane 
    gridDim.coarsen(a_levels[i-1]->refRatio());
    gridDim.scale  (a_levels[i-1]->refRatio());
  
    m_gridsDim[i]   = gridDim;    
    D_DECL(
      m_gridsBoundaries[i][0] = gridDim[0]*dx,
      m_gridsBoundaries[i][1] = gridDim[1]*dx,
      m_gridsBoundaries[i][2] = gridDim[2]*dx);          
  }
  
  Real AU2Meters = 0.01*eos_AU;
  
  // Upper case letters is a notation of Chombo coordinate system
  // Lower case letters is a notation of Jacob's coordinate system
  Real* grids_boundariesX = new Real[m_ngrids+1];
  Real* grids_boundariesY = new Real[m_ngrids+1];
  int*  grids_nX          = new int [m_ngrids];
  int*  grids_nY          = new int [m_ngrids];
  for (i=0;i<m_ngrids;i++)
  {
    grids_boundariesX[i+1] = m_gridsBoundaries[m_ngrids-i-1][0]*AU2Meters;
    grids_boundariesY[i+1] = m_gridsBoundaries[m_ngrids-i-1][1]*AU2Meters;
    grids_nX[i] = m_gridsDim[m_ngrids-i-1][0];
    grids_nY[i] = m_gridsDim[m_ngrids-i-1][1];
  }  
  grids_boundariesX[0] = 0.0;
  grids_boundariesY[0] = 0.0;
  Real p_Xmax = grids_boundariesX[m_ngrids];
  Real p_Ymax = grids_boundariesY[m_ngrids];
  Real p_Xmin = -m_sunXC*AU2Meters;
  
  if (p_Xmax < fabs(p_Xmin))
    MayDay::Error("KineticSources2D::initExternalSC: p_Xmax must be greater than |p_Xmin|");
        
      
  Real  total_neutrals = m_totalNeutrals;
  
  int  restart = (m_restartFile.size() > 0 ? 1 : 0);
  char loadfile[60]={0}; loadfile[0] = 0;
  if (restart == 1) m_restartFile.copy(loadfile,m_restartFile.size());
  
  int verbosity = m_verbosity;
  verbosity = 3;  
    
  FORT_MC_INIT_2D(&p_Ymax, &p_Xmin, &p_Xmax,
                  &m_ngrids,  grids_nY, grids_nX, grids_boundariesY, grids_boundariesX, 
                  &restart,  loadfile, &total_neutrals, 
                  &m_netN, &m_lismV,  &m_lismT, &verbosity, &m_photoionize);
      
  
  delete[] grids_boundariesX;
  delete[] grids_boundariesY;
  delete[] grids_nX;
  delete[] grids_nY;
  
  
  m_gridsBoxes.resize(m_ngrids);  
  
  curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[0]);
  m_gridsBoxes[0] = curLevel->problemDomain().domainBox();
  SunPos.scale(curLevel->refRatio());
  
  for (i=1;i<m_ngrids;i++)
  {
    IntVect LCorner(SunPos),UCorner(SunPos);
    
    LCorner[0]-=m_gridsDim[i][0];
    UCorner[0]+=m_gridsDim[i][0]-1;        
    UCorner[1]+=m_gridsDim[i][1]-1;        
            
    Box& bToCopy = m_gridsBoxes[i];    
    bToCopy.define(LCorner,UCorner);
    
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);    
    //if (i==0) CH_assert(bToCopy == curLevel->problemDomain().domainBox());
        
    // Indices of sun position on next level
    SunPos.scale(curLevel->refRatio());
  }
}


//                                       Check time for source terms calculation
bool KineticSources2D :: checkTime( Real dTime, int iStep )
{
  //return true;
  
  if (iStep == m_cur_step)
    MayDay::Error("KineticSources2D::checkTime must never be called twice for the same step");
  
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

// number of variables that are needed for external source calculator
int KineticSources2D::numVarsForSC()
{
  return 4;
}

//                                             Number of calculated source terms
//                                                    and neutral parameters
int KineticSources2D :: numSourceTerms( void )
{
  return 21;
}



// Using conservative variables a_U, the method prepares data (stored in a_Buf)
// that will be used in an external source calculator.
void KineticSources2D::PrepareDataForSC(const FArrayBox& a_U, FArrayBox& a_Buf) 
{
  const Box& b = a_Buf.box();
  FORT_PREPAREDATAFORKINETICSC2D(
      CHF_CONST_FRA(a_U),
      CHF_FRA(a_Buf),
      CHF_BOX(b));
}

// Send data to zero processor
void KineticSources2D::sendDataToZeroProc(Vector<AMRLevel*> a_levels)
{
#ifdef CH_MPI  
  if( m_verbosity >= 3 )
  {
    pout() << "KineticSources2D::sendDataToZeroProc" << endl;
  }
  
  int igrid;

  int curProc = procID();  
  if (curProc==0) return; // Avoid sending data to itself
  
  // Used to debug broadcasting
  char file_name[50];  
  static bool isWritten = true;
  FILE* tfile = NULL;      
      
  for (igrid = 0; igrid < m_ngrids; igrid++)
  {        
    Box& bToCopy = m_gridsBoxes[igrid];
    AMRLevelIdealMHD* curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[igrid]);
    Real dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,igrid);
    
    if (igrid==0) CH_assert(bToCopy == curLevel->problemDomain().domainBox());
    
    const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
    DataIterator dit = levelDomain.dataIterator();
    
    int BoxDataSize,numBoxes;
    Real* BoxData;
          
    int RealBufferSize;
    Real *RealBuffer;
    
    int nSendVariables = numVarsForSC();                
             
    numBoxes = 0;RealBufferSize=0;
    // Iterator of all grids on this level  
    for (dit.begin(); dit.ok(); ++dit)  
    {    
      const Box& b = levelDomain[dit()];    
      if (!bToCopy.intersects(b)) continue;
      
      numBoxes++;   
         
      BoxDataSize = b.numPts()*nSendVariables;
      if (CH_SPACEDIM == 2)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*nSendVariables));
      if (CH_SPACEDIM == 3)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*b.size(2)*nSendVariables));
      RealBufferSize+=BoxDataSize;            
    }
    
    if (numBoxes==0) 
    {
      Real realBuf1 = 0.0;
      MPI_Send(&realBuf1,1 /*buffer size*/, MPI_CH_REAL, 0 /* processor */, igrid /* tag */, Chombo_MPI::comm);    
      continue;
    }
        
    if (!isWritten)
    {
      sprintf(file_name,"cpu%ilev%i.dat",curProc,igrid);
      tfile = OpenTecplotFile(file_name,"TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"rho\" \"u\" \"v\" \"T\" ");
    }    
        
    RealBuffer=new Real[RealBufferSize];    
    int posInRealBuffer = 0;  
        
    for (dit.begin(); dit.ok(); ++dit)  
    {        
      const Box& b = levelDomain[dit()];
      const FArrayBox& UFab = curLevel->getStateNew()[dit()];	
      if (!bToCopy.intersects(b)) continue;
                
      FArrayBox PhisVars(b,nSendVariables);      
      PrepareDataForSC(UFab,PhisVars);
      
      if (!isWritten) WriteFArrayBoxToTecplotFile(tfile, PhisVars, b, Interval(0, PhisVars.nComp()-1), dx);
      
      #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(PhisVars)
        );
      #endif
          
      BoxDataSize = b.numPts()*nSendVariables;
      if (CH_SPACEDIM == 2) CH_assert(BoxDataSize == (b.size(0)*b.size(1)*nSendVariables));
      if (CH_SPACEDIM == 3) CH_assert(BoxDataSize == (b.size(0)*b.size(1)*b.size(2)*nSendVariables));
      BoxData = PhisVars.dataPtr();
      
      memcpy(RealBuffer+posInRealBuffer,BoxData,sizeof(Real)*BoxDataSize);
      posInRealBuffer+=BoxDataSize;
    }
      
    CH_assert(posInRealBuffer==RealBufferSize);
      
    MPI_Send(RealBuffer,RealBufferSize,MPI_CH_REAL,0 /* processor */ , igrid /* tag */ ,Chombo_MPI::comm);
      
    delete[] RealBuffer;        
  }
  
  if (!isWritten) CloseTecplotFile(tfile);
  isWritten = true;
#endif
}

// Gather data on zero processor. Gathered data are stored in a_Data
void KineticSources2D::GatherDataOnZeroProc(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_Data)
{
  if( m_verbosity >= 3 )
  {
    pout() << "KineticSources2D::GatherDataOnZeroProcessor" << endl;
  }
  
  int curProc = procID();  
  if (curProc!=0) return;  
  
  int nGetVariables = numVarsForSC();
  AMRLevelIdealMHD* curLevel;
  
  int i,igrid;
    
  // Gather data on "0" processor  
#ifdef CH_MPI      
  for (igrid = 0; igrid < m_ngrids; igrid++)
  {       
    Box& bToCopy = m_gridsBoxes[igrid];
    FArrayBox& FBox = (*a_Data[igrid]);    
    
    // On zero level we have the grid that is larger than problem domain.
    // Here we initialize a_Data on this level
    if (igrid == 0)
    {
      FBox.setVal(m_lismN*1e+6, KRHO);
      FBox.setVal(m_lismV*1e-2, KVELX);
      FBox.setVal(0.0         , KVELY);
      FBox.setVal(m_lismT     , KTEMP);
    }
    
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[igrid]);
    const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
    

    int iCount,iBox, posInRealBuffer;
    int **intBuffers=NULL, *nReceivedMessages=NULL, *sizesOfRealBuffers=NULL;
    
    nReceivedMessages  = (numProc()>1 ? new int[numProc()-1] : NULL);
    sizesOfRealBuffers = (numProc()>1 ? new int[numProc()-1] : NULL);
    for (i=0;i<numProc()-1;i++)
    {
      nReceivedMessages[i]=0;
    }
  
    MPI_Status ProbeStatus, RecvStatus;  
    char error_buffer[1000];int actual_length,rescode;
    
    while (1)
    {
      bool stopwhile = true;
      for (i=0;i<numProc()-1;i++) stopwhile = stopwhile && (nReceivedMessages[i]==1);
      if (stopwhile) break;
    
      rescode = MPI_Probe(MPI_ANY_SOURCE, igrid /* tag */, Chombo_MPI::comm, &ProbeStatus);    
      if (rescode!=MPI_SUCCESS)
      {
        MPI_Error_string(rescode,error_buffer,&actual_length);
      }    
      CH_assert(ProbeStatus.MPI_TAG==igrid);
                      
      // Get data 
      CH_assert((nReceivedMessages[ProbeStatus.MPI_SOURCE-1]==0));
      rescode = MPI_Get_count(&ProbeStatus,MPI_CH_REAL,&iCount);
      sizesOfRealBuffers[ProbeStatus.MPI_SOURCE-1]=iCount; 
      Real* RealBuffer=new Real[iCount];
      rescode = MPI_Recv(RealBuffer,
        iCount,
        MPI_CH_REAL,
        ProbeStatus.MPI_SOURCE,
        ProbeStatus.MPI_TAG,
        Chombo_MPI::comm,&RecvStatus);        
        
      nReceivedMessages[ProbeStatus.MPI_SOURCE-1]++;
      CH_assert((ProbeStatus.MPI_SOURCE==RecvStatus.MPI_SOURCE)&&(ProbeStatus.MPI_TAG==RecvStatus.MPI_TAG));
              
      if (iCount == 1)  // no actual data          
      {
        delete[] RealBuffer;      
        continue; 
      }
                              
      posInRealBuffer = 0;
                   
      LayoutIterator lit = levelDomain.layoutIterator();
      for (lit.begin();lit.ok();++lit)        
      if (levelDomain.procID(lit.i()) == RecvStatus.MPI_SOURCE)
      {            
        const Box& b = levelDomain[lit()];                
        if (!bToCopy.intersects(b)) continue;
                   
        Real* BoxData = &RealBuffer[posInRealBuffer];
        FArrayBox fab(b,nGetVariables,BoxData);
        
        #ifndef NDEBUG
        FORT_VIEWBOXDATA(
          CHF_FRA(fab)
          );
        #endif
        
        FBox.copy(fab);
        
        posInRealBuffer+=b.numPts()*nGetVariables;      
      }
      CH_assert(posInRealBuffer==iCount);
      
      delete[] RealBuffer;      
    }    
    if (nReceivedMessages != NULL) delete[] nReceivedMessages;
    if (sizesOfRealBuffers != NULL) delete[] sizesOfRealBuffers;        
  }
#endif
    
  
  // Handle self data       
  for (i=0;i<m_ngrids;i++)
  {    
    Box& bToCopy = m_gridsBoxes[i];
    FArrayBox& FBox = (*a_Data[i]);
    if (FBox.nComp() == 0) FBox.define(bToCopy,nGetVariables);
    
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);
    const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
      
    DataIterator dit = levelDomain.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)  
    {        
      const Box& b = levelDomain[dit()];
      const FArrayBox& UFab = curLevel->getStateNew()[dit()];	            
      if (!bToCopy.intersects(b)) continue;
             
      FArrayBox PhisVars(b,nGetVariables);      
      PrepareDataForSC(UFab,PhisVars);        
      FBox.copy(PhisVars);      
    }
    
  }
  
  static bool isWritten = true;
  if (!isWritten)
  for (i=0;i<m_ngrids;i++)
  {
    char file_name[50];  
    FArrayBox& FBox = (*a_Data[i]);
    AMRLevelIdealMHD* curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);
    Real dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,i);
    sprintf(file_name,"lev%i.dat",i);
    FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"rho\" \"u\" \"v\" \"T\" ");    
    WriteFArrayBoxToTecplotFile(tfile, FBox, FBox.box(), Interval(0,FBox.nComp()-1), dx);
    CloseTecplotFile(tfile);
  }
  isWritten = true;

}


void KineticSources2D::CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep)
{  
  int i,nx,ny,igrid,reg;
  Real dx;  
  AMRLevelIdealMHD* curLevel;
  
  int curProc = procID();  
  
  Real cur_time = a_time;
  int  cur_step = a_curStep;    
  
  if ((abs(m_verbosity)>=1) && (procID() == 0))
  {
    pout() << "KineticSources2D::CalculateSources: step/time for new source terms" << endl;
  }
      
  sendDataToZeroProc(a_levels);
  
  // Boxes that define Kinetic grid. The same as m_gridsBoxes for all levels except base level
  Vector <Box> KineticBoxes(m_gridsBoxes);
  
  IntVect SunPos(D_DECL(m_sunIC,m_sunJC,m_sunKC));
  Box Box0    = m_gridsBoxes[0];  
  IntVect sE0 = Box0.smallEnd();
  IntVect bE0 = Box0.bigEnd();
  // Check if the heliotail has bigger size 
  CH_assert(m_sunIC < bE0[0] - m_sunIC);
  CH_assert(sE0==IntVect::Zero);
  bE0-=SunPos;    
  sE0-=(bE0+1);
  if (CH_SPACEDIM == 2) sE0[1]=0;
  KineticBoxes[0].define(sE0,bE0);  
  KineticBoxes[0]+=SunPos;
  CH_assert(KineticBoxes[0].bigEnd()==m_gridsBoxes[0].bigEnd());
    
  Vector<FArrayBox*> DataForSC;      
  Vector<FArrayBox*> SourceTerms;    
  Vector< BaseFab<int>* > Region;  
  DataForSC.resize(m_ngrids);
  SourceTerms.resize(m_ngrids);
  Region.resize(m_ngrids);
  for (i=0;i<m_ngrids;i++)
  {        
    Box& bToCopy = KineticBoxes[i];    
    DataForSC[i] = new FArrayBox(bToCopy,numVarsForSC());
    SourceTerms[i] = new FArrayBox(bToCopy,numSourceTerms());    
    Region[i] = new BaseFab<int>(bToCopy,1);    
  }  
  GatherDataOnZeroProc(a_levels,DataForSC);

#ifdef CH_MPI        
  int BufferSize;
  for (i=0;i<m_ngrids;i++)
  {
    BufferSize = DataForSC[i]->box().numPts()*DataForSC[i]->nComp();
    MPI_Bcast(DataForSC[i]->dataPtr(0),BufferSize,MPI_CH_REAL,0 /* rank of broadcast root */,Chombo_MPI::comm);
  }
#endif  

  // Base level is the special case
  FArrayBox Region0, DataForSC0;

  for (i=0;i<m_ngrids;i++)
  {        
    Box& bToCopy = KineticBoxes[i];            
    curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);
    dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,i);
    FORT_DEFINE_REGIONS_KINETIC( CHF_CONST_FRA((*DataForSC[i])),
                          CHF_FIA1((*Region[i]),0),
                          CHF_CONST_REAL(dx),
                          CHF_BOX(bToCopy) );        
  }  
      
  igrid = 1;
  for (i=m_ngrids-1;i>=0;i--)
  {
    nx = m_gridsDim[i][0];
    ny = m_gridsDim[i][1];
    CH_assert(DataForSC[i]->box().size(0)   == 2*nx);
    CH_assert(DataForSC[i]->box().size(1)   == ny);    
    FORT_RUN_PASS_PLASMA_GRID( &igrid,  &nx, &ny, DataForSC[i]->dataPtr(0),Region[i]->dataPtr(0));  
    igrid++;
  }
  
  
  Real runtimeYears = m_runtimeYears;    
  int  restart = (m_restartFile.size() > 0 ? 1 : 0);
  if ((restart == 0)&&(m_firstCall)) runtimeYears = m_initialRuntimeYears;
  
  Real min_nchex;
  if ((abs(m_verbosity)>=2) && (procID() == 0))
  {
    pout() << "FORT_MC_NEUTRALS_2D" << endl;
  }
  
  FORT_MC_NEUTRALS_2D(&runtimeYears, &min_nchex);
  
  m_firstCall = false;
  
  
  if (curProc == 0)
  {
    // Get source terms
    igrid = 1;    
    for (i=m_ngrids-1;i>=0;i--)
    {
      nx = m_gridsDim[i][0];
      ny = m_gridsDim[i][1];
      Box& bToCopy = KineticBoxes[i];    
      FArrayBox FAB(bToCopy,NCHEX+1);
      
      FORT_RUN_RETURN_SOURCE_GRID( &igrid,  &nx, &ny, FAB.dataPtr(0));        
      SourceTerms[i]->copy(FAB,0,SRHO,FAB.nComp());
      
      // For debugging only
      //SourceTerms[i]->setVal(0.0);
      //SourceTerms[i]->copy((*DataForSC[i]), 0, 0, DataForSC[i]->nComp());
      
      igrid++;
    }
      
    // Get neutral data      
    for (reg =-1;reg<4;reg++)  
    {
      if (reg == 0) continue; // We don't have region 0 now, but Jacob has.
      igrid = 1;
      int destcomp;
      for (i=m_ngrids-1;i>=0;i--)
      {
        nx = m_gridsDim[i][0];
        ny = m_gridsDim[i][1];
        Box& bToCopy = KineticBoxes[i];    
        FArrayBox FAB(bToCopy,4);
        destcomp = NRHO;
        if (reg>0) destcomp = NRHO + 4*reg;
        FORT_RUN_RETURN_NEUTRAL_GRID( &igrid,  &nx, &ny, &reg, FAB.dataPtr(0));  
        SourceTerms[i]->copy(FAB,0,destcomp,FAB.nComp());
        igrid++;
      }
    }    
    
    static bool isWritten = false;
    if (!isWritten)
    for (i=0;i<m_ngrids;i++)
    {
      char file_name[50];        
      AMRLevelIdealMHD* curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);
      dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,i);
      sprintf(file_name,"st%iorig.dat",i);
      FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"SRHO\" \"SMOMX\" \"SMOMY\" \"SENG\" \"NRHO\" \"NVELX\" \"NVELY\"  \"NTEMP\" ");      
      WriteFArrayBoxToTecplotFile(tfile, (*SourceTerms[i]), SourceTerms[i]->box(), Interval(0,NTEMP), dx);
      CloseTecplotFile(tfile);
    }
              
    for (i=0;i<m_ngrids;i++)
    {            
      ModifySourceTerms((*SourceTerms[i]),SourceTerms[i]->box());    
    }  
    
    if (!isWritten)
    for (i=0;i<m_ngrids;i++)
    {
      char file_name[50];        
      AMRLevelIdealMHD* curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[i]);
      dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,i);
      sprintf(file_name,"st%itrans.dat",i);
      FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"Data for SC\"\n VARIABLES=\"X\" \"Y\" \"SRHO\" \"SMOMX\" \"SMOMY\" \"SENG\" \"NRHO\" \"NVELX\" \"NVELY\"  \"NTEMP\" ");      
      WriteFArrayBoxToTecplotFile(tfile, (*SourceTerms[i]), SourceTerms[i]->box(), Interval(0,NTEMP), dx);
      CloseTecplotFile(tfile);
    }
  }   
   
  BroadcastSourceTerms(a_levels, SourceTerms);   
      
  
  PrepareSourceTermsOnLevels(a_levels,SourceTerms,Region);
  
  for (i=0;i<m_ngrids;i++)
  {            
    delete DataForSC[i];
    delete SourceTerms[i];    
    delete Region[i];
  } 
    
  if (m_lastKineticFiles.size() > 3)
  {    
    std::string fName(m_lastKineticFiles.front());
    //if (procID() == 0) remove(fName.c_str());
    m_lastKineticFiles.pop();    
  }     
  char iter_str[100];iter_str[0] = 0;
  sprintf(iter_str, "%s%06i.hdf5",  m_chk_prefix.c_str(), a_curStep );
  //writeCheckpointFile(iter_str);
  m_lastKineticFiles.push(std::string(iter_str));
   
}

void KineticSources2D::BroadcastSourceTerms(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_Data)
{
#ifdef CH_MPI      
  if( m_verbosity >= 3 )
  {
    pout() << "KineticSources2D::BroadcastSourceTerms" << endl;
  }
  int i;
  
  int curProc = procID();          
   
  int BoxDataSize;             
  int BufferSize = 0;      
  int nSendVariables = numSourceTerms();
             
  for (i=0;i<m_ngrids;i++)
  {    
    const Box& b = m_gridsBoxes[i];
    BoxDataSize = b.numPts()*nSendVariables;
    if (CH_SPACEDIM == 2)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*nSendVariables));
    if (CH_SPACEDIM == 3)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*b.size(2)*nSendVariables));
    BufferSize += BoxDataSize;            
  }
  
  FArrayBox zeroLevelData(m_gridsBoxes[0],numSourceTerms());
  zeroLevelData.copy(*a_Data[0]);      
    
  Real *BoxData;    
  Real *RealBuffer=new Real[BufferSize];    
  int posInRealBuffer = 0;  
  
  // Prepare data on master for broadcasting
  if (curProc == 0)
  {
    for (i=0;i<m_ngrids;i++)
    {        
      const Box& b = m_gridsBoxes[i];
            
      FArrayBox& FBox = (i==0 ? zeroLevelData : (*a_Data[i]));	    
                          
      BoxDataSize = b.numPts()*nSendVariables;
      if (CH_SPACEDIM == 2) CH_assert(BoxDataSize == (b.size(0)*b.size(1)*nSendVariables));
      if (CH_SPACEDIM == 3) CH_assert(BoxDataSize == (b.size(0)*b.size(1)*b.size(2)*nSendVariables));
      BoxData = FBox.dataPtr();
      
      memcpy(RealBuffer+posInRealBuffer,BoxData,sizeof(Real)*BoxDataSize);
      posInRealBuffer+=BoxDataSize;
    }      
    CH_assert(posInRealBuffer==BufferSize);
  }
  
  MPI_Bcast(RealBuffer,BufferSize,MPI_CH_REAL,0 /* rank of broadcast root */,Chombo_MPI::comm);
  
  if (curProc != 0)
  {
    posInRealBuffer=0; 
    for (i=0;i<m_ngrids;i++)
    {        
      const Box& b    = m_gridsBoxes[i];
      FArrayBox& FBox = (*a_Data[i]);	    
      
      BoxData = &RealBuffer[posInRealBuffer];
      FArrayBox TmpFBox(b, numSourceTerms(), BoxData);	    
      FBox.copy(TmpFBox);      
                          
      BoxDataSize = b.numPts()*nSendVariables;
      posInRealBuffer+=BoxDataSize;                            
    }
    CH_assert(posInRealBuffer==BufferSize);
  }        
                
  delete[] RealBuffer;                
#endif                
}


// Modification of source terms calculated by external source calculator.
void KineticSources2D::ModifySourceTerms(FArrayBox& a_SourceTerms, const Box& a_box)
{  
  FORT_MODIFYSOURCEAFTERKINETICSC2D(
        CHF_FRA(a_SourceTerms),
        CHF_BOX(a_box));            

}

// Interpolate source terms from zero level to higher levels.
void KineticSources2D::PrepareSourceTermsOnLevels(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_SCData, Vector< BaseFab<int>* > Region) 
{
  int level,i,j,ivar;
  Vector<AMRLevel*> amrlevels = a_levels;
  
  AMRLevelIdealMHD *amrMHDCoarserPtr,*curLevel;
  
  // Copy to data to Level 0
  curLevel = static_cast<AMRLevelIdealMHD*>(a_levels[0]);  
  LevelData<FArrayBox>& LSCData0 = curLevel->getSCData(); 
  FArrayBox* SCData = a_SCData[0];
  const DisjointBoxLayout& levelDomain0 = curLevel->getStateNew().disjointBoxLayout();  
  Box Sbox = m_gridsBoxes[0];
  DataIterator dit = levelDomain0.dataIterator();    
  for (dit.begin(); dit.ok(); ++dit)  
  {        
    const Box& b = levelDomain0[dit()]; 
    if (!Sbox.intersects(b)) continue;
    FArrayBox& SFab = LSCData0[dit()];               
    SFab.copy(*SCData);                
  }
  
  int nAverageCells = 1;
  
  // Fill data on finer levels  
  for(level = 1; level < amrlevels.size(); ++level)
  if (amrlevels[level]->boxes().size() > 0) 
  {                                         
    curLevel = static_cast<AMRLevelIdealMHD*>(amrlevels[level]);    
          
    amrMHDCoarserPtr = curLevel->getCoarserLevel(); 
    LevelData<FArrayBox>& LSCData = curLevel->getSCData();
    LevelData<FArrayBox>& LSCDataCoarse = amrMHDCoarserPtr->getSCData();
    
    const DisjointBoxLayout& grids = LSCData.disjointBoxLayout();
    
    nAverageCells = amrlevels[level-1]->refRatio();
  
    // Interpolate from coarser level 
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
    
    // Substitute interpolated values with exact values
    if (level < a_SCData.size())
    {
      SCData = a_SCData[level];
      
      // Average data near the axis
      IntVect LCorner=SCData->box().smallEnd();
      IntVect UCorner=SCData->box().bigEnd();
      IntVect aux_IntVect;int region;
      if (level == amrlevels.size()-1)
      for (i=LCorner[0];i<=UCorner[0];i++)
      for (ivar = 0; ivar < SCData->nComp(); ivar++)
      if ((ivar == SMOMX) || (ivar == SMOMY))       
      {
        /*D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=1;,aux_IntVect[2]=0;);
        Region[level]->getVal(&region,aux_IntVect);
        if (region == 2)
        {
          Real slope    = SCData->get(IntVect(D_DECL(i,2,0)),ivar) - SCData->get(IntVect(D_DECL(i,1,0)),ivar);          
          Real newValue = SCData->get(IntVect(D_DECL(i,1,0)),ivar) - slope;
          SCData->set(IntVect(D_DECL(i,0,0)),ivar,newValue);
        }*/
      
        // Simple averaging
        
        //Real Avg = 0.0;
        //for (j=0; j<nAverageCells; j++)
        //{
        //  D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=0;);
        //  Avg += SCData->get(aux_IntVect,ivar);
        //}
        //Avg /= nAverageCells;
        //for (j=0; j<nAverageCells; j++)
        //{
        //  D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=0;);
        //  SCData->set(aux_IntVect,ivar,Avg);
        //}
        
        // Filling first value with averaging of 2nd and 3rd
        //for (j=1; j<3; j++)
        //{
        //  D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=0;);
        //  Avg += SCData->get(aux_IntVect,ivar);
        //}
        //Avg /= 2.0;
        //for (j=0; j<3; j++)
        //{
        //  D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=0;);
        //  SCData->set(aux_IntVect,ivar,Avg);
        //}
      }
      
      const DisjointBoxLayout& levelDomain = curLevel->getStateNew().disjointBoxLayout();  
      Sbox = m_gridsBoxes[level];
      dit = levelDomain.dataIterator();    
      for (dit.begin(); dit.ok(); ++dit)  
      {        
        const Box& b = levelDomain[dit()]; 
        if (!Sbox.intersects(b)) continue;
        FArrayBox& SFab = LSCData[dit()];               
        SFab.copy(*SCData);                
      }
    }
  
    // Filling ghost cells
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
}

//                                                     Add external source terms
void KineticSources2D :: addExternalSources(       FArrayBox & a_U,
                                           const FArrayBox & a_S,
                                           const FArrayBox & a_W,
                                           const Real      & a_dt,
                                           const Real      & a_dx,
                                           const Box       & a_box)
{
  FORT_ADDKINETICSOURCE2D( CHF_FRA(a_U),
                         CHF_CONST_FRA(a_S),
                         CHF_CONST_REAL(a_dt),
                         CHF_BOX(a_box) );
}

void KineticSources2D::writeCheckpointFile(const char* chkFileName) const
{
  //return;
  FORT_MC_OUTPUT_RAW_2D(chkFileName);
}

// Number additional variables for writing to plot file  
int KineticSources2D::numPlotVars()
{
  return numSourceTerms();
}
    
 // Names of the additional variables for writing to plot file 
Vector<std::string> KineticSources2D::plotNames()
{
  Vector<std::string> retval;
  
  retval.push_back("srho");
  retval.push_back("x-smom");
  retval.push_back("y-smom");
  retval.push_back("seng");
  retval.push_back("nchex");
  retval.push_back("nrho");
  retval.push_back("x-nvel");
  retval.push_back("y-nvel");
  retval.push_back("ntemp");
  
  retval.push_back("nrho1");
  retval.push_back("x-nvel1");
  retval.push_back("y-nvel1");
  retval.push_back("ntemp1");
  
  retval.push_back("nrho2");
  retval.push_back("x-nvel2");
  retval.push_back("y-nvel2");
  retval.push_back("ntemp2");
  
  retval.push_back("nrho3");
  retval.push_back("x-nvel3");
  retval.push_back("y-nvel3");
  retval.push_back("ntemp3");
    
  return retval;
}