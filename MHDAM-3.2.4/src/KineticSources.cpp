#include "LayoutIterator.H"
#include <PiecewiseLinearFillPatch.H> 

#include "KineticSources.H"
#include "AMRLevelIdealMHD.H"
#include "PatchIdealMHDF_F.H"
#include "LGintegrator.H"
#include "DebugF_F.H"
#include "EosCommon.H"
#include "tecplotF_F.H"




                                                                   // Costructor
KineticSources::KineticSources()
{
}
                                                                   // Destructor
KineticSources::~KineticSources()
{
}
                              // Factory method - this object is its own factory
                              //  Dummy method. No sense for abstract classes.
SourceCalculator * KineticSources :: new_SourceCalculator( void )
{
  CH_assert(0);
//  SourceCalculator* retval = new KineticSources();
//  return retval;
  return NULL;
}
                                                             // Input parameters


// Send all data to zero processor 
void KineticSources::sendDataToZeroProcessor(const AMRLevelIdealMHD* aLevel0)
{
#ifdef CH_MPI  
  int curProc = procID();  
  if (curProc==0) return; // Avoid sending data to itself

    
  const AMRLevelIdealMHD* Level0=aLevel0;

  //PatchMHDAM* curPatch = Level0->getpatchMHDAM();
  //assert(curPatch!=NULL);
  
  //PhysProblem* pProblem = curPatch->getPhysProblem();  
  //if (pProblem->physicalModel()!=PhysProblem::PP_2FluidPM) return;
  
  const DisjointBoxLayout& levelDomain = Level0->getStateNew().disjointBoxLayout();  
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
    numBoxes++;   
    const Box& b = levelDomain[dit()];    
       
    BoxDataSize = b.numPts()*nSendVariables;
    if (CH_SPACEDIM == 2)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*nSendVariables));
    if (CH_SPACEDIM == 3)   CH_assert(BoxDataSize == (b.size(0)*b.size(1)*b.size(2)*nSendVariables));
    RealBufferSize+=BoxDataSize;            
  }
  
  if (numBoxes==0) 
  {
    Real realBuf1 = 0.0;
    MPI_Send(&realBuf1,1 /*buffer size*/, MPI_CH_REAL, 0 /* processor */, 2 /* tag */, Chombo_MPI::comm);    
    return;
  }
  
  
  RealBuffer=new Real[RealBufferSize];
  
  int posInRealBuffer = 0;  
      
  for (dit.begin(); dit.ok(); ++dit)  
  {        
    const Box& b = levelDomain[dit()];
    const FArrayBox& UFab = Level0->getStateNew()[dit()];	
        
    
    FArrayBox PhisVars(b,nSendVariables);
    
    PrepareDataForSC(UFab,PhisVars);
    
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
    
  MPI_Send(RealBuffer,RealBufferSize,MPI_CH_REAL,0 /* processor */ ,2 /* tag */ ,Chombo_MPI::comm);
    
  delete[] RealBuffer;
#endif
}

void KineticSources::GatherDataOnZeroProcessor(const AMRLevelIdealMHD* aLevel0, FArrayBox& Level0Data)
{
  if( aLevel0->verbosity() >= 2 )
  {
    pout() << "KineticSources::GatherDataOnZeroProcessor" << endl;
  }
  
  int curProc = procID();  
  if (curProc!=0) return;  
  
  int nGetVariables = numVarsForSC();
  
  const AMRLevelIdealMHD* Level0=aLevel0;
 
  
  const DisjointBoxLayout& levelDomain = Level0->getStateNew().disjointBoxLayout();  
    
  // Gather data on "0" processor  
#ifdef CH_MPI      

  int i,iCount,iBox, posInRealBuffer;
  int **intBuffers=NULL, *nReceivedMessages=NULL, *sizesOfRealBuffers=NULL;
  
  nReceivedMessages = (numProc()>1 ? new int[numProc()-1] : NULL);
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
  
    rescode = MPI_Probe(MPI_ANY_SOURCE, 2 /* tag */, Chombo_MPI::comm, &ProbeStatus);    
    if (rescode!=MPI_SUCCESS)
    {
      MPI_Error_string(rescode,error_buffer,&actual_length);
    }    
    CH_assert(ProbeStatus.MPI_TAG==2);
                    
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
                 
      Real* BoxData = &RealBuffer[posInRealBuffer];
      FArrayBox fab(b,nGetVariables,BoxData);
      
      #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(fab)
        );
      #endif
      
      Level0Data.copy(fab);
      
      posInRealBuffer+=b.numPts()*nGetVariables;      
    }
    CH_assert(posInRealBuffer==iCount);
    
    delete[] RealBuffer;      
  }    
  
  if (nReceivedMessages != NULL) delete[] nReceivedMessages;
  if (sizesOfRealBuffers != NULL) delete[] sizesOfRealBuffers;
#endif
  
  // Handle self data     
  DataIterator dit = levelDomain.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)  
  {        
    const Box& b = levelDomain[dit()];
    const FArrayBox& UFab = Level0->getStateNew()[dit()];	
           
    FArrayBox PhisVars(b,nGetVariables);
    
    PrepareDataForSC(UFab,PhisVars);
      
    Level0Data.copy(PhisVars);      
  }
  
  //int NITEMS = 3, cur_step = 0;  
  //FORT_OPENTECFILE(CHF_CONST_INT(cur_step),CHF_CONST_INT(NITEMS));
  // Real dx0 = Level0->dx();
  //FORT_OUTPUTTECPLOT(   CHF_CONST_FRA(Level0Data),
  //                      CHF_CONST_REAL(dx0),
  //                      CHF_BOX(Level0Data.box()) );        
  //FORT_CLOSETECFILE(CHF_CONST_INT(cur_step));
}


void KineticSources::BroadcastSourceTerms(AMRLevelIdealMHD* aLevel0, FArrayBox& a_SourceTerms)
{
  if( aLevel0->verbosity() >= 2 )
  {
    pout() << "KineticSources::BroadcastSourceTerms" << endl;
  }

  int i;
  int iProc, BoxDataSize, iCount, posInRealBuffer;   
  
  int curProc = procID();  
  
  
  AMRLevelIdealMHD* Level0=aLevel0;    
  const DisjointBoxLayout& levelDomain = Level0->getStateNew().disjointBoxLayout();  
  
  int numST = a_SourceTerms.nComp(); // Number of source terms
  
  // Zero processor sends to others theirs parts of source terms.
  if (curProc==0)
  {
#ifdef CH_MPI      
    int* sizesOfRealBuffers = new int[numProc()];
    Real** RealBuffers=new Real*[numProc()];
    for (iProc=0;iProc<numProc();iProc++) 
    {
      sizesOfRealBuffers[iProc] = 0;
      RealBuffers[iProc] = NULL;
    }
    
    LayoutIterator lit = levelDomain.layoutIterator();
    for (lit.begin();lit.ok();++lit)        
    {
      const Box& b = levelDomain[lit()];
      BoxDataSize = b.numPts()*a_SourceTerms.nComp();
      sizesOfRealBuffers[levelDomain.procID(lit.i())] += BoxDataSize;
    }    
    
    MPI_Request* request_array = (numProc()>1 ? new MPI_Request[numProc()-1] : NULL);
    MPI_Status* status_array = (numProc()>1 ? new MPI_Status[numProc()-1] : NULL);
        
    for (iProc=1;iProc<numProc();iProc++) 
    {
      iCount = sizesOfRealBuffers[iProc]; 
      CH_assert(iCount>=0);
      
      if (iCount == 0)  // no actual data
      {        
        Real realBuf1 = 0.0;          
        MPI_Isend(&realBuf1,1 /*buffer size*/,MPI_CH_REAL,iProc,0 /* tag */,Chombo_MPI::comm, &request_array[iProc-1]);
        continue;
      }
      
      RealBuffers[iProc]=new Real[iCount];                  
      posInRealBuffer = 0;
      
      LayoutIterator lit = levelDomain.layoutIterator();
      for (lit.begin();lit.ok();++lit)        
      if (levelDomain.procID(lit.i()) == iProc)
      {
        const Box& b = levelDomain[lit()];
        
        FArrayBox fab(b,numSourceTerms());
        fab.copy(a_SourceTerms);
                        
        #ifndef NDEBUG
        FORT_VIEWBOXDATA(
          CHF_FRA(fab)
          );
        #endif          
                
        Real* BoxData = fab.dataPtr();        
        BoxDataSize = b.numPts()*fab.nComp();
        memcpy(RealBuffers[iProc]+posInRealBuffer,BoxData,sizeof(Real)*BoxDataSize);
        posInRealBuffer+=BoxDataSize;       
      }
      CH_assert(posInRealBuffer==iCount);      
      MPI_Isend(RealBuffers[iProc],iCount,MPI_CH_REAL,iProc,0 /* tag */,Chombo_MPI::comm, &request_array[iProc-1]);      
    }
#endif          
                
    LevelData<FArrayBox>& MySourceTerms = Level0->getSCData();    
    DataIterator dit = levelDomain.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)  
    {        
      const Box& b = levelDomain[dit()]; 
      FArrayBox& SFab = MySourceTerms[dit()];               
      SFab.copy(a_SourceTerms);            
      ModifySourceTerms(SFab,b);        
    }            
    
#ifdef CH_MPI          
    if (numProc()>1)
    {
      int rescode = MPI_Waitall(numProc()-1,request_array,status_array);
      if (rescode!=MPI_SUCCESS) MayDay::Error("KineticSources::BroadcastSourceTerms. MPI_Waitall communication error");
      delete[] request_array;
      delete[] status_array;
    }
        
    for (i=0;i<numProc();i++)  
    if (RealBuffers[i]!=NULL) delete[] RealBuffers[i];
    delete[] RealBuffers;  
    
    delete[] sizesOfRealBuffers;
#endif              
    
  } // if (curProc==0)
 
#ifdef CH_MPI           
  // Other processors receive source terms from zero processor.
  if (curProc>0)
  {
    MPI_Status ProbeStatus, RecvStatus;
    IntVect LCorner, UCorner;    
    char error_buffer[1000];int actual_length,rescode;
            
    rescode = MPI_Probe(0 /*source*/ , 0 /*tag*/, Chombo_MPI::comm, &ProbeStatus);    
    if (rescode!=MPI_SUCCESS)    
      MPI_Error_string(rescode,error_buffer,&actual_length);
      
    rescode = MPI_Get_count(&ProbeStatus,MPI_CH_REAL,&iCount);    
        
    Real* RealBuffer=new Real[iCount];
    rescode = MPI_Recv(RealBuffer,iCount,MPI_CH_REAL,0 /*source*/,0 /*tag*/,Chombo_MPI::comm,&RecvStatus);  
    CH_assert((ProbeStatus.MPI_SOURCE==RecvStatus.MPI_SOURCE)&&(ProbeStatus.MPI_TAG==RecvStatus.MPI_TAG));
        
    if (iCount > 1)  
    {
      LevelData<FArrayBox>& MySourceTerms = Level0->getSCData();
      posInRealBuffer=0;
      DataIterator dit = levelDomain.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)  
      {        
        const Box& b = levelDomain[dit()]; 
        FArrayBox& SFab = MySourceTerms[dit()];               
        Real* BoxData = &RealBuffer[posInRealBuffer];
        FArrayBox fab(b,numSourceTerms(),BoxData);                       
                         
        SFab.copy(fab);      
        ModifySourceTerms(SFab,b);
        
        posInRealBuffer+=b.numPts()*fab.nComp();      
      }            
      CH_assert(posInRealBuffer==iCount);
    }
    
    delete[] RealBuffer;
  }
#endif              

}

void KineticSources::CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep)
{
  if( m_verbosity >= 2 )
  {
    pout() << "KineticSources::CalculateSources" << endl;
  }
  
  int curProcId = procID();  
  
  Real cur_time = a_time;
  int cur_step  = a_curStep;
  
  if (checkTime(cur_time, cur_step) == false) return;  
  
  Vector<AMRLevel*>& amrlevels = a_levels;
  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(amrlevels[0]);
  
  const Box& Level0Box = amrlevels[0]->problemDomain().domainBox();
  FArrayBox Level0Data(Level0Box, numVarsForSC() );
  
  sendDataToZeroProcessor(Level0);
  GatherDataOnZeroProcessor(Level0,Level0Data);
  
  FArrayBox SourceTerms;    
  if (curProcId==0) SourceTerms.define(Level0Box,numSourceTerms()); 
  
  Real dx = Level0->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,0);
  
  CallExternalSC(Level0Data, Level0->dt(), dx, SourceTerms);
  
  BroadcastSourceTerms(Level0, SourceTerms);
   
      
  // Interpolate source terms on finer levels  
  InterpolateSourceTerms(a_levels);
   
}

// Interpolate source terms from zero level to higher levels.
void KineticSources::InterpolateSourceTerms(Vector<AMRLevel*> a_levels) 
{
  int level;
  Vector<AMRLevel*> amrlevels = a_levels;
  
  AMRLevelIdealMHD *amrMHDCoarserPtr,*curLevel;
  
  
  // Interpolate source terms on finer levels  
  for(level = 0; level < (int)amrlevels.size(); ++level)
  if (amrlevels[level]->boxes().size() > 0) 
  {                                         
    curLevel = static_cast<AMRLevelIdealMHD*>(amrlevels[level]);    
          
    amrMHDCoarserPtr = curLevel->getCoarserLevel(); 
    LevelData<FArrayBox>& SCData = curLevel->getSCData();
    LevelData<FArrayBox>& SCDataCoarse = amrMHDCoarserPtr->getSCData();
    
    const DisjointBoxLayout& grids = SCData.disjointBoxLayout();
  
    // Interpolate from coarser level 
    if (amrMHDCoarserPtr!=NULL)
    {
      int nRefCrse = amrMHDCoarserPtr->refRatio();
      
      FineInterp fineInterpST;  
      fineInterpST.define(grids,
                          SCData.nComp(),
                          nRefCrse,
                          curLevel->problemDomain());
                          
      fineInterpST.interpToFine(SCData,SCDataCoarse);
    }
  
    // Filling ghost cells
    if (SCData.ghostVect() != IntVect::Zero)
    {
      int numGhost = SCData.ghostVect()[0];
      if (amrMHDCoarserPtr!=NULL)
      {      
        const DisjointBoxLayout& levelDomain = SCData.disjointBoxLayout();
        
        PiecewiseLinearFillPatch pwl(levelDomain,
                                 SCDataCoarse.disjointBoxLayout(),
                                 SCData.nComp(),
                                 amrMHDCoarserPtr->problemDomain(),
                                 amrMHDCoarserPtr->refRatio(),
                                 numGhost);
  
        pwl.fillInterp(SCData,
                   SCDataCoarse,
                   SCDataCoarse,
                   1.0,
                   0,
                   0,
                   SCData.nComp());
      }
      
      SCData.exchange(Interval(0,SCData.nComp()-1));
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
