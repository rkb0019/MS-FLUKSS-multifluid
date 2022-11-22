#include "AMRLevelIdealMHD.H"


#include <list>
#include <map>
#include <vector>


#include "CH_HDF5.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "FArrayBox.H"
#include "DisjointBoxLayout.H"
#include "LoHiCenter.H"
#include <BoxIterator.H> 
#include <BaseFab.H> 

#include "PatchMHDAMF_F.H"
#include "PatchIdealMHDF_F.H"
#include <AuxVectors.h>
#include "SourceCalculator.H"
#include "tecplotF_F.H"
#include "LGintegrator.H"
#include "PiecewiseLinearFillPatchMHDAM.H"
#include "EquationSystem.H"
#include "EqSysMHDMF.H"
#include "MHDAMDefs.H"
#include "AMRLevelIdealMHD.H"
#include "DebugF_F.H"
#include "PatchMHDMF.H"

#ifdef CH_USE_HDF5

// Write checkpoint header
void AMRLevelIdealMHD::writeCheckpointHeader(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writeCheckpointHeader" << endl;
  }

  // Setup the number of components
  HDF5HeaderData header;
  header.m_int["num_components"] = m_numStates;

  // Setup the component names
  char compStr[30];
  for (int comp = 0; comp < m_numStates; ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = m_stateNames[comp];
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }
}

// Write checkpoint data for this level
void AMRLevelIdealMHD::writeCheckpointLevel(HDF5Handle& a_handle) const
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writeCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;
  
  Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);

  a_handle.setGroup(label);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]       = m_ref_ratio;
  header.m_int ["tag_buffer_size"] = m_tagBufferSize;
  header.m_real["dx"]              = dx;
  header.m_real["dt"]              = m_dt;
  header.m_box ["prob_domain"]     = m_problem_domain.domainBox();

  // Setup the periodicity info
  D_TERM(
         if (m_problem_domain.isPeriodic(0))
           header.m_int ["is_periodic_0"] = 1;
         else
           header.m_int ["is_periodic_0"] = 0; ,

         if (m_problem_domain.isPeriodic(1))
           header.m_int ["is_periodic_1"] = 1;
         else
           header.m_int ["is_periodic_1"] = 0; ,

         if (m_problem_domain.isPeriodic(2))
           header.m_int ["is_periodic_2"] = 1;
         else
           header.m_int ["is_periodic_2"] = 0; );

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  // Write the data for this level
  write(a_handle,m_UNew.boxLayout());
  write(a_handle,m_UNew,"data");
    
}

// Read checkpoint header
void AMRLevelIdealMHD::readCheckpointHeader(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::readCheckpointHeader" << endl;
  }

  // Reader the header
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the number of components
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointHeader: checkpoint file does not have num_components");
  }

  int numStates = header.m_int["num_components"];
  if (numStates < m_numStates)
  {
    //MayDay::Error("AMRLevelIdealMHD::readCheckpointHeader: num_components in checkpoint file does not match solver");
  }

  // Get the component names
  std::string stateName;
  char compStr[60];
  for (int comp = 0; comp < MIN(numStates,m_numStates); ++comp)
  {
    sprintf(compStr,"component_%d",comp);
    if (header.m_string.find(compStr) == header.m_string.end())
    {
      MayDay::Error("AMRLevelIdealMHD::readCheckpointHeader: checkpoint file does not have enough component names");
    }

    stateName = header.m_string[compStr];
    if (stateName != m_stateNames[comp])
    {
      pout() << "level " << m_level << ", stateName = " << stateName << ", m_stateNames[comp] = " << m_stateNames[comp];
      MayDay::Error("AMRLevelIdealMHD::readCheckpointHeader: state_name in checkpoint does not match solver");
    }
  }
}

// Read checkpoint data for this level
void AMRLevelIdealMHD::readCheckpointLevel(HDF5Handle& a_handle)
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::readCheckpointLevel" << endl;
  }

  // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  // Read the header for this level
  a_handle.setGroup(label);
  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << "hdf5 header data:" << endl;
    pout() << header << endl;
  }

  // Get the refinement ratio
  if (header.m_int.find("ref_ratio") == header.m_int.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain ref_ratio");
  }
  
  // SB 2/21/2011 I need to comment this. For the finest level in the checkpoint we must have ability to change ref_ratio since next level does not exist.
  // Otherwise it leads to troubles when we add new levels.
  // AMR_MHDAM::setupForRestart checks that input file is consistent with the checkpoint
  
  // m_ref_ratio = header.m_int["ref_ratio"];

  if (s_verbosity >= 2)
  {
    pout() << "read ref_ratio = " << m_ref_ratio << endl;
  }

  // Get the tag buffer size
  //if (header.m_int.find("tag_buffer_size") == header.m_int.end())
  //{
  //  MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain tag_buffer_size");
  //}
  //m_tagBufferSize=  header.m_int["tag_buffer_size"];

  if (s_verbosity >= 2)
  {
    pout() << "read tag_buffer_size = " << m_tagBufferSize << endl;
  }

  // Get dx
  if (header.m_real.find("dx") == header.m_real.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain dx");
  }

  if (m_csh->constStep(0))
  if (header.m_real["dx"]!=m_csh->dx(0,m_level))
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: checkpoint file inconsistent with inputs to define, 'dx' is mismatch");
  }
  

  if (s_verbosity >= 2)
  {
    //pout() << "read dx = " << m_dudvdw[0] << endl;
  }

  // Get dt
  if (header.m_real.find("dt") == header.m_real.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain dt");
  }
  m_dt = header.m_real["dt"];

  if (s_verbosity >= 2)
  {
    pout() << "read dt = " << m_dt << endl;
  }

  // Get the problem domain
  if (header.m_box.find("prob_domain") == header.m_box.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain prob_domain");
  }

  Box domainBox = header.m_box["prob_domain"];

  // Get the periodicity info -- this is more complicated than it really
  // needs to be in order to preserve backward compatibility
  bool isPeriodic[SpaceDim];
  D_TERM(if (!(header.m_int.find("is_periodic_0") == header.m_int.end()))
           isPeriodic[0] =  (header.m_int["is_periodic_0"] == 1);
         else
           isPeriodic[0] = false; ,

         if (!(header.m_int.find("is_periodic_1") == header.m_int.end()))
           isPeriodic[1] =  (header.m_int["is_periodic_1"] == 1);
         else
           isPeriodic[1] = false; ,

         if (!(header.m_int.find("is_periodic_2") == header.m_int.end()))
           isPeriodic[2] =  (header.m_int["is_periodic_2"] == 1);
         else
           isPeriodic[2] = false;);

  if (domainBox!=m_problem_domain.domainBox())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: checkpoint file inconsistent with inputs to define, problem domains are mismatch");
  }
  m_problem_domain = ProblemDomain(domainBox,isPeriodic);

  // Get the grids
  Vector<Box> grids;
  const int gridStatus = read(a_handle,grids);
  if (gridStatus != 0)
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain a Vector<Box>");
  }

#if CHOMBO_VERSION_MAJOR < 4    
  initialGrid(grids);
#else
  int size;
  ParmParse pp;
  pp.get("max_grid_size", size);
  defineDBLFromVBox(m_grids, grids, m_problem_domain, size);
  initialGrid(m_grids);
#endif  

  if (s_verbosity >= 4) if (0)
  {
    pout() << "read level domain: " << endl;
    LayoutIterator lit = m_grids.layoutIterator();
    for (lit.begin(); lit.ok(); ++lit)
    {
      const Box& b = m_grids[lit()];
#if CHOMBO_VERSION_MAJOR < 4          
      pout() << lit().intCode() << ": " << b << endl;
#else   
      pout() << ": " << b << endl;
#endif     
    }
    pout() << endl;
  }

  // num_components should be read from root
  a_handle.setGroup("/");
  header.readFromFile(a_handle);
  if (header.m_int.find("num_components") == header.m_int.end())
  {
    MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain num_components");
  }
  int num_components = header.m_int["num_components"];

  // Return to level header
  a_handle.setGroup(label);
  header.readFromFile(a_handle);

  /*if (s_verbosity >= 4) pout() << "define LevelData<FArrayBox> m_UNew" << endl;
  IntVect ivGhost = m_numGhost*IntVect::Unit;
  // Reshape state with new grids
  m_UNew.define(m_grids,m_numStates,ivGhost);*/

  if (num_components == m_numStates)
  {
    if (s_verbosity >= 4) pout() << "read LevelData<FArrayBox> m_UNew" << endl;
    const int dataStatus = read<FArrayBox>(a_handle,
                                           m_UNew,
                                           "data",
                                           m_grids,
                                           Interval(0,m_numStates-1),
                                           false);
    if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain state data");
    }
  } else
  {    
    // Create temporal buffer
    LevelData<FArrayBox> UBuffer(m_grids,num_components);
    const int dataStatus = read<FArrayBox>(a_handle,
                                           UBuffer,
                                           "data",
                                           m_grids);
    if (dataStatus != 0)
    {
      MayDay::Error("AMRLevelIdealMHD::readCheckpointLevel: file does not contain state data");
    }
    
    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
    DataIterator dit = m_UNew.dataIterator();
    
    if (num_components < m_numStates)
    {
      for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& U          = m_UNew[dit()];
        const FArrayBox& UBuf = UBuffer[dit()];
        const Box& b          = levelDomain[dit()];
        U.copy(UBuf,0,0,num_components);
      }    
      
      Interval extraComp(num_components,m_numStates-1);
      
      PhysProblem* PhysProblemPtr = m_patchMHDAM->getPhysProblem();
      PhysProblemPtr->initialize(m_UNew,extraComp);            
    } else
    {    
      for (dit.begin(); dit.ok(); ++dit)
      {
        FArrayBox& U          = m_UNew[dit()];
        const FArrayBox& UBuf = UBuffer[dit()];
        const Box& b          = levelDomain[dit()];
        U.copy(UBuf,0,0,m_numStates);
      }      
    }
    
  }
  
//Temporary change of the magnetic filed sign
//      DataIterator dit = m_UNew.dataIterator();
//      for (dit.begin(); dit.ok(); ++dit)
//      {
//        FArrayBox& U = m_UNew[dit()];
//        U.mult(-1.0,UBX,3);
//      }    


  // I do not trust 'dt' stored in a checkpoint. It is safer to recalculate it.  
  m_dt = computeNewDt(m_cfl);
  
  if ((s_verbosity >= 2) || ((s_verbosity <= -2)&&(procID() == 0)))
  {
    pout() << "readCheckpointLevel, level = " << m_level << " dt = " << m_dt << endl;
  }
  
  
  

  /*if (s_verbosity >= 4) pout() << "define LevelData<FArrayBox> m_UOld" << endl;
  m_UOld.define(m_grids,m_numStates,ivGhost);

  eDivergenceCleaning dc = m_patchMHDAM->getDivergenceCleaning();
  if (dc == DC_ProjHancock)
  {
    m_divB.define(m_grids,1,IntVect::Zero);
    m_phi.define(m_grids,1,IntVect::Unit);
  }
  if ((dc == DC_CT_BS)||(dc == DC_CT_GS0)||(dc == DC_CT_GS1))
  {
    m_E   .define( m_grids, 1 );
    m_EAcc.define( m_grids, 1 );
  }

  IntVect STGhost = m_numSCGhost*IntVect::Unit;
  if (m_patchMHDAM->getSourceCalculator()!=NULL)
  if (m_patchMHDAM->getSourceCalculator()->numSourceTerms()>0)
  {
    m_SCData.define(m_grids,m_patchMHDAM->getSourceCalculator()->numSourceTerms(), STGhost);
    for(DataIterator dit = m_SCData.dataIterator();dit.ok(); ++dit) m_SCData[dit()].setVal(0.0);
  }  
  
  if (getFinerLevel()!=NULL)
  {
    m_invVolumes.define(m_grids,1,IntVect::Zero);
  }

  // Set up data structures
  levelSetup();*/
}

// Write plotfile header
void AMRLevelIdealMHD::writePlotHeader(HDF5Handle& a_handle) const
{
  CH_TIME("writePlotHeader"); 
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writePlotHeader" << endl;
  }

  int numVars = m_eqSys->numPrimitives() +
                m_patchMHDAM->getPhysProblem()->numPlotVars() +
                m_patchMHDAM->getPhysProblem()->getSourceCalculator()->numPlotVars();

  Vector<string> varNames;
  varNames.reserve(numVars);


    // Setup the number of components -- include space for error
  HDF5HeaderData header;
  header.m_int["num_components"] = numVars;

  // Setup the component names
  char compStr[30];int comp = 0,i;
  varNames = m_eqSys->primitiveNames();
  for (i = 0; i < varNames.size(); ++i)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = varNames[i];
    comp++;
  }
  varNames = m_patchMHDAM->getPhysProblem()->plotNames();
  for (i = 0; i < varNames.size(); ++i)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = varNames[i];
    comp++;
  }
  varNames = m_patchMHDAM->getPhysProblem()->getSourceCalculator()->plotNames();
  for (i = 0; i < varNames.size(); ++i)
  {
    sprintf(compStr,"component_%d",comp);
    header.m_string[compStr] = varNames[i];
    comp++;
  }

  // Write the header
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }


}

// Write plotfile data for this level
void AMRLevelIdealMHD::writePlotLevel(HDF5Handle& a_handle) const
{
  CH_TIME("writePlotLevel"); 
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writePlotLevel" << endl;
  }
  
  {
  CH_TIME("writePlotLevel_otherstuff"); 

    // Setup the level string
  char levelStr[20];
  sprintf(levelStr,"%d",m_level);
  const std::string label = std::string("level_") + levelStr;

  a_handle.setGroup(label);
  
  Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);

  // Setup the level header information
  HDF5HeaderData header;

  header.m_int ["ref_ratio"]   = m_ref_ratio;
  header.m_real["dx"]          = dx;
  header.m_real["dt"]          = m_dt;
  header.m_box ["prob_domain"] = m_problem_domain.domainBox();

  // Write the header for this level
  header.writeToFile(a_handle);

  if (s_verbosity >= 3)
  {
    pout() << header << endl;
  }

  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();

  int numVars = m_eqSys->numPrimitives() +
                PhPr->numPlotVars() +
                SC->numPlotVars();

  // Write the data for this level
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  LevelData<FArrayBox> plotData(levelDomain, numVars);

  int iOffset;  

  DataIterator dit = plotData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    iOffset = 0;
    const FArrayBox& U       = m_UNew[dit()];
    FArrayBox& plotFab = plotData[dit()];
    const Box&       plotBox = levelDomain[dit()];

#ifndef NDEBUG
    int LCorner[CH_SPACEDIM]={D_DECL(plotBox.smallEnd()[0], plotBox.smallEnd()[1], plotBox.smallEnd()[2])};
    int UCorner[CH_SPACEDIM]={D_DECL(plotBox.bigEnd()[0],   plotBox.bigEnd()[1],   plotBox.bigEnd()[2])};
#endif


    FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
    m_eqSys->stateToPrim(W, U, plotBox);
    iOffset += m_eqSys->numPrimitives();

    if (PhPr->numPlotVars() > 0)
    {
      //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
      PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
    }    
    
    iOffset += PhPr->numPlotVars();

    if (SC->numPlotVars() > 0)
    {
      CH_assert(SC->numPlotVars() == m_SCData.nComp());
      const FArrayBox& SCData = m_SCData[dit()];
      plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
    }
    
    PhPr->primForPlot(W,plotBox);
    if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);

  }
  
  


  {
    CH_TIME("writeboxLayout");
    write(a_handle,m_UNew.boxLayout());
  }
  
  {
    CH_TIME("writedata");
    write(a_handle,plotData,"data");  
  }
  
  }
    

}
#endif

#ifdef CH_MPI  
void AMRLevelIdealMHD::writeLevelforTecPlot(MPI_File tecplot_file)
#else
void AMRLevelIdealMHD::writeLevelforTecPlot(FILE* tecplot_file)
#endif  
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writeLevelforTecPlot" << endl;
  }
  
    
  int D_DECL(i,j,k), ivar; char buffer[1000];
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();
  
  int numVars = m_eqSys->numPrimitives() + 
                PhPr->numPlotVars() +
                SC->numPlotVars();
                
  if (m_output_divB == true)
  {
    numVars++;
    //Interval comps(UBX,UBZ);    
    //fillGhostCellsCons(comps);
  }
  
  if (SpaceDim == 3)
  {
    Interval comps(0,m_eqSys->numStates()-1);    
    fillGhostCellsCons(comps);
  }
                
  
  if (m_level == 0)
  {
  #ifdef CH_MPI  
    MPI_Barrier(Chombo_MPI::comm);
    if (procID() == 0)
  #endif  
    {
      if (SpaceDim == 2) 
        strcpy(buffer,"VARIABLES=\"X\" \"Y\"");
      else if (SpaceDim == 3) 
        strcpy(buffer,"VARIABLES=\"X\" \"Y\" \"Z\"");
      else MayDay::Error("AMRLevelIdealMHD::writeLevelforTecPlot unsuported dimension");
                  
      Vector<string> varNames;
      varNames.reserve(numVars);        
                            
      varNames = m_eqSys->primitiveNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");      
      }
      varNames = m_patchMHDAM->getPhysProblem()->plotNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");      
      }
      varNames = m_patchMHDAM->getPhysProblem()->getSourceCalculator()->plotNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");      
      }            
      if (m_output_divB == true)
      {
        strcat(buffer," \"divb\"");        
      }
      strcat(buffer,"\n");      
      
  #ifdef CH_MPI  
      MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
  #else
      fprintf(tecplot_file,"%s",buffer);     
  #endif    
    }
  #ifdef CH_MPI    
    MPI_Barrier(Chombo_MPI::comm);
  #endif      
  }
  
    
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  LevelData<FArrayBox> plotData(levelDomain, numVars);

  int iOffset;
  //FArrayBox localFAB;

  Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);
    

  DataIterator dit = plotData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    iOffset = 0;
    const FArrayBox& U       = m_UNew[dit()];
    Box              plotBox = levelDomain[dit()];
    
    if (SpaceDim == 3)
    {
      // I want smooth isosurfaces without cracks
      plotBox.grow(1); // 1 ghost cell
      plotBox &= m_problem_domain;
    }
    
    FArrayBox plotFab(plotBox,numVars);
    
    // Calculation of variables that are necessary for plot
    FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
    m_eqSys->stateToPrim(W, U, plotBox);
    iOffset += m_eqSys->numPrimitives();

    if (PhPr->numPlotVars() > 0)
    {
      //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
      PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
    }
    iOffset += PhPr->numPlotVars();

    if (SC->numPlotVars() > 0)
    {
      CH_assert(SC->numPlotVars() == m_SCData.nComp());
      const FArrayBox& SCData = m_SCData[dit()];
      plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
    }
    
    PhPr->primForPlot(W,plotBox);
    if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);
      
    if (m_output_divB == true)
    {
      iOffset = plotFab.nComp() - 1;
      FArrayBox divB(plotBox,1);
      (static_cast<PatchMHDMF*>(m_patchMHDAM))->computeDivB(divB, U, plotBox);                            
      plotFab.copy(divB, plotBox, 0, plotBox, iOffset, 1);                              
    }
    
    
    std::string buffer_main;
    buffer_main.reserve(15*(numVars+CH_SPACEDIM)*plotBox.numPts()+100); // Reserve big enough chunk of memory    
    
    Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);
      
    D_TERM(  
      int nx=plotBox.size(0);,
      int ny=plotBox.size(1);,
      int nz=plotBox.size(2););

    IntVect LCorner=plotBox.smallEnd();
    IntVect UCorner=plotBox.bigEnd();
      

  #if CH_SPACEDIM == 2   
    sprintf(buffer,"ZONE T=L%i, I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([3-%i]=CELLCENTERED)\n", m_level, nx+1, ny+1, numVars+CH_SPACEDIM);
  #endif  
  #if CH_SPACEDIM == 3   
    sprintf(buffer,"ZONE T=L%i, I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([4-%i]=CELLCENTERED)\n", m_level, nx+1, ny+1, nz+1, numVars+CH_SPACEDIM);
  #endif
    buffer_main+=buffer;
    
    
    Real D_DECL(x,y,z);RealVect xyz;
    
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    {
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        x = xyz[0];
        sprintf(buffer,"%.6e ",x);                  
        buffer_main+=buffer;
      }  
      buffer_main+="\n";
    }
        
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif  
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    {
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        y = xyz[1];
        sprintf(buffer,"%.6e ",y);          
        buffer_main+=buffer;
      }
      buffer_main+="\n";          
    }

  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    for (i=LCorner[0];i<=UCorner[0]+1;i++)
    {          
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        z = xyz[2];
        sprintf(buffer,"%.6e ",z);          
        buffer_main+=buffer;
      }    
      buffer_main+="\n";          
    }
  #endif    
      
      
    IntVect aux_IntVect; Real value;
  
    for (ivar=0;ivar<numVars;ivar++)
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2];k++)
  #endif  
    for (j=LCorner[1];j<=UCorner[1];j++)               
    {
      for (i=LCorner[0];i<=UCorner[0];i++)
      { 
        D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
        value=plotFab.get(aux_IntVect,ivar);      
        sprintf(buffer,"%.6e ",value);           
        buffer_main+=buffer;
      }
      buffer_main+="\n";
    }
    
#ifdef CH_MPI
    char *buffer_to_write = new char[buffer_main.size()+1];
    buffer_main.copy(buffer_to_write, std::string::npos);
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(tecplot_file,"%s",buffer_main.c_str());
    fprintf(tecplot_file,"# Level %i ended\n",m_level);
#endif  
  }
  
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
  if (procID() == 0)
  {
    sprintf(buffer,"# Level %i ended\n",m_level);				
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
  }  
  MPI_Barrier(Chombo_MPI::comm);
#endif

}


#ifdef CH_MPI  
void AMRLevelIdealMHD::writeLevel0forTecPlot(MPI_File tecplot_file)
#else
void AMRLevelIdealMHD::writeLevel0forTecPlot(FILE* tecplot_file)
#endif  
{
  if (s_verbosity >= 3)
  {
    pout() << "AMRLevelIdealMHD::writeLevel0forTecPlot" << endl;
  }
  
#ifndef CH_MPI    
  writeLevelforTecPlot(tecplot_file);
  return;
#endif  

  const int maxNumbersPerLine = 100;
  int       iNumbersInLine;
  
#ifdef CH_MPI        
  MPI_Offset D_DECL(i,j,k), ivar; MPI_Offset offs;  
  char buffer[1000];
  
  std::string buffer_main;char *buffer_to_write;
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();
  
  Box pd_box = m_problem_domain.domainBox();
  
  int numVars = m_eqSys->numPrimitives() + 
                PhPr->numPlotVars() +
                SC->numPlotVars();
                
  if (m_output_divB == true)
  {
    numVars++;
    Interval comps(UBX,UBZ);    
    fillGhostCellsCons(comps);
  }
                                    
  
  MPI_Barrier(Chombo_MPI::comm);
  
  if (procID() == 0)
  {
    if (SpaceDim == 2) 
      strcpy(buffer,"VARIABLES=\"X\" \"Y\"");
    else if (SpaceDim == 3) 
      strcpy(buffer,"VARIABLES=\"X\" \"Y\" \"Z\"");
    else MayDay::Error("AMRLevelIdealMHD::writeLevelforTecPlot unsuported dimension");
                
    Vector<string> varNames;
    varNames.reserve(numVars);        
                          
    varNames = m_eqSys->primitiveNames();
    for (i = 0; i < (int)varNames.size(); ++i)          
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    varNames = m_patchMHDAM->getPhysProblem()->plotNames();
    for (i = 0; i < (int)varNames.size(); ++i)
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    varNames = m_patchMHDAM->getPhysProblem()->getSourceCalculator()->plotNames();
    for (i = 0; i < (int)varNames.size(); ++i)
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    if (m_output_divB == true)
    {
      strcat(buffer," \"divb\"");        
    }
    strcat(buffer,"\n");      

    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);

    // Output header and coordinates

    buffer_main.reserve(15*pd_box.size(0)*pd_box.size(1)+100); // Reserve big enough chunk of memory    
            
    D_TERM(  
      int nx=pd_box.size(0);,
      int ny=pd_box.size(1);,
      int nz=pd_box.size(2););

    IntVect LCorner=pd_box.smallEnd();
    IntVect UCorner=pd_box.bigEnd();
    
  #if CH_SPACEDIM == 2       
    if (numVars == 1)
      sprintf(buffer,"ZONE T=L0 I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([3]=CELLCENTERED)\n", nx+1, ny+1);
    else 
      sprintf(buffer,"ZONE T=L0 I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([3-%i]=CELLCENTERED)\n", nx+1, ny+1, numVars+CH_SPACEDIM);
  #endif  
  #if CH_SPACEDIM == 3   
    if (numVars == 1)
      sprintf(buffer,"ZONE T=L0 I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([4]=CELLCENTERED)\n", nx+1, ny+1, nz+1);
    else 
      sprintf(buffer,"ZONE T=L0 I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([4-%i]=CELLCENTERED)\n", nx+1, ny+1, nz+1, numVars+CH_SPACEDIM);
  #endif
    buffer_main+=buffer;
    
    
    Real D_DECL(x,y,z);RealVect xyz;
    
    iNumbersInLine = 0;
  
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif
    for (j=LCorner[1];j<=UCorner[1]+1;j++)    
    {
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        x = xyz[0];
        sprintf(buffer,"%.6e ",x);                  
        buffer_main+=buffer;
        iNumbersInLine++;
        if (iNumbersInLine > maxNumbersPerLine)
        {
          iNumbersInLine = 0;
          buffer_main+="\n";
        }        
      }        
      buffer_main+="\n";
    }
    buffer_to_write =  const_cast<char*> (buffer_main.c_str());            
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);    
    buffer_main.clear();

    iNumbersInLine = 0;    
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif  
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    {
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        y = xyz[1];
        sprintf(buffer,"%.6e ",y);          
        buffer_main+=buffer;
        iNumbersInLine++;
        if (iNumbersInLine > maxNumbersPerLine)
        {
          iNumbersInLine = 0;
          buffer_main+="\n";
        }
      }
      buffer_main+="\n";          
    }

    buffer_to_write =  const_cast<char*> (buffer_main.c_str());            
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);    
    buffer_main.clear();

    iNumbersInLine = 0;
#if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    for (i=LCorner[0];i<=UCorner[0]+1;i++)
    {          
      for (i=LCorner[0];i<=UCorner[0]+1;i++)
      {
        m_csh->getNodeCoordsCartesian(xyz, IntVect(D_DECL(i,j,k)), m_level);
        z = xyz[2];
        sprintf(buffer,"%.6e ",z);          
        buffer_main+=buffer;
        iNumbersInLine++;
        if (iNumbersInLine > maxNumbersPerLine)
        {
          iNumbersInLine = 0;
          buffer_main+="\n";
        }
      }    
      buffer_main+="\n";          
    }
    buffer_to_write =  const_cast<char*> (buffer_main.c_str());            
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);    
    buffer_main.clear();
#endif    
  }
              
  MPI_Barrier(Chombo_MPI::comm);
  MPI_File_seek(tecplot_file,0,MPI_SEEK_END);
  MPI_File_get_position(tecplot_file, &offs);
  MPI_Barrier(Chombo_MPI::comm);
      
    
  D_TERM(  
        MPI_Offset nxpd=pd_box.size(0);,
        MPI_Offset nypd=pd_box.size(1);,
        MPI_Offset nzpd=pd_box.size(2););
    
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  LevelData<FArrayBox> plotData(levelDomain, numVars);

  int iOffset; MPI_Offset blocksize = 15;
  //FArrayBox localFAB;
  
  iNumbersInLine = 0;
    
  // Prepare data for output
  DataIterator dit = plotData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    iOffset = 0;
    const FArrayBox& U       = m_UNew[dit()];
    const Box&       plotBox = levelDomain[dit()];
    FArrayBox & plotFab = plotData[dit()];
    
    // Calculation of variables that are necessary for plot
    FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
    m_eqSys->stateToPrim(W, U, plotBox);
    iOffset += m_eqSys->numPrimitives();

    if (PhPr->numPlotVars() > 0)
    {
      //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
      PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
    }
    iOffset += PhPr->numPlotVars();

    if (SC->numPlotVars() > 0)
    {
      CH_assert(SC->numPlotVars() == m_SCData.nComp());
      const FArrayBox& SCData = m_SCData[dit()];
      plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
    }    
    PhPr->primForPlot(W,plotBox);
    if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);
      
    if (m_output_divB == true)
    {
      iOffset = plotFab.nComp() - 1;
      FArrayBox divB(plotBox,1);
      (static_cast<PatchMHDMF*>(m_patchMHDAM))->computeDivB(divB, U, plotBox);                            
      plotFab.copy(divB, plotBox, 0, plotBox, iOffset, 1);                              
    }

    
    //  This works, but TOO slow. 
  
/*    IntVect LCorner=plotBox.smallEnd();
    IntVect UCorner=plotBox.bigEnd();
    
    MPI_Offset blocksize = 14;
    
    IntVect aux_IntVect; Real value;
    
#if CH_SPACEDIM == 3        
    for (ivar=0;ivar<numVars;ivar++)
    {  
      for (k=LCorner[2];k<=UCorner[2];k++)
      for (j=LCorner[1];j<=UCorner[1];j++)               
      {
        for (i=LCorner[0];i<=UCorner[0];i++)
        { 
          D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
          value=plotFab.get(aux_IntVect,ivar);      
          sprintf(buffer,"%13.6e ",value);           
          buffer_main+=buffer;          
        }        
        buffer_main[buffer_main.length()-1]='\n';
        buffer_to_write = const_cast<char*>(buffer_main.c_str());                            
              
        MPI_File_write_at(tecplot_file, offs + 
            ivar*nxpd*nypd*nzpd*blocksize + k*nxpd*nypd*blocksize + j*nxpd*blocksize + LCorner[0]*blocksize, 
            buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);        
        buffer_main.clear();
      }                  
    }
#endif      

#if CH_SPACEDIM == 2
    for (ivar=0;ivar<numVars;ivar++)
    {        
      for (j=LCorner[1];j<=UCorner[1];j++)               
      {
        for (i=LCorner[0];i<=UCorner[0];i++)
        { 
          D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
          value=plotFab.get(aux_IntVect,ivar);      
          sprintf(buffer,"%13.6e ",value);           
          buffer_main+=buffer;          
        }        
        buffer_main[buffer_main.length()-1]='\n';
        buffer_to_write = const_cast<char*>(buffer_main.c_str());                            
               
        MPI_File_write_at(tecplot_file, offs + ivar*nxpd*nypd*blocksize + j*nxpd*blocksize + LCorner[0]*blocksize, 
            buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);        
        buffer_main.clear();
      }                  
    }
#endif      
  */      
  }
    
            
  int iLoDomain = m_problem_domain.domainBox().smallEnd()[SpaceDim-1];
  int iHiDomain = m_problem_domain.domainBox().bigEnd()[SpaceDim-1];
  
  Vector<Box> vstrips;
  
  for (j=iLoDomain;j<=iHiDomain;++j)
  {
    Box b = pd_box;
    b.setSmall(SpaceDim-1,j);
    b.setBig(SpaceDim-1,j);
    vstrips.push_back(b);
  }
  
  Vector<int> procs;
  LoadBalance(procs, vstrips);
  DisjointBoxLayout strips(vstrips, procs, m_problem_domain);  
  LevelData<FArrayBox> kdata(strips,1);     
                      
  if (SpaceDim == 2)
    buffer_main.reserve(nxpd*blocksize+20);      // Reserve big enough chunk of memory            
  else    
    buffer_main.reserve(nypd*nxpd*blocksize+20); // Reserve big enough chunk of memory              
    
  IntVect aux_IntVect; Real value; 

  for (ivar=0;ivar<numVars;ivar++)  
  {
    plotData.copyTo(Interval(ivar,ivar),kdata,Interval(0,0));    
                            
    DataIterator dit = kdata.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      const FArrayBox& plotFab = kdata[dit()];
      const Box&       plotBox = plotFab.box();
      int kind = plotBox.smallEnd(SpaceDim - 1);
      
      IntVect LCorner=plotBox.smallEnd();
      IntVect UCorner=plotBox.bigEnd();

#if CH_SPACEDIM == 3                      
      for (j=LCorner[1];j<=UCorner[1];j++)               
#elif CH_SPACEDIM == 2
      j = kind;
#endif
      {      
        for (i=LCorner[0];i<=UCorner[0];i++)
        { 
          D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=kind;);
          value=plotFab.get(aux_IntVect,0);      
          sprintf(buffer,"%14.6e ",value);           
          buffer_main+=buffer;    
          iNumbersInLine++;
          if (iNumbersInLine > maxNumbersPerLine)
          {
            iNumbersInLine = 0;
            buffer_main[buffer_main.length()-1]='\n';
          }      
        }       
        buffer_main[buffer_main.length()-1]='\n';
      }      
      
      buffer_to_write = const_cast<char*>(buffer_main.c_str());  
            
#if CH_SPACEDIM == 3                            
      MPI_Offset offs_b = ivar*nzpd*nypd*nxpd*blocksize + kind*nypd*nxpd*blocksize;
#elif CH_SPACEDIM == 2
      MPI_Offset offs_b = ivar*nypd*nxpd*blocksize + kind*nxpd*blocksize;
#endif

//      pout() << "offs = " << offs << ", offs_b=" << offs_b << ", length = " << buffer_main.length() << " " << buffer_main;
      
      MPI_File_write_at(tecplot_file, offs + offs_b, buffer_to_write, buffer_main.length(), MPI_CHAR, MPI_STATUS_IGNORE);        
      buffer_main.clear();
    }    
  }    
  

  MPI_Barrier(Chombo_MPI::comm);
  MPI_File_set_view(tecplot_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
  MPI_File_seek_shared(tecplot_file,0,MPI_SEEK_END);
  
#endif

}



// LZ: Level per Zone (one level = one zone)
void AMRLevelIdealMHD::writeLevelforTecPlot_LZ(FILE* tecplot_file,int a_nComp)
{

    if (s_verbosity >= 3)
    {
        pout() << "AMRLevelIdealMHD::writeLevel for TecPlot LZ" << endl;
    }

    Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);

    int ivar;
    int D_DECL(i,j,k);

    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
    DataIterator it = levelDomain.dataIterator();


    if (m_hasCoarser)
    {
        const AMRLevelIdealMHD* amrMHDCoarserPtr = getCoarserLevel();

        if (m_output_density_gradient)
        {
            PiecewiseLinearFillPatchMHDAM pwl_main(levelDomain,
                                     amrMHDCoarserPtr->m_UNew.disjointBoxLayout(),
                                     m_numStates,
                                     amrMHDCoarserPtr->m_problem_domain,
                                     amrMHDCoarserPtr->m_ref_ratio,
                                     m_level-1,
                                     1,
                                     m_csh);

            pwl_main.fillInterp(m_UNew,
                       amrMHDCoarserPtr->m_UNew,
                       amrMHDCoarserPtr->m_UNew,
                       1.0,
                       0,
                       0,
                       m_numStates);
        }


    }

    if (m_output_density_gradient)   m_UNew.exchange(Interval(0,m_numStates-1));


    bool hasValidFiner=false;
    if (m_hasFiner)
    {
        hasValidFiner = getFinerLevel()->m_UNew.isDefined();
    }

    std::list<IntVect> notOverlappedCells;
    std::list<IntVect>::iterator list_iter;

    std::list<Real> output_data_list;
    std::vector<Real> output_data_vector;

    IntVect aux_IntVect, tmp_iv; Real value;

    int num_variables=m_eqSys->numStates();
    if (m_output_density_gradient)   num_variables++;
    if (m_output_B_gradient) num_variables++;

    num_variables = a_nComp;

    for( it.begin(); it.ok(); ++it)
    {
        const FArrayBox& U = m_UNew[it()];        
        const Box& box = levelDomain[it()];
        
        FArrayBox data(box, m_eqSys->numPrimitives());        
        m_eqSys->stateToPrim(data, U, box);

        int D_DECL(nx=box.size(0), ny=box.size(1), nz=box.size(2));

        IntVect LCorner=box.smallEnd();
        IntVect UCorner=box.bigEnd();

        D_TERM(
        CH_assert(nx==(UCorner[0]-LCorner[0]+1));,
        CH_assert(ny==(UCorner[1]-LCorner[1]+1));,
        CH_assert(nz==(UCorner[2]-LCorner[2]+1)););

        FArrayBox* gradDensityMagFab=NULL;
        FArrayBox* gradBMagFab=NULL;

        if( m_output_density_gradient )
        {
            FArrayBox gradDensityFab(box,SpaceDim);
            const FArrayBox& UFab = m_UNew[it()];

            for (int dir = 0; dir < SpaceDim; ++dir)
            {
                const Box bCenter = box & grow(m_problem_domain,-BASISV(dir));
                const Box bLo     = box & adjCellLo(bCenter,dir);
                const int hasLo = ! bLo.isEmpty();
                const Box bHi     = box & adjCellHi(bCenter,dir);
                const int hasHi = ! bHi.isEmpty();

                    FORT_GETRELGRAD(CHF_FRA1(gradDensityFab,dir),
                      CHF_CONST_FRA1(UFab,0),
                      CHF_CONST_INT(dir),
                      CHF_BOX(bLo),
                      CHF_CONST_INT(hasLo),
                      CHF_BOX(bHi),
                      CHF_CONST_INT(hasHi),
                      CHF_BOX(bCenter));									  						
            }
            gradDensityMagFab=new FArrayBox(box,1);			
            FORT_MAGNITUDE(CHF_FRA1((*gradDensityMagFab),0),
                           CHF_CONST_FRA(gradDensityFab),
                           CHF_BOX(box));
        }

/*        if( m_output_B_gradient )
        {
            FArrayBox gradBFab(box,SpaceDim);
            const FArrayBox& UExtraFab = m_UExtra[it()];

            const FArrayBox& UFab = m_UNew[it()];

            for (int dir = 0; dir < SpaceDim; ++dir)
            {
                const Box bCenter = box & grow(m_problem_domain,-BASISV(dir));
                const Box bLo     = box & adjCellLo(bCenter,dir);
                const int hasLo = ! bLo.isEmpty();
                const Box bHi     = box & adjCellHi(bCenter,dir);
                const int hasHi = ! bHi.isEmpty();

                //Box loBox, hiBox, centerBox, entireBox;
                //int hasLo, hasHi;

                //loHiCenter( loBox, hasLo, hiBox, hasHi, centerBox, entireBox,
                //box, m_problem_domain, dir );

               // Compute each component of the divergence of the Magnetic Field
                //FORT_GETGRAD(CHF_FRA1(gradBFab,dir),
                // CHF_CONST_FRA1(UFab,UBX+dir),
                // CHF_CONST_INT(dir),
                // CHF_BOX(loBox),
                // CHF_CONST_INT(hasLo),
                // CHF_BOX(hiBox),
                // CHF_CONST_INT(hasHi),
                // CHF_BOX(centerBox));



//                    FORT_GETBRELGRAD(CHF_FRA1(gradBFab,dir),
//                      CHF_CONST_FRA1(UExtraFab,0),
//                      CHF_CONST_INT(dir),
//                      CHF_BOX(bLo),
//                      CHF_CONST_INT(hasLo),
//                      CHF_BOX(bHi),
//                      CHF_CONST_INT(hasHi),
//                      CHF_BOX(bCenter));

                    FORT_GETGRAD(CHF_FRA1(gradBFab,dir),
                      CHF_CONST_FRA1(UFab,UBX + dir),
                      CHF_CONST_INT(dir),
                      CHF_BOX(bLo),
                      CHF_CONST_INT(hasLo),
                      CHF_BOX(bHi),
                      CHF_CONST_INT(hasHi),
                      CHF_BOX(bCenter));

            }
            gradBMagFab=new FArrayBox(box,1);	
//            FORT_MAGNITUDE(CHF_FRA1((*gradBMagFab),0),
//                CHF_CONST_FRA(gradBFab),
//                CHF_BOX(box));

            FORT_GETDIVBSTEP2(CHF_FRA1((*gradBMagFab),0),
                CHF_CONST_FRA(gradBFab),
                CHF_BOX(box),
                CHF_CONST_REAL(dx));
        }*/
        FArrayBox* data_additional[2]={gradDensityMagFab,gradBMagFab};

        if (hasValidFiner==false)
        {
            // It is simple structured mesh. All boxes are visible.
            D_TERM(
            for (i=LCorner[0];i<=UCorner[0];i++),
            for (j=LCorner[1];j<=UCorner[1];j++),
            for (k=LCorner[2];k<=UCorner[2];k++))
            {
                D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
                notOverlappedCells.push_back(aux_IntVect);

                for (ivar=0;ivar<num_variables;ivar++)
                if (ivar<m_eqSys->numStates())
                {
                    value=data.get(aux_IntVect,ivar);
                    output_data_list.push_back(value);
                }

                // Writing gradients if needed
                for (ivar=0;ivar<2;ivar++)
                if (data_additional[ivar]!=NULL)
                {
                    value=data_additional[ivar]->get(aux_IntVect,0);
                    output_data_list.push_back(value);
                }
            }
        } else
        {
            Box scaledBox = box;
            scaledBox.refine(refRatio());

            std::list<Box> finerBoxes; // Boxes of finer level that have nonzero intersection with this box
            std::list<Box>::iterator bIter;

            DataIterator it_internal;
            for(it_internal = getFinerLevel()->m_UNew.dataIterator(); it_internal.ok(); ++it_internal)
            {
                const Box& box_internal = getFinerLevel()->m_UNew.box(it_internal());
                if (scaledBox.intersects(box_internal))
                  finerBoxes.push_back(box_internal);
            }				

            IntVect shifted_IntVect;
            Box bInternal;

            D_TERM(
            for (i=LCorner[0];i<=UCorner[0];i++),
            for (j=LCorner[1];j<=UCorner[1];j++),
            for (k=LCorner[2];k<=UCorner[2];k++))
            {
                D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
                shifted_IntVect = aux_IntVect;
                shifted_IntVect.scale(refRatio());

                bool CellIsOverlapped=false;

                for (bIter=finerBoxes.begin();bIter!=finerBoxes.end();++bIter)
                {
                  bInternal = (*bIter);
                  if (bInternal.contains(shifted_IntVect))
                  {					
                    CellIsOverlapped=true;
                    break;
                  }				
                }

                //for(it_internal = getFinerLevel()->m_UNew.dataIterator(); it_internal.ok(); ++it_internal)
                //{
                //    const Box& box_internal = getFinerLevel()->m_UNew.box(it_internal());
                //    if (box_internal.contains(shifted_IntVect))
                //    {					
                //        CellIsOverlapped=true;
                //        break;
                //    }				
                //}	
			
                if (CellIsOverlapped) continue;

                notOverlappedCells.push_back(aux_IntVect);
                for (ivar=0;ivar<num_variables;ivar++)
                if (ivar<m_eqSys->numStates())
                {
                    value=data.get(aux_IntVect,ivar);
                    output_data_list.push_back(value);
                }

                // Writing gradients if needed
                for (ivar=0;ivar<2;ivar++)
                if (data_additional[ivar]!=NULL)
                {
                    value=data_additional[ivar]->get(aux_IntVect,0);
                    output_data_list.push_back(value);
                }
            }
        }
        if (gradDensityMagFab!=NULL) delete gradDensityMagFab;
        if (gradBMagFab!=NULL) delete gradBMagFab;
    }

    // Copying all data to "output_data_vector"
    // Of course it would be possible to use "output_data_vector" as a list so auxiliary variable output_data_list wouldn't require.
    // But std::vector makes a lot of reallocations. So "output_data_list" was introduced to avoid numerous reallocations.
    std::list<Real>::iterator Real_iter;
    output_data_vector.resize(output_data_list.size()+1,0.0);
    i=0;
    for (Real_iter=output_data_list.begin();Real_iter!=output_data_list.end();++Real_iter)
    {
        output_data_vector[i]=(*Real_iter);
        i++;
    }

    // Get unique map of points from the cells that form mesh on this level.
    int_vector3d o_vec3d(0,0,0);

    std::map<int_vector3d, int> map_of_points;
    std::map<int_vector3d, int>::iterator map_of_points_iter;

    std::list<Real> D_DECL(x_coords,y_coords,z_coords);
    RealVect xyz;        

    if (notOverlappedCells.size()>0)
    {
        IntVect D_DECL(iv0(IntVect::Zero),iv1(IntVect::Zero),iv2(IntVect::Zero));
        D_TERM(iv0[0]=1;,iv1[1]=1;,iv2[2]=1);
        
        int TecPlot_PointIndex=0;
        for (list_iter=notOverlappedCells.begin();list_iter!=notOverlappedCells.end();++list_iter)
#if CH_SPACEDIM == 3
        for (k=0;k<2;k++)
#endif
        {
          aux_IntVect=(*list_iter);
          tmp_iv=IntVect(D_DECL(aux_IntVect[0],aux_IntVect[1],aux_IntVect[2]+k));
          
          o_vec3d.set_data(D_DECL(tmp_iv[0],tmp_iv[1],tmp_iv[2]));
          map_of_points_iter=map_of_points.find(o_vec3d);
          if (map_of_points_iter==map_of_points.end())
          {
            map_of_points[o_vec3d]=++TecPlot_PointIndex; //eq. "map_of_points[o_vec3d]=TecPlot_PointIndex+1";
            m_csh->getNodeCoordsCartesian(xyz, tmp_iv, m_level);                              
            D_TERM(
              x_coords.push_back(xyz[0]);,
              y_coords.push_back(xyz[1]);,					
              z_coords.push_back(xyz[2]));					
          }
          
          tmp_iv=IntVect(D_DECL(aux_IntVect[0]+1,aux_IntVect[1],aux_IntVect[2]+k));
          o_vec3d.set_data(D_DECL(tmp_iv[0],tmp_iv[1],tmp_iv[2]));
          map_of_points_iter=map_of_points.find(o_vec3d);
          if (map_of_points_iter==map_of_points.end())
          {
            map_of_points[o_vec3d]=++TecPlot_PointIndex;
            m_csh->getNodeCoordsCartesian(xyz, tmp_iv, m_level);                              
            D_TERM(
              x_coords.push_back(xyz[0]);,
              y_coords.push_back(xyz[1]);,					
              z_coords.push_back(xyz[2]));					

          }

          tmp_iv=IntVect(D_DECL(aux_IntVect[0]+1,aux_IntVect[1]+1,aux_IntVect[2]+k));
          o_vec3d.set_data(D_DECL(tmp_iv[0],tmp_iv[1],tmp_iv[2]));
          map_of_points_iter=map_of_points.find(o_vec3d);
          if (map_of_points_iter==map_of_points.end())
          {
              map_of_points[o_vec3d]=++TecPlot_PointIndex;
              m_csh->getNodeCoordsCartesian(xyz, tmp_iv, m_level);                              
              D_TERM(
                x_coords.push_back(xyz[0]);,
                y_coords.push_back(xyz[1]);,					
                z_coords.push_back(xyz[2]));					
          }

          tmp_iv=IntVect(D_DECL(aux_IntVect[0],aux_IntVect[1]+1,aux_IntVect[2]+k));
          o_vec3d.set_data(D_DECL(tmp_iv[0],tmp_iv[1],tmp_iv[2]));
          map_of_points_iter=map_of_points.find(o_vec3d);
          if (map_of_points_iter==map_of_points.end())
          {
              map_of_points[o_vec3d]=++TecPlot_PointIndex;
              m_csh->getNodeCoordsCartesian(xyz, tmp_iv, m_level);                              
              D_TERM(
                x_coords.push_back(xyz[0]);,
                y_coords.push_back(xyz[1]);,					
                z_coords.push_back(xyz[2]));					
          }
        }

        CH_assert(x_coords.size()==TecPlot_PointIndex);
        CH_assert(y_coords.size()==TecPlot_PointIndex);

        const int num_numbers_in_row=30;
        int numbers_in_row=0;	

        std::string zonetype;
        if (CH_SPACEDIM == 2) zonetype = "ZONETYPE=FEQUADRILATERAL";
        if (CH_SPACEDIM == 3) zonetype = "ZONETYPE=FEBRICK";		

        if (num_variables == 1 )
          fprintf(tecplot_file,"ZONE N=%i,E=%i, DATAPACKING=BLOCK, %s, VARLOCATION=([%i]=CELLCENTERED)\n",
            TecPlot_PointIndex,
            (int)(notOverlappedCells.size()),
            zonetype.c_str(),
            CH_SPACEDIM+1);		
        else
          fprintf(tecplot_file,"ZONE N=%i,E=%i, DATAPACKING=BLOCK, %s, VARLOCATION=([%i-%i]=CELLCENTERED)\n",
            TecPlot_PointIndex,
            (int)(notOverlappedCells.size()),
            zonetype.c_str(),
            CH_SPACEDIM+1,
            num_variables+CH_SPACEDIM);		

        for (Real_iter=x_coords.begin();Real_iter!=x_coords.end();++Real_iter)
        {
          value=(*Real_iter);
          fprintf(tecplot_file,"%.6e ",value);
          numbers_in_row++;
          if (numbers_in_row>num_numbers_in_row)
          {
            numbers_in_row=0;
            fprintf(tecplot_file,"\n");
          }			
        }	   			
        if (numbers_in_row!=0) fprintf(tecplot_file,"\n");
        for (Real_iter=y_coords.begin();Real_iter!=y_coords.end();++Real_iter)
        {
          value=(*Real_iter);
          fprintf(tecplot_file,"%.6e ",value);
          numbers_in_row++;
          if (numbers_in_row>num_numbers_in_row)
          {
            numbers_in_row=0;
            fprintf(tecplot_file,"\n");
          }			
        }	
        if (numbers_in_row!=0) fprintf(tecplot_file,"\n");

#if CH_SPACEDIM == 3
        for (Real_iter=z_coords.begin();Real_iter!=z_coords.end();++Real_iter)
        {
          value=(*Real_iter);
          fprintf(tecplot_file,"%.6e ",value);
          numbers_in_row++;
          if (numbers_in_row>num_numbers_in_row)
          {
            numbers_in_row=0;
            fprintf(tecplot_file,"\n");
          }			
        }	
        if (numbers_in_row!=0) fprintf(tecplot_file,"\n");
#endif


        // Writing data
        numbers_in_row=0;			
        for (ivar=0;ivar<num_variables;ivar++)
        {
          i=0;
          for (list_iter=notOverlappedCells.begin();list_iter!=notOverlappedCells.end();++list_iter)
          {
            value=output_data_vector[i*num_variables+ivar];
            fprintf(tecplot_file,"%.6e ",value);				
            numbers_in_row++;
            if (numbers_in_row>num_numbers_in_row)
            {
              numbers_in_row=0;
              fprintf(tecplot_file,"\n");
            }				
            i++;
          }
        }			

        if (numbers_in_row!=0) fprintf(tecplot_file,"\n");

        int ind1,ind2,ind3,ind4;
        // Writing Cells
        for (list_iter=notOverlappedCells.begin();list_iter!=notOverlappedCells.end();++list_iter)
        {
#if CH_SPACEDIM == 3
          for (k=0;k<2;k++)
#endif
          {
            aux_IntVect=(*list_iter);	

            o_vec3d.set_data(D_DECL(aux_IntVect[0],aux_IntVect[1],aux_IntVect[2]+k));
            map_of_points_iter=map_of_points.find(o_vec3d);
            CH_assert(map_of_points_iter!=map_of_points.end());
            ind1=map_of_points_iter->second;

            o_vec3d.set_data(D_DECL(aux_IntVect[0]+1,aux_IntVect[1],aux_IntVect[2]+k));
            map_of_points_iter=map_of_points.find(o_vec3d);
            CH_assert(map_of_points_iter!=map_of_points.end());
            ind2=map_of_points_iter->second;

            o_vec3d.set_data(D_DECL(aux_IntVect[0]+1,aux_IntVect[1]+1,aux_IntVect[2]+k));
            map_of_points_iter=map_of_points.find(o_vec3d);
            CH_assert(map_of_points_iter!=map_of_points.end());
            ind3=map_of_points_iter->second;

            o_vec3d.set_data(D_DECL(aux_IntVect[0],aux_IntVect[1]+1,aux_IntVect[2]+k));
            map_of_points_iter=map_of_points.find(o_vec3d);
            CH_assert(map_of_points_iter!=map_of_points.end());
            ind4=map_of_points_iter->second;

            fprintf(tecplot_file,"%i %i %i %i",ind1,ind2,ind3,ind4);
            #if CH_SPACEDIM == 3
              if (k==0) fprintf(tecplot_file," ");
            #endif
          }
          fprintf(tecplot_file,"\n");
        }
    }// if (notOverlappedCells.size()>0)

	
	
//#endif
}



/// Creates slice
/**
    plane = 'X', x = sValue
    plane = 'Y', y = sValue
    plane = 'Z', z = sValue
 */
 

#ifdef CH_MPI
void AMRLevelIdealMHD::writeSlice(MPI_File tecplot_file, char plane, Real sValue)
#else
void AMRLevelIdealMHD::writeSlice(FILE* tecplot_file, char plane, Real sValue)
#endif
{
#if CH_SPACEDIM == 3
  int i,j; char buffer[1000];

  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();

  int numVars = m_eqSys->numPrimitives() +
                PhPr->numPlotVars() +
                SC->numPlotVars();

  if (m_output_divB == true)
  {
    numVars++;
    Interval comps(UBX,UBZ);
    fillGhostCellsCons(comps);
  }
    
  //bool slice2D = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? true : false);


  if (m_level == 0)
  {
  #ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
    if (procID() == 0)
  #endif
    {
      strcpy(buffer,"VARIABLES=\"X\" \"Y\" \"Z\"");
      
      //if (slice2D == true)
      //{
      //  if (plane == 'X') strcpy(buffer,"VARIABLES=\"Y\" \"Z\"");
      //  if (plane == 'Y') strcpy(buffer,"VARIABLES=\"X\" \"Z\"");
      //  if (plane == 'Z') strcpy(buffer,"VARIABLES=\"X\" \"Y\"");
      //} 

      Vector<string> varNames;
      varNames.reserve(numVars);

      varNames = m_eqSys->primitiveNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");
      }
      varNames = m_patchMHDAM->getPhysProblem()->plotNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");
      }
      varNames = m_patchMHDAM->getPhysProblem()->getSourceCalculator()->plotNames();
      for (i = 0; i < (int)varNames.size(); ++i)
      {
        strcat(buffer," \"");
        strcat(buffer,varNames[i].c_str());
        strcat(buffer,"\"");
      }
      if (m_output_divB == true)
      {
        strcat(buffer," \"divb\"");
      }
      strcat(buffer,"\n");
      std::string auxdata;
      PhPr->auxDataTecplot(auxdata,m_time,0);
      if (auxdata.size()>0)
      {
        strcat(buffer,auxdata.c_str());
        strcat(buffer,"\n");
      }
  #ifdef CH_MPI
      MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
  #else
      fprintf(tecplot_file,"%s",buffer);
  #endif
    }
  #ifdef CH_MPI
    MPI_Barrier(Chombo_MPI::comm);
  #endif
  }

  // Write the data for this level
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  LevelData<FArrayBox> plotData(levelDomain, numVars);

  int iOffset,dirn,dirt1,dirt2, slice_index;
  //FArrayBox localFAB;
  
  if (plane == 'X')  {     dirn = 0; dirt1 = 1; dirt2 = 2; }
  else if (plane == 'Y') { dirn = 1; dirt1 = 0; dirt2 = 2; }
  else if (plane == 'Z') { dirn = 2; dirt1 = 0; dirt2 = 1; }
      
  IntVect iv_off(IntVect::Zero);
  iv_off[dirn]=1;
  
  IntVect iv_t1(IntVect::Zero);
  iv_t1[dirt1]=1;
  
  IntVect iv_t2(IntVect::Zero);
  iv_t2[dirt2]=1;
  
  Real n1,n2;
  
  DataIterator dit;
  DisjointBoxLayout coarsened_fine_domain;  
  
  LevelData<BaseFab<int> > coarsened_fine_tags;
  LevelData<BaseFab<int> >           this_tags;


  if (m_hasFiner)
  if (getFinerLevel()->m_UNew.isDefined())  
  {  
    this_tags.define(m_grids, 1);
    dit = this_tags.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      BaseFab<int>& tags = this_tags[dit()];
      tags.setVal(0);
    }

    coarsen(coarsened_fine_domain, getFinerLevel()->m_grids, m_ref_ratio);
    
    coarsened_fine_tags.define(coarsened_fine_domain,1);    
    dit = coarsened_fine_tags.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {
      BaseFab<int>& tags = coarsened_fine_tags[dit()];
      tags.setVal(1);
    }

    coarsened_fine_tags.copyTo(this_tags);

  }


  dit = plotData.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {
    iOffset = 0;
    const FArrayBox& U       = m_UNew[dit()];
    const Box&       plotBox = levelDomain[dit()];
    BaseFab<int> * tags = (this_tags.isDefined() ? &this_tags[dit()] : NULL);
    
    
    FArrayBox plotFab(plotBox,numVars);

    int LCorner[CH_SPACEDIM]={D_DECL(plotBox.smallEnd()[0], plotBox.smallEnd()[1], plotBox.smallEnd()[2])};
    int UCorner[CH_SPACEDIM]={D_DECL(plotBox.bigEnd()[0],   plotBox.bigEnd()[1],   plotBox.bigEnd()[2])};
    
    Box nCoordsBox(plotBox.smallEnd()*iv_off,plotBox.bigEnd()*iv_off);
    nCoordsBox.surroundingNodes(dirn);

    FArrayBox nCoords(nCoordsBox,1);
    m_csh->getNodeCoords(nCoords, nCoords.box(), dirn, m_level);
    
    Real nmin = nCoords.get(nCoords.smallEnd(),0);
    Real nmax = nCoords.get(nCoords.bigEnd()  ,0);
    
    if ( !((nmin <= sValue) && (sValue < nmax)) ) continue;
    
    slice_index = -1000;
    Box cCoordsBox(plotBox.smallEnd()*iv_off,plotBox.bigEnd()*iv_off);
    BoxIterator bit( cCoordsBox );    
    for( bit.begin(); bit.ok(); ++bit )
    {
      IntVect iv = bit();
      n1 = nCoords.get(iv,0);
      iv[dirn]++;
      n2 = nCoords.get(iv,0);
      if ((n1 <= sValue) && (sValue < n2)) 
      {
        slice_index = iv[dirn]-1;
        break;
      }
    }
    CH_assert(slice_index != -1000);
      
     
    Box Bslice(plotBox);
    Bslice.setSmall(dirn,slice_index);
    Bslice.setBig  (dirn,slice_index);
    
    int nCoveredCells = 0;
  
    if (tags!=NULL)
    {      
      bit.define( Bslice );
      for( bit.begin(); bit.ok(); ++bit )
      {
        IntVect iv = bit();
        nCoveredCells+=(*tags)(iv,0);        
      }      
      if (nCoveredCells == Bslice.numPts())
        continue; // This box is fully covered by finer level, ignore it
        
      #ifndef NDEBUG      
      if (nCoveredCells > 0)
      FORT_VIEWFIABOX(CHF_CONST_FIA((*tags)));     
      #endif      

    }
    

    // Calculation of variables that are necessary for plot
    FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
    m_eqSys->stateToPrim(W, U, plotBox);
    iOffset += m_eqSys->numPrimitives();

    if (PhPr->numPlotVars() > 0)
    {
      //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
      PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
    }
    iOffset += PhPr->numPlotVars();

    if (SC->numPlotVars() > 0)
    {
      CH_assert(SC->numPlotVars() == m_SCData.nComp());
      const FArrayBox& SCData = m_SCData[dit()];
      plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
    }
    
    PhPr->primForPlot(W,plotBox);
    if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);

    if (m_output_divB == true)
    {
      iOffset = plotFab.nComp() - 1;
      FArrayBox divB(plotBox,1);
      (static_cast<PatchMHDMF*>(m_patchMHDAM))->computeDivB(divB, U, plotBox);
      plotFab.copy(divB, plotBox, 0, plotBox, iOffset, 1);
    }

    FArrayBox Fslice(Bslice, numVars);

    // Copy data to Fslice
    Fslice.copy(plotFab);
    

    std::string buffer_main;
    buffer_main.reserve(15*(numVars+CH_SPACEDIM)*plotBox.numPts()+100); // Reserve big enough chunk of memory    

    // It is simple structured mesh
    
    Box nBoxT1(plotBox.smallEnd()*iv_t1,plotBox.bigEnd()*iv_t1);
    Box nBoxT2(plotBox.smallEnd()*iv_t2,plotBox.bigEnd()*iv_t2);    
    nBoxT1.surroundingNodes(dirt1);
    nBoxT2.surroundingNodes(dirt2);
    FArrayBox nCoordsT1(nBoxT1,1);
    FArrayBox nCoordsT2(nBoxT2,1);
    m_csh->getNodeCoords(nCoordsT1, nCoordsT1.box(), dirt1, m_level);    
    m_csh->getNodeCoords(nCoordsT2, nCoordsT2.box(), dirt2, m_level);
    
    
    if (nCoveredCells > 0)
    {
      // Unstructured mesh
      
      int Nsize = (Bslice.size(dirt1)+1)*(Bslice.size(dirt2)+1);      
      int Esize = (Bslice.size(dirt1))*(Bslice.size(dirt2)) - nCoveredCells;      
      
      sprintf(buffer,"ZONE T=L%i N=%i,E=%i, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL, VARLOCATION=([4-%i]=CELLCENTERED)\n", m_level, Nsize, Esize, numVars+3);		
      buffer_main+=buffer;
      
    } else
    {    				  			
      // Structured mesh
      sprintf(buffer,"ZONE T=L%i DATAPACKING=BLOCK, VARLOCATION=([4-%i]=CELLCENTERED), I=%i, J=%i\n", m_level, numVars+3, Bslice.size(dirt1)+1, Bslice.size(dirt2)+1);
      buffer_main+=buffer;
    }
          
    IntVect iv; RealVect cv_coord,cs_coord;
    iv[dirn] = slice_index;
  
    for (j=LCorner[dirt2];j<=UCorner[dirt2]+1;j++)
    {
      for (i=LCorner[dirt1];i<=UCorner[dirt1]+1;i++)
      {
        cv_coord[dirt1] = nCoordsT1.get(i*iv_t1,0);
        cv_coord[dirt2] = nCoordsT2.get(j*iv_t2,0);
        cv_coord[dirn]  = sValue;
        
        m_csh->transCurvCoordsToCartesian(cs_coord, cv_coord);
                
        sprintf(buffer,"%.6e ",cs_coord[0]);
        buffer_main+=buffer;
      }
      buffer_main+="\n";
    }				
    for (j=LCorner[dirt2];j<=UCorner[dirt2]+1;j++)
    {
      for (i=LCorner[dirt1];i<=UCorner[dirt1]+1;i++)      
      {
        cv_coord[dirt1] = nCoordsT1.get(i*iv_t1,0);
        cv_coord[dirt2] = nCoordsT2.get(j*iv_t2,0);
        cv_coord[dirn]  = sValue;
        
        m_csh->transCurvCoordsToCartesian(cs_coord, cv_coord);
        
        sprintf(buffer,"%.6e ",cs_coord[1]);        
        buffer_main+=buffer;
      }
      buffer_main+="\n";
    }
    for (j=LCorner[dirt2];j<=UCorner[dirt2]+1;j++)
    {
      for (i=LCorner[dirt1];i<=UCorner[dirt1]+1;i++)      
      {
        cv_coord[dirt1] = nCoordsT1.get(i*iv_t1,0);
        cv_coord[dirt2] = nCoordsT2.get(j*iv_t2,0);
        cv_coord[dirn]  = sValue;
        
        m_csh->transCurvCoordsToCartesian(cs_coord, cv_coord);
        
        sprintf(buffer,"%.6e ",cs_coord[2]);        
        buffer_main+=buffer;
      }
      buffer_main+="\n";
    }    
        
    int numbers_in_row = 0,num_numbers_in_row = UCorner[dirt1] - LCorner[dirt1] + 1;
    
    Real value; IntVect aux_IntVect;int ivar;
    for (ivar=0;ivar<numVars;ivar++)
    {		                
      for (j=LCorner[dirt2];j<=UCorner[dirt2];j++)
      {
        for (i=LCorner[dirt1];i<=UCorner[dirt1];i++)      
        {
          aux_IntVect[dirt1]=i;aux_IntVect[dirt2]=j;aux_IntVect[dirn]=slice_index;		
           
          if (tags!=NULL)
          {
            int tg = (*tags)(aux_IntVect,0);
            if (tg == 1)  continue; // This cell is fully covered by finer level, ignore it              
          }
                                               
          value=Fslice.get(aux_IntVect,ivar);
          sprintf(buffer,"%.6e ",value);				
          buffer_main+=buffer;
          numbers_in_row++;
        }			        
        
        if (numbers_in_row>=num_numbers_in_row)
        {
          numbers_in_row=0;
          buffer_main+="\n";
        }        
      }      
    }    
    
    //if (nCoveredCells > 0)
    //{      
    //  pout() << check << " q " << nCoveredCells << endl;      
    //}
    
    if (nCoveredCells > 0)
    {
      // Write connectivity list
      CH_assert(tags!=NULL);
      
      int nx =  UCorner[dirt1] - LCorner[dirt1] + 2;
      int ny =  UCorner[dirt2] - LCorner[dirt2] + 2;
      
      ny = 1;
          
      for (j=LCorner[dirt2];j<=UCorner[dirt2];j++)      
      for (i=LCorner[dirt1];i<=UCorner[dirt1];i++)      
      {
        aux_IntVect[dirt1]=i;aux_IntVect[dirt2]=j;aux_IntVect[dirn]=slice_index;		
                 
        int tg = (*tags)(aux_IntVect,0);
        if (tg == 1) continue; // This cell is fully covered by finer level, ignore it
        
        int n1 = 1+(i  -LCorner[dirt1])*ny + nx*(j-LCorner[dirt2]);
        int n2 = 1+(i+1-LCorner[dirt1])*ny + nx*(j-LCorner[dirt2]);
        int n3 = 1+(i+1-LCorner[dirt1])*ny + nx*(j-LCorner[dirt2]+1);
        int n4 = 1+(i  -LCorner[dirt1])*ny + nx*(j-LCorner[dirt2]+1);
        
        sprintf(buffer,"%i %i %i %i\n",n1,n2,n3,n4);				
        buffer_main+=buffer;                
      }                          
    }
        
    

#ifdef CH_MPI
    char *buffer_to_write = new char[buffer_main.size()+1];
    buffer_main.copy(buffer_to_write, std::string::npos);
    buffer_to_write[buffer_main.size()] = 0;
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(tecplot_file,"%s",buffer_main.c_str());
    fprintf(tecplot_file,"# Level %i ended\n",m_level);
#endif
  }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
  if (procID() == 0)
  {
    sprintf(buffer,"# Level %i ended\n",m_level);				
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
  }  
  MPI_Barrier(Chombo_MPI::comm);
#endif

#endif  // #if CH_SPACEDIM == 3
}



int AMRLevelIdealMHD::writeOneLine(FILE * tecplot_file, Box a_box, int a_dirLine, LevelData<FArrayBox>* a_lineData)
{  
  
  AMRLevelIdealMHD* amrMHDFinerPtr   = getFinerLevel();
  
  if (amrMHDFinerPtr!=NULL)
  if (amrMHDFinerPtr->m_grids.size() == 0) amrMHDFinerPtr = NULL;
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();

  int numVars = m_eqSys->numPrimitives() +
                PhPr->numPlotVars() +
                SC->numPlotVars();
  
  DataIterator dit = a_lineData[m_level].dataIterator();
  
  int returnValue = -1;
  
  IntVect iv; int iOffset; int i;
  
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& U = a_lineData[m_level][dit()];
    Box b = a_box & U.box();
    if (b.isEmpty()) continue;
    
    returnValue = 1;
    
    BoxIterator bit( b );    
    for( bit.begin(); bit.ok(); ++bit )
    {
      iv = bit();
      Box b_fine(iv,iv);
      b_fine.refine(refRatio());
      
      bool overlapped =  false;
      if (amrMHDFinerPtr!=NULL)
        overlapped = (amrMHDFinerPtr->writeOneLine(tecplot_file,b_fine,a_dirLine,a_lineData) == 1);
      
      if (overlapped) continue;
      
      
      
      Box plotBox(iv,iv);        
      FArrayBox plotFab(plotBox,numVars);
  
      iOffset = 0;
      // Calculation of variables that are necessary for plot
      FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
      m_eqSys->stateToPrim(W, U, plotBox);
      iOffset += m_eqSys->numPrimitives();

      if (PhPr->numPlotVars() > 0)
      {
        //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
        PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
      }
      iOffset += PhPr->numPlotVars();

      if (SC->numPlotVars() > 0)
      {
        CH_assert(SC->numPlotVars() == m_SCData.nComp());
        const FArrayBox& SCData = m_SCData[dit()];
        plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
      }
  
      PhPr->primForPlot(W,plotBox);
      if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);
      
      RealVect cv_coords,cs_coords;
      
      m_csh->getCellCenter( cv_coords, iv, m_level);
      //fprintf(tecplot_file,"level %i (%i,%i) ",m_level, iv[0],iv[1]);                
      fprintf(tecplot_file,"%.6e ",cv_coords[a_dirLine]);                
      m_csh->transCurvCoordsToCartesian(cs_coords, cv_coords);
      fprintf(tecplot_file,"%.6e %.6e ",cs_coords[0],cs_coords[1]);                
      if (SpaceDim == 3) fprintf(tecplot_file,"%.6e ",cs_coords[2]);
      
      for (i=0;i<numVars;i++)
      {
        Real value = plotFab.get(iv,i);
        fprintf(tecplot_file,"%.6e ", value);
      }
      fprintf(tecplot_file,"\n");                                   
    }      
          
  }  
  
  return returnValue;
  
}


int AMRLevelIdealMHD::writeProbeData(const ProbeFilesInfo::ProbeInfo & a_probe)
{
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();

  int numVars = m_eqSys->numPrimitives() + 
                PhPr->numPlotVars() +
                SC->numPlotVars();
  
  int i;  char buffer[1000],buf1[40];  
    
  Real time = PhPr->getPhysTime(m_time);
      
                
  RealVect pos;
  a_probe.get1Dposition(pos, time);
                    
  IntVect ivAU;
  m_csh->getClosestCell(ivAU, pos, m_level); 
  
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
                 
  int local_found = 0;
    
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {  
    const FArrayBox& U = m_UNew[dit()];
    const Box&       b = levelDomain[dit()];
    
    if (b.contains(ivAU)==false) continue;
    
    Box plotBox(ivAU,ivAU);
    plotBox.grow(1);
    FArrayBox plotFab(plotBox,numVars);
    
    int iOffset = 0;
    // Calculation of variables that are necessary for plot
    FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
    m_eqSys->stateToPrim(W, U, plotBox);
    iOffset += m_eqSys->numPrimitives();

    if (PhPr->numPlotVars() > 0)
    {
      //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
      PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
    }
    iOffset += PhPr->numPlotVars();

    if (SC->numPlotVars() > 0)
    {
      CH_assert(SC->numPlotVars() == m_SCData.nComp());
      const FArrayBox& SCData = m_SCData[dit()];
      plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
    }
    
    PhPr->primForPlot(W,plotBox);
    if (m_output_vecCS == false)
    m_csh->transCartesianVectToCurv(W,plotBox,m_level);
    
    Real lon = 180.0*pos[1]/d_PI;
    Real lat = 180.0*(d_PI_2-pos[2])/d_PI;
    
    sprintf(buffer,"%.12e %.12e %.12e %.12e",time, pos[0], lon, lat);
    //for (i=0;i<numVars;i++) 
    //{
    //  Real v = W.get(ivAU,i);
    //  sprintf(buf1," %.12e ",v);
    //  strcat(buffer,buf1);
    //}      
    
    RealVect rvAU;
    m_csh->getCellCenter(rvAU, ivAU, m_level);
    Real slopes[SpaceDim]; 
    Real v,v1,slope;
    Real * Wrv = new Real[numVars];
    
    for (i=0;i<numVars;i++) 
    {
      Real v = W.get(ivAU,i);
      Wrv[i] = v;
    }      
    
    for (int dir = 0; dir < SpaceDim; ++dir)
    {
      IntVect ivAU1(ivAU); RealVect rvAU1;
      if (pos[dir]<rvAU[dir])
      {        
        ivAU1[dir]-=1;
        m_csh->getCellCenter(rvAU1, ivAU1, m_level);
        for (i=0;i<numVars;i++) 
        {
          v1 = W.get(ivAU1,i);
          v  = W.get(ivAU,i);
          slope = (v-v1)/(rvAU[dir]-rvAU1[dir]);      
          Wrv[i] += slope*(pos[dir] - rvAU[dir]);
        }
      } else
      {        
        ivAU1[dir]+=1;
        m_csh->getCellCenter(rvAU1, ivAU1, m_level);
        for (i=0;i<numVars;i++) 
        {
          v1 = W.get(ivAU1,i);
          v  = W.get(ivAU,i);
          slope = (v1-v)/(rvAU1[dir]-rvAU[dir]);          
          Wrv[i]+= slope*(pos[dir] - rvAU[dir]);
        }
      }        
    }      
               
    //interpolateData(rvAU,ivAU,W,Wrv);  
                
    for (i=0;i<numVars;i++) 
    {                
      v  = Wrv[i];
      
      sprintf(buf1," %.12e",v);
      strcat(buffer,buf1);
    }
    delete[] Wrv;
    
    strcat(buffer,"\n");
        
    
    a_probe.writeString(buffer);
    
//#ifdef CH_MPI          
//    MPI_File_write_shared(a_pfi.m_files[iProbe], buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
//#else      
//    fprintf(a_pfi.m_files[iProbe],"%s",buffer);      
//#endif      
    local_found = 1;
        
    break;
        
  }
  
  int found;
  MPI_Allreduce(&local_found, &found, 1, MPI_INT, MPI_SUM, Chombo_MPI::comm);
  
  CH_assert((found == 0) || (found == 1));
  
  return found;
    
    //if (found == false)
    //{
      //char err_str[100];
      //if (SpaceDim == 2) sprintf(err_str,"probe point p=(%g,%g) ij=(%i,%i) was not found in mesh",pos[0],pos[1],ivAU[0],ivAU[1]);
      //if (SpaceDim == 3) sprintf(err_str,"probe point p=(%g,%g,%g) ijk=(%i,%i,%i) was not found in mesh",pos[0],pos[1],pos[2],ivAU[0],ivAU[1],ivAU[2]);      
      //MayDay::Error(err_str);
    //}
  
  
}

// Not used any more
void AMRLevelIdealMHD::write1DProbe(ProbeFilesInfo& a_pfi)
{
/*
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();
  SourceCalculator* SC = PhPr->getSourceCalculator();

  int numVars = m_eqSys->numPrimitives() + 
              PhPr->numPlotVars() +
              SC->numPlotVars();
  
  int i;            
  Vector<RealVect> probePoints;  
  
  PhPr->probe1DPoints(probePoints,m_time);
  
  CH_assert(probePoints.size() == a_pfi.m_files.size());
  
  m_UNew.exchange();  
  
  char buffer[1000],buf1[40];  
  
  //FArrayBox localFAB;  
  Real time = PhPr->getPhysTime(m_time);
  
  for (int iProbe = 0;iProbe < probePoints.size();++iProbe)
  {    
      
    RealVect pAU = probePoints[iProbe];    
    
    IntVect ivAU;
    
    m_csh->getClosestCell(ivAU, pAU, m_level);    
            
    
    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
    
    bool found = false;
    
    DataIterator dit = m_UNew.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {  
      const FArrayBox& U = m_UNew[dit()];
      const Box&       b = levelDomain[dit()];
      
      if (b.contains(ivAU)==false) continue;
      
      Box plotBox(ivAU,ivAU);
      plotBox.grow(1);
      FArrayBox plotFab(plotBox,numVars);
      
      int iOffset = 0;
      // Calculation of variables that are necessary for plot
      FArrayBox  W(plotBox, numVars, plotFab.dataPtr(iOffset));
      m_eqSys->stateToPrim(W, U, plotBox);
      iOffset += m_eqSys->numPrimitives();

      if (PhPr->numPlotVars() > 0)
      {
        //localFAB.define(plotBox, numVars - m_eqSys->numPrimitives(), plotFab.dataPtr(iOffset));
        PhPr->calcPlotVars(plotFab, m_eqSys->numPrimitives(), W, plotBox);
      }
      iOffset += PhPr->numPlotVars();

      if (SC->numPlotVars() > 0)
      {
        CH_assert(SC->numPlotVars() == m_SCData.nComp());
        const FArrayBox& SCData = m_SCData[dit()];
        plotFab.copy(SCData, plotBox, 0, plotBox, iOffset, m_SCData.nComp());
      }
      
      PhPr->primForPlot(W,plotBox);
      if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);
      
      Real lon = 180.0*pAU[1]/d_PI;
      Real lat = 180.0*(d_PI_2-pAU[2])/d_PI;
      
      sprintf(buffer,"%.12e %.12e %.12e %.12e",time, pAU[0], lon, lat);
      //for (i=0;i<numVars;i++) 
      //{
      //  Real v = W.get(ivAU,i);
      //  sprintf(buf1," %.12e ",v);
      //  strcat(buffer,buf1);
      //}      
      
      RealVect rvAU;
      m_csh->getCellCenter(rvAU, ivAU, m_level);
      Real slopes[SpaceDim]; 
      Real v,v1,slope;
      Real * Wrv = new Real[numVars];
      
      for (i=0;i<numVars;i++) 
      {
        Real v = W.get(ivAU,i);
        Wrv[i] = v;
      }      
      
      for (int dir = 0; dir < SpaceDim; ++dir)
      {
        IntVect ivAU1(ivAU); RealVect rvAU1;
        if (pAU[dir]<rvAU[dir])
        {        
          ivAU1[dir]-=1;
          m_csh->getCellCenter(rvAU1, ivAU1, m_level);
          for (i=0;i<numVars;i++) 
          {
            v1 = W.get(ivAU1,i);
            v  = W.get(ivAU,i);
            slope = (v-v1)/(rvAU[dir]-rvAU1[dir]);      
            Wrv[i] += slope*(pAU[dir] - rvAU[dir]);
          }
        } else
        {        
          ivAU1[dir]+=1;
          m_csh->getCellCenter(rvAU1, ivAU1, m_level);
          for (i=0;i<numVars;i++) 
          {
            v1 = W.get(ivAU1,i);
            v  = W.get(ivAU,i);
            slope = (v1-v)/(rvAU1[dir]-rvAU[dir]);          
            Wrv[i]+= slope*(pAU[dir] - rvAU[dir]);
          }
        }        
      }      
                 
      //interpolateData(rvAU,ivAU,W,Wrv);  
                  
      for (i=0;i<numVars;i++) 
      {                
        v  = Wrv[i];
        
        sprintf(buf1," %.12e",v);
        strcat(buffer,buf1);
      }
      delete[] Wrv;
      
      strcat(buffer,"\n");
      
      //FILE* plot_file=fopen(filenames[iProbe].c_str(),"a");
      //fprintf(plot_file,"%s\n",buffer);
      //fclose(plot_file);
#ifdef CH_MPI          
      MPI_File_write_shared(a_pfi.m_files[iProbe], buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else      
      fprintf(a_pfi.m_files[iProbe],"%s",buffer);      
#endif      
      found = true;
          
      break;
          
    }
    
    //if (found == false)
    //{
      //char err_str[100];
      //if (SpaceDim == 2) sprintf(err_str,"probe point p=(%g,%g) ij=(%i,%i) was not found in mesh",pAU[0],pAU[1],ivAU[0],ivAU[1]);
      //if (SpaceDim == 3) sprintf(err_str,"probe point p=(%g,%g,%g) ijk=(%i,%i,%i) was not found in mesh",pAU[0],pAU[1],pAU[2],ivAU[0],ivAU[1],ivAU[2]);      
      //MayDay::Error(err_str);
    //}
  }  
  */      
}

void AMRLevelIdealMHD::writeAxisXData(int a_step, Real a_time)
{


  int i,li;

  Real* point_data;

  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  DataIterator it = levelDomain.dataIterator();

  //int num_layers = 2;
  int num_layers = m_problem_domain.domainBox().size(1);

  std::map<int,Real*> * sorted_data = new std::map<int,Real*>[num_layers];

  IntVect aux_IntVect;

  for( it.begin(); it.ok(); ++it)
  {
      const FArrayBox& data = m_UNew[it()];
      const Box& box = levelDomain[it()];

      IntVect LCorner=box.smallEnd();
      IntVect UCorner=box.bigEnd();

      if (LCorner[1]>(num_layers-1)) continue;

      //for (li=0;li<num_layers;li++)
      for (li=LCorner[1];li<=UCorner[1];li++)
      for (i=LCorner[0];i<=UCorner[0];i++)
      {
        aux_IntVect[0]=i;aux_IntVect[1]=li;		
        if (sorted_data[li].find(i)!=sorted_data[li].end()) continue;

        point_data=new Real[6];
        point_data[0]=data.get(aux_IntVect,URHO);
        point_data[1]=data.get(aux_IntVect,UMOMX)/point_data[0];
        point_data[2]=data.get(aux_IntVect,UMOMY)/point_data[0];
        point_data[3]=data.get(aux_IntVect,UENG);
        point_data[4]=data.get(aux_IntVect,UBX);
        point_data[5]=data.get(aux_IntVect,UBY);

        sorted_data[li][i]=point_data;
      }
  }

  for (li=0;li<num_layers;li++)
  CH_assert(sorted_data[li].size()==m_problem_domain.domainBox().size(0));


  char file_name[1000];
  sprintf(file_name,"X%06i.dat",a_step);
  FILE* plot_file=fopen(file_name,"w");
  fprintf(plot_file,"T= %.6E\n",a_time);
  fprintf(plot_file,"Density W U E BZ BX\n");

  std::map<int,Real*>::iterator map_iter;

  for (li=0;li<num_layers;li++)
  {
    fprintf(plot_file,"row: %i\n",li);
    for (map_iter=sorted_data[li].begin();map_iter!=sorted_data[li].end();++map_iter)
    {
      point_data=map_iter->second;
      fprintf(plot_file,"%i %.15E %.15E %.15E %.15E %.15E %.15E\n",map_iter->first,point_data[0],point_data[1],
         point_data[2],point_data[3],point_data[4],point_data[5]);
      delete[] point_data;
    }
    fprintf(plot_file," \n");
  }
  fclose(plot_file);

  delete[] sorted_data;

}

void AMRLevelIdealMHD::loadfort55()
{
  int nx = m_problem_domain.domainBox().size(0);
  int ny = m_problem_domain.domainBox().size(1);

  //int nx55 = nx+4;
  //int ny55 = ny+4;

  IntVect LCorner=m_problem_domain.domainBox().smallEnd();
  IntVect UCorner=m_problem_domain.domainBox().bigEnd();

  Box WBox = m_problem_domain.domainBox();
  WBox.grow(0,2);
  WBox.grow(1,2);

  FArrayBox W(WBox,m_eqSys->numStates());
  FArrayBox U(m_problem_domain.domainBox(),m_eqSys->numStates());


  W.setVal(0.0);
  U.setVal(0.0);

  FORT_READFORT55(
     CHF_FRA(W));

  /*IntVect aux_IntVect;

  int read_items,i,ind,pos,indy;
  double time;

  double* r = new double[ny55];
  double* pg = new double[ny55];
  double* u = new double[ny55];
  double* w = new double[ny55];
  double* bx = new double[ny55];
  double* bz = new double[ny55];

  FILE* cnp_file = fopen("fort.55","rb");

  read_items=fread(&time, sizeof(double),1,cnp_file);	

  for(i=0;i<nx55;i++)
  {
		read_items=fread(r, sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);
    read_items=fread(pg,sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);
    read_items=fread(u, sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);
    read_items=fread(w, sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);
    read_items=fread(bx,sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);
    read_items=fread(bz,sizeof(double),ny55,cnp_file);	
    CH_assert(read_items==ny55);

    if ((i<2) || (i>(nx55 - 2))) continue;

    aux_IntVect[0]= LCorner[0]+i-2;

    for (indy=LCorner[1];indy<=UCorner[1];indy++)
    {
      aux_IntVect[1]=indy;	
      pos=indy-LCorner[1]+2;
      W.set(aux_IntVect,WRHO,  r[pos]);
      W.set(aux_IntVect,WVELX, w[pos]);
      W.set(aux_IntVect,WVELY, u[pos]);
      W.set(aux_IntVect,WPRES, pg[pos]);
      W.set(aux_IntVect,WBX,   bz[pos]);
      W.set(aux_IntVect,WBY,   bx[pos]);
    }
		
  }
  fclose(cnp_file);

  delete[] r;
  delete[] pg;
  delete[] u;
  delete[] w;
  delete[] bx;
  delete[] bz;*/

  FORT_PRIMTOCONS(CHF_FRA(U),
                  CHF_CONST_FRA(W),
                  CHF_BOX(m_problem_domain.domainBox()));

  //m_patchMHDAM->primToCons(U,W,m_problem_domain.domainBox());


  DataIterator dit = m_UNew.boxLayout().dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {

    FArrayBox& Udata = m_UNew[dit()];

    // Box of current grid
    Box uBox = Udata.box();
    uBox &= m_problem_domain.domainBox();

    Udata.copy(U, uBox);
  }


}

// write level 0 in text format
void AMRLevelIdealMHD::writeLevel0SF(int a_step)
{
  D_TERM(
    int nx = m_problem_domain.domainBox().size(0);,
    int ny = m_problem_domain.domainBox().size(1);,
    int nz = m_problem_domain.domainBox().size(2););
  
  
  int D_DECL(i,j,k);
  
  Real dx = (m_csh->constStep(0) ? m_csh->dx(0,m_level) : 1.0);

  IntVect LCorner=m_problem_domain.domainBox().smallEnd();
  IntVect UCorner=m_problem_domain.domainBox().bigEnd();
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();    
  

  
  char file_name[100],buf[400];
  sprintf(file_name,"LZ%06i.dat",a_step);
  
#ifdef CH_MPI
  MPI_File data_file;
  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &data_file);
  MPI_File_close(&data_file);

  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &data_file);
  MPI_File_set_view(data_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
#else
  FILE*     data_file=fopen(file_name,"w");
#endif

  std::string big_buf;
  big_buf.reserve(15*(D_TERM(UCorner[0],+UCorner[1],+UCorner[2])));
  

  if (procID() == 0)
  {    
#if CH_SPACEDIM==2
    sprintf(buf,"%i %in",nx,ny);                  
#else
    sprintf(buf,"%i %i %i\n",nx,ny,nz);
#endif      
    big_buf+=buf;

    for (i=LCorner[0];i<=UCorner[0];i++){sprintf(buf,"%g ",(i+0.5)*dx);big_buf.append(buf);}
    big_buf.append("\n");
    for (j=LCorner[1];j<=UCorner[1];j++){sprintf(buf,"%g ",(j+0.5)*dx);big_buf.append(buf);}
    big_buf.append("\n");
#if CH_SPACEDIM==3    
    for (k=LCorner[2];k<=UCorner[2];k++){sprintf(buf,"%g ",(k+0.5)*dx);big_buf.append(buf);}
    big_buf.append("\n");
#endif    
  }
  
#ifdef CH_MPI
    char *buffer_to_write = new char[big_buf.size()+1];
    big_buf.copy(buffer_to_write, std::string::npos);
    MPI_File_write_shared(data_file, buffer_to_write, big_buf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(data_file,"%s",big_buf.c_str());
#endif   

  int numVars = m_eqSys->numPrimitives() + 
                  PhPr->numPlotVars();
                
      
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
  
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {  
    const FArrayBox& U       = m_UNew[dit()];
    const Box&       plotBox = levelDomain[dit()];
    
    FArrayBox W       (plotBox,m_eqSys->numPrimitives());
    FArrayBox plotVars(plotBox,PhPr   ->numPlotVars());
        
    m_eqSys->stateToPrim(W, U, plotBox);    
    PhPr->calcPlotVars(plotVars, 0, W, plotBox);           
    PhPr->primForPlot(W,plotBox);    
    
    
    big_buf.reserve(15*(numVars+CH_SPACEDIM)*plotBox.numPts()+100); // Reserve big enough chunk of memory            
              
    BoxIterator bit(plotBox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit(); 

#if CH_SPACEDIM==2
      sprintf(buf,"%i %i ",iv[0],iv[1]);                  
#else
      sprintf(buf,"%i %i %i ",iv[0],iv[1],iv[2]);
#endif      
      big_buf+=buf;
      
      sprintf(buf,"%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n",
        W.get(iv,WRHO),W.get(iv,WVELX),W.get(iv,WVELY),W.get(iv,WVELZ),W.get(iv,WPRES),W.get(iv,WBX),W.get(iv,WBY),W.get(iv,WBZ));
        
      big_buf+=buf;
    }      
       
#ifdef CH_MPI
    char *buffer_to_write = new char[big_buf.size()+1];
    big_buf.copy(buffer_to_write, std::string::npos);
    MPI_File_write_shared(data_file, buffer_to_write, big_buf.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(data_file,"%s",buffer_main.c_str());
#endif  
  }


#ifdef CH_MPI
  MPI_File_close(&data_file);
#else
  fclose(data_file);
#endif

  
}
 
void AMRLevelIdealMHD::MedvedevOutput(int a_step) const
{
  if (m_level!=1) return;
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  
  int iWRHO1 = eqSys->densityIndexPrim(1);
  int iWRHO2 = eqSys->densityIndexPrim(2);
  int iWRHO3 = eqSys->densityIndexPrim(3);  
  

  char file_name[100],buffer[1000];
  sprintf(file_name,"M%06i.dat",a_step);

#ifdef CH_MPI
  MPI_File tecplot_file;
  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &tecplot_file);
  MPI_File_close(&tecplot_file);

  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &tecplot_file);
  MPI_File_set_view(tecplot_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
#else
  FILE*     tecplot_file=fopen(file_name,"w");
#endif
 
  int imax = 864; // ~600AU

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
  if (procID() == 0)
#endif
  {
    int nx = m_problem_domain.size(0);
    int ny = m_problem_domain.size(1);
    
    nx = imax+1;
    Real time = PhPr->getPhysTime(m_time);
    
    sprintf(buffer,"%.12e\n%i %i\n",time,nx,ny);
          
#ifdef CH_MPI
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else
    fprintf(tecplot_file,"%s",buffer);
#endif
  }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif  


  
    
  int numVars = m_eqSys->numPrimitives() + 
                  PhPr->numPlotVars();
                
      
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
  
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {  
    const FArrayBox& U       = m_UNew[dit()];
    const Box&       plotBox = levelDomain[dit()];
    
    FArrayBox W       (plotBox,m_eqSys->numPrimitives());
    FArrayBox plotVars(plotBox,PhPr   ->numPlotVars());
        
    m_eqSys->stateToPrim(W, U, plotBox);    
    PhPr->calcPlotVars(plotVars, 0, W, plotBox);           
    PhPr->primForPlot(W,plotBox);    
    
    std::string buffer_main;
    buffer_main.reserve(15*(numVars+CH_SPACEDIM)*plotBox.numPts()+100); // Reserve big enough chunk of memory            
          
    RealVect cv_coord;
    BoxIterator bit(plotBox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit(); 
      if (iv[0]>imax) continue;
      m_csh->getCellCenter(cv_coord, iv, m_level);
      Real nH = W.get(iv,iWRHO1)+W.get(iv,iWRHO2)+W.get(iv,iWRHO3);
      Real vx = W.get(iv,WVELX)*1e+5;
      Real vy = W.get(iv,WVELY)*1e+5;
      sprintf(buffer,"%i %i %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",iv[0],iv[1],cv_coord[0],cv_coord[1],W.get(iv,WRHO),nH,vx,vy,plotVars.get(iv,0));
      //sprintf(buffer,"%i %i\n",iv[0],iv[1]);                        
      buffer_main+=buffer;
    }      
       
#ifdef CH_MPI
    char *buffer_to_write = new char[buffer_main.size()+1];
    buffer_main.copy(buffer_to_write, std::string::npos);
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(tecplot_file,"%s",buffer_main.c_str());
#endif  
  }


#ifdef CH_MPI
  MPI_File_close(&tecplot_file);
#else
  fclose(tecplot_file);
#endif

}

void AMRLevelIdealMHD::V1V2TSOutput(int a_step) const
{
  if (m_level!=1) return;
  
  
  Real V2TS[3]={86,0.6302118173999999850565246,2.207726999000000134287802};
  Real dangle = 8.0*d_PI_180;
  Real dr     = 6.0;
  
  RealVect rLo(D_DECL(V2TS[0]-dr,V2TS[1]-dangle,V2TS[2]-dangle));
  RealVect rHi(D_DECL(V2TS[0]+dr,V2TS[1]+dangle,V2TS[2]+dangle));
  
  IntVect ivLo;IntVect ivHi;
  
  m_csh->getClosestCell(ivLo,rLo,m_level);
  m_csh->getClosestCell(ivHi,rHi,m_level);
  
  Box b(ivLo,ivHi);
  
  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
  

  char file_name[100],buffer[1000];
  sprintf(file_name,"V2TS%06i.dat",a_step);

#ifdef CH_MPI
  MPI_File tecplot_file;
  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &tecplot_file);
  MPI_File_close(&tecplot_file);

  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &tecplot_file);
  MPI_File_set_view(tecplot_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
#else
  FILE*     tecplot_file=fopen(file_name,"w");
#endif
   

#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
  if (procID() == 0)
#endif
  {    
    Real time = PhPr->getPhysTime(m_time);
    
    sprintf(buffer,"%.12e\nbox (%i %i %i) - (%i %i %i)\ni j k R phi theta density VR VPhi VTheta BR BPhi BTheta T\n",
            time,ivLo[0],ivLo[1],ivLo[2],ivHi[0],ivHi[1],ivHi[2]);
          
#ifdef CH_MPI
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else
    fprintf(tecplot_file,"%s",buffer);
#endif
  }
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif  


  
    
  int numVars = m_eqSys->numPrimitives() + 
                  PhPr->numPlotVars();
                
      
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
  
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {  
    const FArrayBox& U       = m_UNew[dit()];
    const Box&       plotBox = levelDomain[dit()];
    
    if (b.intersects(plotBox)==false) continue;
    
    FArrayBox W       (plotBox,m_eqSys->numPrimitives());
    FArrayBox plotVars(plotBox,PhPr   ->numPlotVars());
        
    m_eqSys->stateToPrim(W, U, plotBox);    
    PhPr->calcPlotVars(plotVars, 0, W, plotBox);           
    PhPr->primForPlot(W,plotBox);    
    
    std::string buffer_main;
    buffer_main.reserve(15*(numVars+CH_SPACEDIM)*plotBox.numPts()+100); // Reserve big enough chunk of memory            
          
    RealVect cv_coord;
    BoxIterator bit(plotBox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit(); 
      if (b.contains(iv) == false) continue;      
      
      m_csh->getCellCenter(cv_coord, iv, m_level);      
      
      sprintf(buffer, "%i %i %i %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e\n",
              iv[0],iv[1],iv[2],cv_coord[0],cv_coord[1],cv_coord[2],             
              W.get(iv,WRHO),W.get(iv,WVELR),W.get(iv,WVELP),W.get(iv,WVELT),W.get(iv,WBR),W.get(iv,WBP),W.get(iv,WBT),
              plotVars.get(iv,0));
      //sprintf(buffer,"%i %i\n",iv[0],iv[1]);                        
      buffer_main+=buffer;
    }      
       
#ifdef CH_MPI
    char *buffer_to_write = new char[buffer_main.size()+1];
    buffer_main.copy(buffer_to_write, std::string::npos);
    MPI_File_write_shared(tecplot_file, buffer_to_write, buffer_main.size(), MPI_CHAR, MPI_STATUS_IGNORE);
    delete[] buffer_to_write;
#else
    fprintf(tecplot_file,"%s",buffer_main.c_str());
#endif  
  }


#ifdef CH_MPI
  MPI_File_close(&tecplot_file);
#else
  fclose(tecplot_file);
#endif

}


void AMRLevelIdealMHD::wrtieLevel0H5() const
{
  if (m_level!=0) return;
      
  hid_t dataspace,attr,dataset,memspace;
  herr_t status;        
      
  hid_t plist_file    = H5P_DEFAULT;
  hid_t plist_dataset = H5P_DEFAULT;
  
#ifdef CH_MPI        
  // Set up file access property list with parallel I/O access     
  MPI_Info info = MPI_INFO_NULL;
  plist_file = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plist_file, Chombo_MPI::comm, info);
  
  // Create property list for independent dataset write.  
  //plist_dataset = H5Pcreate(H5P_DATASET_XFER);
  //H5Pset_dxpl_mpio(plist_dataset, H5FD_MPIO_INDEPENDENT);  
  //H5Pset_dxpl_mpio(plist_dataset,H5FD_MPIO_COLLECTIVE);
#endif

  
  PhysProblem* PhPr    = m_patchMHDAM->getPhysProblem();    
  
  int numVars = m_eqSys->numPrimitives() + 
                PhPr->numPlotVars();
  
  const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();
  LevelData<FArrayBox> plotData(levelDomain, numVars);
  
  DataIterator dit = m_UNew.dataIterator();
  for (dit.begin(); dit.ok(); ++dit)
  {  
    const FArrayBox& U = m_UNew[dit()];
          FArrayBox& W = plotData[dit()];
    const Box& plotBox = levelDomain[dit()];      
                      
        
    m_eqSys->stateToPrim(W, U, plotBox);          
    PhPr->calcPlotVars(W, m_eqSys->numPrimitives(), W, plotBox);
    PhPr->primForPlot(W,plotBox);        
    if (m_output_vecCS == false)
      m_csh->transCartesianVectToCurv(W,plotBox,m_level);
  }
  
  hid_t * h5files = new hid_t[numVars];
  
  Vector<string> varNames;
  varNames.reserve(numVars);
  
  int D_DECL(i,j,k); int ivar,dim;

  varNames = m_eqSys->primitiveNames();  
  for (i = 0; i < varNames.size(); ++i)
  {
    varNames[i]+=".h5";    
    h5files[i] = H5Fcreate(varNames[i].c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file);
  }  
  int shift = varNames.size();
  varNames = PhPr->plotNames();
  for (i = 0; i < varNames.size(); ++i)
  {
    varNames[i]+=".h5";    
    h5files[i+shift] = H5Fcreate(varNames[i].c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_file);
  }
  
  hsize_t dimsf[SpaceDim];
  hsize_t offset[SpaceDim];
  hsize_t count[SpaceDim];
  for (dim = 0; dim < SpaceDim; ++dim) dimsf[dim] = m_problem_domain.domainBox().size(dim);
  
  
  for (ivar = 0; ivar < numVars; ++ivar)
  //for (ivar = numVars - 1; ivar < numVars; ++ivar)
  {        
  
    // Writing mesh
    //if (procID() == 0)
    { 
      char const *dataspace_names[3] = {"r","phi","theta"};
                               
      for (int dir = 0;dir < SpaceDim; ++dir)
      { 
        IntVect iv_off(IntVect::Zero);
        iv_off[dir]=1;  
        Box b(m_problem_domain.domainBox().smallEnd()*iv_off, 
              m_problem_domain.domainBox().bigEnd()  *iv_off);  
                  
        FArrayBox spacing(b,1);                     
        m_csh->getCellCenters(spacing, b, dir, 0);            
        
        hsize_t dimsf = b.size(dir);
        
        dataspace = H5Screate_simple(1, &dimsf, NULL);

#ifdef H516   
  dataset   = H5Dcreate(h5files[ivar], dataspace_names[dir], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);  
#else
  dataset   = H5Dcreate2(h5files[ivar], dataspace_names[dir], H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);  
#endif

        if (procID() == 0)
          status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, spacing.dataPtr());
          
        H5Sclose(dataspace);
        H5Dclose(dataset);                          
      }                
    }

  
    dataspace = H5Screate_simple(SpaceDim, dimsf, NULL); 

#ifdef H516
  dataset   = H5Dcreate(h5files[ivar], "data", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
  dataset   = H5Dcreate2(h5files[ivar], "data", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif            

    const DisjointBoxLayout& levelDomain = m_UNew.disjointBoxLayout();  
  
    DataIterator dit = m_UNew.dataIterator();
    for (dit.begin(); dit.ok(); ++dit)
    {        
      const FArrayBox& W = plotData[dit()];
      const Box& plotBox = levelDomain[dit()];      
            
      for (dim = 0; dim < SpaceDim; ++dim)
      {
        offset[dim] = plotBox.smallEnd(dim);
        count[dim]  = plotBox.size(dim);
      }
      
      Real * buf = new Real[plotBox.numPts()];
      int    ind;
      
      
    #if CH_SPACEDIM == 2      
      for (j=0;j<count[1];j++)    
      for (i=0;i<count[0];i++)
      {
        ind = j+count[1]*i;
                
        buf[ind] = W.get(IntVect(D_DECL(i+offset[0],j+offset[1],k+offset[2])),ivar);        
      }
    #endif  
      
    #if CH_SPACEDIM == 3 
      for (k=0;k<count[2];k++)    
      for (j=0;j<count[1];j++)    
      for (i=0;i<count[0];i++)
      {
        ind = k + count[2]*(j+count[1]*i);
                
        buf[ind] = W.get(IntVect(D_DECL(i+offset[0],j+offset[1],k+offset[2])),ivar);        
      }
    #endif  
      
      hid_t memdataspace = H5Screate_simple(SpaceDim, count, NULL);
      
      //memspace = H5Dget_space(dataset);
      status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET,  offset, NULL, count, NULL);      
      //status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, memdataspace, memspace, plist_dataset,  buf);
      status = H5Dwrite (dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace, H5P_DEFAULT,  buf);
      
      H5Sclose(memdataspace);
                  
      //H5Sclose(memspace);
      
      delete[]  buf;
    }
    H5Dclose (dataset);
    H5Sclose (dataspace);
    H5Fclose(h5files[ivar]);        
  }
  
#ifdef CH_MPI  
  H5Pclose(plist_file);
  //H5Pclose(plist_dataset);
#endif    

  delete[] h5files;
  

}
