#include <iostream>
using std::cout;
using std::cin;
using std::cerr;
using std::endl;
#include <iomanip>
#include <cmath>
#include <cstdio>
#include <string>
#include <sys/time.h>
#include <libgen.h>
#include <limits>


#include "CHOMBO_VERSION.H"

#include "REAL.H"
#include "Box.H"
#include "Vector.H"
#include "IntVectSet.H"
#include "BRMeshRefine.H"
#include "DisjointBoxLayout.H"
#include "CH_HDF5.H"
#include "AMRIO.H"
#include "NodeAMRIO.H"
#include "AMRLevel.H"
#include "parstream.H"
#include "Tuple.H"
#include "BoxIterator.H"

#include "AMRLevelIdealMHD.H"
#include "tecplotF_F.H"
#include "PatchMHDAM.H"
#include "PoissonBCMHD.H"
#include "MHDAMDefs.H"
#include "SourceCalculator.H"
#include "CSHandler.H"
#include "EquationSystem.H"
#include "AMR_MHDAM.H"

AMR_MHDAM::AMR_MHDAM()
{
  setDefaultValues();
  m_pSourCalc = NULL;
  m_maxChkFiles = -1;
  m_steps_since_Poisson = 0;
  m_sliceInterval = -1;
  m_numSlices = 0;
}

AMR_MHDAM::~AMR_MHDAM()
{
  clearMemory();
}

void
AMR_MHDAM::define( int                          a_max_level,
                   const Vector<int>&           a_ref_ratios,
                   const Box&                   a_prob_domain,
                   const AMRLevelFactory* const a_amrLevelFact,
                   SourceCalculator* a_SourCalc)
{
  ProblemDomain physdomain(a_prob_domain);
  
  define( a_max_level, a_ref_ratios, physdomain, a_amrLevelFact, a_SourCalc);  
}


void
AMR_MHDAM::define( int                          a_max_level,
                   const Vector<int>&           a_ref_ratios,
                   const ProblemDomain&         a_prob_domain,
                   const AMRLevelFactory* const a_amrLevelFact,
                   SourceCalculator* a_SourCalc)
{
  AMR::define( a_max_level, a_ref_ratios, a_prob_domain, a_amrLevelFact );
  m_step_factor.resize(m_max_level+1,1);
  m_input_step_factor.resize(m_max_level+1,-1);
  m_pSourCalc = a_SourCalc;
}

void AMR_MHDAM::input( ParmParse & a_parser, int a_verbosity )
{
  int i;
  this->verbosity(a_verbosity);
  
  Real fillRatio = 0.75;
  a_parser.query("fill_ratio",fillRatio);
  
  int gridBufferSize =1;
  a_parser.query("grid_buffer_size",gridBufferSize);
  
  // Limit the time step growth
  Real maxDtGrowth = 1.1;
  a_parser.query("max_dt_growth",maxDtGrowth);
  
  // Let the time step grow by this factor above the "maximum" before
  // reducing it
  Real dtToleranceFactor = 1.1;
  a_parser.query("dt_tolerance_factor",dtToleranceFactor);
  
  Vector<int> step_factor;
  if( a_parser.contains( "step_factor" ) )
  {
    a_parser.queryarr("step_factor",step_factor,0,m_max_level+1);
  }  
  
  int PoissonInterval = 0;
  a_parser.query("poisson_interval",PoissonInterval);    
  
  int ilsreinit = -1;
  a_parser.query("lsreinit_interval",ilsreinit); 
  
  a_parser.query("probe_dt", m_probeDt);       

  
  // Set up checkpointing
  int checkpointInterval = 0;
  a_parser.query("checkpoint_interval",checkpointInterval);
  int maxChkFiles = -1;
  a_parser.query("max_chk_files",maxChkFiles);

  // Set up plot file writing
  int plotInterval = 0;
  a_parser.query("plot_interval",plotInterval);
  
  int tecplotInterval = 0;
  a_parser.query("tecplot_interval",tecplotInterval);
  
  // Read which slices should be written  
  int sliceIntreval = -1;
  int numSlices     = 0;
  a_parser.query("slices",numSlices);
  a_parser.query("slice_interval",sliceIntreval);
  if (numSlices < 0) numSlices = 0;
  std::vector<std::string> slice_planes;
  std::vector<Real>        slice_values;
  if (numSlices > 0)
  {
    a_parser.queryarr("slice_planes",slice_planes,0,numSlices);
    a_parser.queryarr("slice_values",slice_values,0,numSlices);
  } 
  
  m_numLines     = 0;
  m_lineInterval = -1;
  a_parser.query("lines",m_numLines);
  a_parser.query("lines_interval",m_lineInterval);
  std::vector<Real> coords; RealVect v;
  for (i=0;i<m_numLines;i++)
  {
    char buffer[20];
    sprintf(buffer,"line%i_coords",i);
    a_parser.queryarr(buffer,coords,0,2*CH_SPACEDIM);
    D_TERM(v[0] = coords[0];,v[1] = coords[1];,v[2] = coords[2];);
    m_lineBgn.push_back(v);
    D_TERM(v[0] = coords[CH_SPACEDIM];,v[1] = coords[CH_SPACEDIM+1];,v[2] = coords[CH_SPACEDIM+2];);
    m_lineEnd.push_back(v);
  }
    
  
  
  // Number of coarse time steps from one regridding to the next
  std::vector<int> regridIntervals;
  a_parser.queryarr("regrid_interval",regridIntervals,0,m_max_level);
    
  
  // Set up output files
  if (a_parser.contains("plot_prefix"))
    {
      std::string prefix;
      a_parser.query("plot_prefix",prefix);
      plotPrefix(prefix);
    }

  if (a_parser.contains("chk_prefix"))
    {
      std::string prefix;
      a_parser.query("chk_prefix",prefix);
      checkpointPrefix(prefix);
    }
    
  int outputChkMapFiles = 1;
  a_parser.query("chk_map_files",outputChkMapFiles);  
  m_outputChkMapFiles = (outputChkMapFiles == 1);
    
  

  
  if (( m_verbosity >= 3 ) && (procID() == 0))
  {
    pout() << "maximum dt growth = " << maxDtGrowth << endl;
    
    pout() << "fill ratio = " << fillRatio << endl;
    pout() << "dt tolerance factor = " << dtToleranceFactor << endl;
    pout() << "regrid interval = ";
    for (size_t i = 0; i < regridIntervals.size(); ++i) pout() << regridIntervals[i] << " ";
    pout() << endl;
    pout() << "checkpoint interval = " << checkpointInterval << endl;
    pout() << "plot interval = " << plotInterval << endl;
  }

  #if CHOMBO_VERSION_MAJOR < 4
  this->fillRatio(fillRatio);
  #endif
  this->regridIntervals(regridIntervals);      
  this->gridBufferSize(gridBufferSize);
  m_PoissonInterval = PoissonInterval;
  m_ilsReinitInterval = ilsreinit;
  
  // Set time step parameters
  maxDtGrow(maxDtGrowth);
  this->dtToleranceFactor(dtToleranceFactor);
  setStepFactor(step_factor);  
  
  // Set output parameters  
  
  this->checkpointInterval(checkpointInterval);
  this->plotInterval(plotInterval);  
  m_tecplot_interval = tecplotInterval;
  m_maxChkFiles      = maxChkFiles;
  
  setSliceParameters(sliceIntreval, numSlices,slice_planes, slice_values);
  
}



void AMR_MHDAM::setStepFactor(const Vector<int>& step_factor)
{
  if (step_factor.size() == 0) return;
  if (m_input_step_factor.size()>step_factor.size())
  {
    MayDay::Error("step_factor array is less than number of levels");
  }
  for (int i=0;i<m_input_step_factor.size();i++)
    m_input_step_factor[i]=step_factor[i];   
}


#ifdef CH_USE_HDF5

void AMR_MHDAM::setupForRestart(HDF5Handle& a_handle,
                                int  a_fixedGridSize)
{
  CH_assert(m_isDefined);
  m_isSetUp = true;

  if(m_verbosity >= 3)
  {
    pout() << "AMR::restart" << endl;
  }
    

  HDF5HeaderData header;
  header.readFromFile(a_handle);

  if(m_verbosity >= 3)
  {
    pout() << "hdf5 header data: " << endl;
    pout() << header << endl;
  }

  // read max level
  if(header.m_int.find("max_level") == header.m_int.end())
  {
    MayDay::Error("AMR::restart: checkpoint file does not contain max_level");
  }
  // note that this check should result in a warning rather than an error,
  // because you should be able to restart with a different number of levels
  // (DFM 2/27/02)
  int max_level_check = header.m_int ["max_level"];
  if(max_level_check != m_max_level)
  {        
    if (procID() == 0)
    {
      MayDay::Warning("Warning. AMR::restart: checkpoint file inconsistent with inputs to define");
      MayDay::Warning("         See pout files for details.");
      //pout() << "Warning. AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      //pout() << "         max level input to define = " << m_max_level << endl;
      //pout() << "         max level in checkpoint = " << max_level_check << endl;
    }  
  }

  if(m_verbosity >= 2)
  {
    pout() << "read max_level = " << m_max_level << endl;
  }
  // read finest level
  if(header.m_int.find("num_levels") == header.m_int.end())
  {
    MayDay::Error("AMR::restart: checkpoint file does not contain num_levels");
  }

  int num_levels = header.m_int ["num_levels"];

  if(m_verbosity >= 2)
  {
    pout() << "read num_levels = " << num_levels << endl;
  }

  m_finest_level = num_levels - 1;
  m_finest_level_old = m_finest_level;
  if(m_finest_level > m_max_level)
  {
    if (procID() == 0) pout() << "         number of levels will be decreased" << endl;
    m_finest_level = m_max_level;// Restart program with decreased number of levels
    m_finest_level_old = m_finest_level;

    //pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
    //pout() << "numlevels input to define = " << m_max_level + 1<< endl;
    //pout() << "numlevels in checkpoint = " << num_levels << endl;
    //MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
  }

  if(m_verbosity >= 2)
  {
    pout() << "set finest_level = " << m_finest_level << endl;
  }

  if(m_finest_level > m_max_level)
  {
    MayDay::Error("AMR::restart: finest_level > max_level");
  }

  if(header.m_int.find("iteration") == header.m_int.end())
  {
    MayDay::Error("AMR::restart: checkpoint file does not contain iteration");
  }

  m_cur_step = header.m_int ["iteration"];
  if(m_verbosity >= 2)
  {
    pout() << "read cur_step = " << m_cur_step << endl;
  }
  m_restart_step = m_cur_step;
  if(m_verbosity >= 2)
  {
    pout() << "set restart_step = " << m_restart_step << endl;
  }
  
  #if CHOMBO_VERSION_MAJOR >= 4
  m_fixedGridSize     = a_fixedGridSize;
  s_step = m_cur_step;
  #endif

  if(header.m_real.find("time") == header.m_real.end())
  {
    MayDay::Error("AMR::restart: checkpoint file does not contain time");
  }
  m_cur_time = header.m_real["time"];
  if(m_verbosity >= 2)
  {
    pout() << "read cur_time = " << m_cur_time << endl;
  }

  // All regrid intervals are taken from input file.

  // read regrid intervals
  //if(header.m_int.find("regrid_interval_0") == header.m_int.end())
  //  {
  //    MayDay::Error("AMR::restart: Error: regrid intervals are not in the checkpoint file");
  //  }
  //else
  //  {
  //    //hey if max_level > 10^98, we have a problem
  //    char level_str[100];
  //    //don't go all the way to max level because you don't have to
  //    for(int level = 0; level < m_max_level; ++level)
  //      {
  //        sprintf(level_str, "%d", level);
  //
  //        const std::string label = std::string("regrid_interval_") + level_str;
  //        if(header.m_int.find(label) == header.m_int.end())
  //          {
  //            pout() << "checkpoint file does not have " << label << endl;
  //		  pout() << "    It will be taken from the input file. regrid_intervals[" << level <<"] = " << m_regrid_intervals[level] <<  endl;			
  //            //MayDay::Error("not enough regrid intervals in amr::restart");
  //          }
  //        else
  //          {
  //            m_regrid_intervals[level] = header.m_int [label];
  //          }
  //      }
  //  }

  // read physics class header data
  for(int level = 0; level <= m_finest_level; ++level)
  {
    m_amrlevels[level]->readCheckpointHeader(a_handle);
    m_amrlevels[level]->readCheckpointLevel(a_handle);
  }

  for(int level = 0; level < m_finest_level; ++level)
  {
    int refratio_test = m_amrlevels[level]->refRatio();
    if(refratio_test != m_ref_ratios[level])
    {
      pout() << "AMR::restart: checkpoint file inconsistent with inputs to define " << endl;
      pout() << "for level " << level << endl;
      pout() << "refratio input to define = " << m_ref_ratios[level] << endl;
      pout() << "refratio in checkpoint = "  << refratio_test << endl;
      MayDay::Error("AMR::restart: checkpoint file inconsistent with inputs to define ");
    }
  }

  // maintain time steps
  m_dt_new.resize(m_max_level+1);
  m_dt_cur.resize(m_max_level+1);
  for(int level = 0; level <= m_finest_level; ++level)
  {
    m_dt_new[level] = m_amrlevels[level]->dt();
    m_dt_cur[level] = m_dt_new[level];
  }
  assignDt();
  // maintain steps since regrid
  m_steps_since_regrid.resize(m_max_level+1,0);
  //restart cell updates(we could also output them to the chk file)
  m_cell_updates.resize(m_max_level+1,0);

  #if CHOMBO_VERSION_MAJOR <= 3  
  // final thing to do -- call initialGrid and initialData on undefined levels
  // (just in case there are setup things which need to be done there
  // (DFM 2/27/02)
  for (int level= m_finest_level+1; level <= m_max_level; ++level)
  {
    m_amrlevels[level]->initialGrid(Vector<Box>());
    m_amrlevels[level]->initialData();
  }
  #endif

}
#endif

void
AMR_MHDAM::writePlotFile() const
{
  CH_TIMERS("writePlotFile");
  CH_assert(m_isDefined);
  if(m_verbosity >= 3)
    {
      pout() << "AMR_MHDAM::writePlotFile" << endl;
    }
    
  timeval time1,time2;
  gettimeofday(&time1,NULL);  
  
#ifdef CH_USE_HDF5
  char iter_str[80];
  if(m_cur_step < 10){
    sprintf(iter_str,
            "%s000%d.%dd.hdf5",
            m_plotfile_prefix.c_str(), m_cur_step, SpaceDim);
  }
  else if( m_cur_step < 100){
    sprintf(iter_str,
            "%s00%d.%dd.hdf5",
            m_plotfile_prefix.c_str(), m_cur_step, SpaceDim);
  }
  else if( m_cur_step < 1000){
    sprintf(iter_str,
            "%s0%d.%dd.hdf5",
            m_plotfile_prefix.c_str(), m_cur_step, SpaceDim);
  }
  else
    {
      sprintf(iter_str,
              "%s%d.%dd.hdf5",
              m_plotfile_prefix.c_str(), m_cur_step, SpaceDim);
    }
  if(m_verbosity >= 2)
    {
      pout() << "plot file name = " << iter_str << endl;
    }
    
  CH_TIMER("hdf5_inputs", t1);
  CH_START(t1);  

  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if(m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // Writing input parameters
  /*H5E_auto_t efunc; void* edata; // turn auto error messaging off

#ifdef H516
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  H5Gunlink(handle.groupID(), m_inputFileName.c_str()); //removes a pre-existing dataset.
  H5Eset_auto(efunc, edata);
#else
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);      
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(handle.groupID(), m_inputFileName.c_str(),H5P_DEFAULT); //removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif
*/
  
  
  hid_t s_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(s_type, m_inputFileContent.length()); //extra requirement for strings
  hid_t aid  = H5Screate(H5S_SCALAR);

  char* strbuffer = new char[m_inputFileName.size()+1];
  strcpy(strbuffer,m_inputFileName.c_str());
  char* baseFileName = basename(strbuffer);

#ifdef H516
  hid_t problemdataset   = H5Dcreate(handle.groupID(), baseFileName,  s_type, aid, H5P_DEFAULT);
#else
  hid_t problemdataset   = H5Dcreate2(handle.groupID(), baseFileName,  s_type, aid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  if(problemdataset < 0) MayDay::Error("AMR::writePlotFile: Problems with writing input file to HDF5");

  delete[] strbuffer;

  char* tmp = (char*)m_inputFileContent.c_str();
  if (procID() == 0)
    H5Dwrite(problemdataset, s_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

  H5Sclose(aid);
  H5Tclose(s_type);
  H5Dclose(problemdataset);
  
  
  CH_STOP(t1);  

  // write physics class header data
  m_amrlevels[0]->writePlotHeader(handle);

  // write physics class per-level data
  for(int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writePlotLevel(handle);
    }

  CH_TIMER("closing_file", t2);
  CH_START(t2);  
  handle.close();
  CH_STOP(t2);  
  
  std::string mapped_prefix(iter_str); 
  std::string::size_type hdf5pos = mapped_prefix.find(".hdf5");
  mapped_prefix.erase(hdf5pos,strlen(".hdf5"));      
  writeMappedGeometry(mapped_prefix);
    
      
  
  if ((abs(m_verbosity) >= 1) &&(procID() == 0)) 
  {
    gettimeofday(&time2,NULL);
    char buf[50];
    
    double stepTime = (double)(time2.tv_sec-time1.tv_sec) + 1e-6*((double)(time2.tv_usec-time1.tv_usec));
    
    sprintf(buf,"%.4g",stepTime);
    pout()  << "plot file writing time: " << buf << " sec"  << endl;    
  }
  
    


#endif
}

int AMR_MHDAM::getCurrentStep() const
{
  return m_cur_step;
}

void AMR_MHDAM::writeCheckpointFile() const
{
  CH_assert(m_isDefined);
  if(m_verbosity >= 3)
    {
      pout() << "AMR::writeCheckpointFile" << endl;
    }
    
  timeval time1,time2;
  gettimeofday(&time1,NULL);

  char iter_str[100];iter_str[0] = 0;
  //sprintf(iter_str,
  //        "%s%06i.%dd.hdf5",
  //        m_checkpointfile_prefix.c_str(), m_cur_step, SpaceDim );

  sprintf(iter_str,
          "%s%d.%dd.hdf5",
          m_checkpointfile_prefix.c_str(), m_cur_step, SpaceDim );

  if(m_verbosity >= 2)
    {
      pout() << "checkpoint file name = " << iter_str << endl;
    }

#ifdef CH_USE_HDF5
  HDF5Handle handle(iter_str, HDF5Handle::CREATE);

  // write amr data
  HDF5HeaderData header;
  header.m_int ["max_level"]  = m_max_level;
  header.m_int ["num_levels"] = m_finest_level + 1;
  header.m_int ["iteration"]  = m_cur_step;
  header.m_real["time"]       = m_cur_time;

  for(int level = 0; level < m_regrid_intervals.size(); ++level)
    {
      char headername[100];
      sprintf(headername, "regrid_interval_%d", level);
      header.m_int[headername] = m_regrid_intervals[level];
    }
  // should steps since regrid be in the checkpoint file?
  header.writeToFile(handle);

  if(m_verbosity >= 3)
    {
      pout() << header << endl;
    }

  // Writing input parameters
  //H5E_auto_t efunc; void* edata; // turn auto error messaging off

#ifdef H516
  H5E_auto_t efunc; void* edata; // turn auto error messaging off
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  H5Gunlink(handle.groupID(), m_inputFileName.c_str()); //removes a pre-existing dataset.
  H5Eset_auto(efunc, edata);
#else
  H5E_auto2_t efunc; void* edata; // turn auto error messaging off
  H5Eget_auto2(H5E_DEFAULT,&efunc, &edata);      
  H5Eset_auto2(H5E_DEFAULT,NULL, NULL);
  H5Ldelete(handle.groupID(), m_inputFileName.c_str(),H5P_DEFAULT); //removes a pre-existing dataset.
  H5Eset_auto2(H5E_DEFAULT,efunc, edata);
#endif

  hid_t s_type = H5Tcopy(H5T_C_S1);
  H5Tset_size(s_type, m_inputFileContent.length()); //extra requirement for strings
  hid_t aid  = H5Screate(H5S_SCALAR);

  char* strbuffer = new char[m_inputFileName.size()+1];
  strcpy(strbuffer,m_inputFileName.c_str());
  char* baseFileName = basename(strbuffer);

#ifdef H516
  hid_t problemdataset   = H5Dcreate(handle.groupID(), baseFileName,  s_type, aid, H5P_DEFAULT);
#else
  hid_t problemdataset   = H5Dcreate2(handle.groupID(), baseFileName,  s_type, aid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

  if(problemdataset < 0) MayDay::Error("AMR::writeCheckpointFile: Problems with writing input file to HDF5");

  delete[] strbuffer;

  char* tmp = (char*)m_inputFileContent.c_str();
  if (procID() == 0)
    H5Dwrite(problemdataset, s_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

  H5Sclose(aid);
  H5Tclose(s_type);
  H5Dclose(problemdataset);
  
  CoordinateSystemHandler * csh = static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->getpatchMHDAM()->getPhysProblem()->coordinateSystem();
  csh->writeGeomInfo(handle.fileID());


  // write physics class data
  for(int level = 0; level <= m_finest_level; ++level)
    {
      m_amrlevels[level]->writeCheckpointHeader(handle);
      m_amrlevels[level]->writeCheckpointLevel(handle);
    }
  handle.close();
#endif

  if (m_outputChkMapFiles == true)
  {
    std::string mapped_prefix(iter_str); 
    std::string::size_type hdf5pos = mapped_prefix.find(".hdf5");
    mapped_prefix.erase(hdf5pos,strlen(".hdf5"));      
    writeMappedGeometry(mapped_prefix);
  }


  if ((m_maxChkFiles > 0) && ((int)m_lastChkFiles.size() > m_maxChkFiles))
  {
    std::string fName(m_lastChkFiles.front());
    if (procID() == 0) remove(fName.c_str());
    m_lastChkFiles.pop();

    //fName = m_lastKineticFiles.front();
    //if (procID() == 0) remove(fName.c_str());
    //m_lastKineticFiles.pop();
  }
  m_lastChkFiles.push(std::string(iter_str));
  
  gettimeofday(&time2,NULL);
  
  char buf[50];double stepTime;
  
  if ((abs(m_verbosity) >= 1) &&(procID() == 0)) 
  {
    gettimeofday(&time2,NULL);    
    
    stepTime = (double)(time2.tv_sec-time1.tv_sec) + 1e-6*((double)(time2.tv_usec-time1.tv_usec));
    
    sprintf(buf,"%.4g",stepTime);
    pout()  << "chombo checkpoint file writing time: " << buf << " sec"  << endl;    
  }

  //iter_str[0] = 0;
  sprintf(iter_str, "n%06i.h5",  m_cur_step );
  if (m_pSourCalc != NULL) 
  {
    time1 = time2;
    
    m_pSourCalc->writeCheckpointFile(iter_str);
    gettimeofday(&time2,NULL);
    
    if ((abs(m_verbosity) >= 1) &&(procID() == 0)) 
    {
      gettimeofday(&time2,NULL);    
      
      stepTime = (double)(time2.tv_sec-time1.tv_sec) + 1e-6*((double)(time2.tv_usec-time1.tv_usec));
      
      sprintf(buf,"%.4g",stepTime);
      pout()  << "source calculator file writing time: " << buf << " sec"  << endl;    
    }
            
    //m_lastKineticFiles.push(std::string(iter_str));
  }


}

void AMR_MHDAM::writeMappedGeometry(const string& a_fileRoot, IntVect a_outputGhost) const
{
  CoordinateSystemHandler * csh = static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->getpatchMHDAM()->getPhysProblem()->coordinateSystem();
  
  if ((csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian) ||
      (csh->coordinateSystem() == CoordinateSystemHandler::CS_Axisymmetric) ||
      (csh->coordinateSystem() == CoordinateSystemHandler::CS_Cylindrical)) return;
      
  Vector<DisjointBoxLayout> vectGrids;
    
  
  // now create node-centered data for geometric info
  Vector<LevelData<NodeFArrayBox>* > vectNodeLoc(m_finest_level+1, NULL);
  for(int level = 0; level <= m_finest_level; ++level)
    {
      const LevelData<FArrayBox> & levelData = static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->getStateNew();
      const DisjointBoxLayout& levelGrids = levelData.disjointBoxLayout();
      vectGrids.push_back(levelGrids);
      
      // use same ghosting as cell-centered data used
      //IntVect ghostVect = levelData.ghostVect();
      vectNodeLoc[level] = new LevelData<NodeFArrayBox>(levelGrids,
                                                        SpaceDim,
                                                        a_outputGhost);

      LevelData<NodeFArrayBox>& levelNodeData = *vectNodeLoc[level];      
      DataIterator dit = levelGrids.dataIterator();
      for (dit.begin(); dit.ok(); ++dit)
        {          
          NodeFArrayBox& thisNodeFAB = levelNodeData[dit];
          
          csh->getNodeCoordsCartesian(thisNodeFAB, level);         
        } // end loop over grids

    } // end loop over levels

  // create names 
  Vector<string> locationNames(SpaceDim);
  D_EXPR(locationNames[0] = "x",
         locationNames[1] = "y",
         locationNames[2] = "z");

  char iter_str[80];
  sprintf(iter_str, "%s.map.hdf5", a_fileRoot.c_str());
  string gridInfoFileName(iter_str);
  
  const ProblemDomain& domain = m_amrlevels[0]->problemDomain();
  
  Real dx = 1.0, dt = 1.0;int numLevels = m_finest_level+1;
 
  // now call nodal WriteAMRHierarchy function...
  WriteAMRHierarchyHDF5(gridInfoFileName,
                        vectGrids,
                        vectNodeLoc,
                        locationNames,
                        domain.domainBox(),
                        dx, 
                        dt, 
                        m_cur_time,
                        m_ref_ratios,
                        numLevels);

  // clean up after ourselves here
  for (int level=0; level<=m_finest_level; level++)
    {
      if (vectNodeLoc[level] != NULL)
        {
          delete vectNodeLoc[level];
          vectNodeLoc[level] = NULL;
        }
    }
}

void AMR_MHDAM::writeTecPlotFile()
{
  CH_assert(m_isDefined);
  if (m_verbosity >= 3)
  {
    pout() << "AMR_MHDAM::writeTecPlotFile" << endl;
  }
  
  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(m_amrlevels[0]);
  
/*  if (numProc() == 1)
  {
    int i;
    string variables_string;
    if (CH_SPACEDIM == 2) variables_string = "VARIABLES=\"X\" \"Y\"";
    if (CH_SPACEDIM == 3) variables_string = "VARIABLES=\"X\" \"Y\" \"Z\"";
    
    EquationSystem  * eqSys = Level0->getpatchMHDAM()->getPhysProblem()->equationSystem();

    int nComp   = eqSys->numStates();
    int numPrim = eqSys->numPrimitives();
    //nComp = 1;

    for (i=0;i<nComp;i++)
    if (i<numPrim)
      variables_string+=(" \""+eqSys->primitiveNames()[i]+"\"");

    if (Level0->output_density_gradient()) variables_string+=(" \"grad_rho\"");
    if (Level0->output_B_gradient()) variables_string+=(" \"grad_B\"");

    char file_name[1000];
    sprintf(file_name,"%s%06i.dat",m_plotfile_prefix.c_str(),m_cur_step);
    FILE* tecplot_file=fopen(file_name,"w");
    fprintf(tecplot_file,"TITLE = \"t=%.25g step=%i\"\n",m_cur_time,m_cur_step);
    fprintf(tecplot_file,"%s\n",variables_string.c_str());
    // write physics class per-level data
    for(int level = 0; level <= m_finest_level; ++level)
     {              
#ifdef CH_MPI       
       MPI_File tf;
       MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &tf);
       MPI_File_set_view(tf, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
       static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->writeLevelforTecPlot(tf);
       MPI_File_close(&tf);
#else 
       static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->writeLevelforTecPlot_LZ(tecplot_file,nComp);
#endif      
     }

    fclose(tecplot_file);    
    return;  
  }*/
  
  // MPI version of output.
  
  char file_name[100],buffer[1000];
  sprintf(file_name,"%s%06i.dat",m_plotfile_prefix.c_str(),m_cur_step);

#ifdef CH_MPI
  MPI_File tecplot_file;
  
  // Delete old file if ut exists
  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &tecplot_file);
  MPI_File_close(&tecplot_file);

  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &tecplot_file);
  //MPI_File_set_view(tecplot_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
  MPI_File_set_view(tecplot_file, 0, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
#else
  FILE*     tecplot_file=fopen(file_name,"w");
#endif

  sprintf(buffer,"TITLE = \"t=%.25g step=%i\"\n",m_cur_time,m_cur_step);
  
  // write physics class per-level data
  for(int level = 0; level <= m_finest_level; ++level)
   {
     int nPatches = static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->getStateNew().disjointBoxLayout().size();
     sprintf(buffer+strlen(buffer),"# Level %i, number of patches: %i\n",level, nPatches);
   }

#ifdef CH_MPI
  if (procID() == 0)
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else
  fprintf(tecplot_file,"%s",buffer);
#endif

  // write physics class per-level data
  static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->writeLevel0forTecPlot(tecplot_file);
  for(int level = 1; level <= m_finest_level; ++level)
   {
     static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->writeLevelforTecPlot(tecplot_file);
   }

#ifdef CH_MPI
  MPI_File_close(&tecplot_file);
#else
  fclose(tecplot_file);
#endif

  //static_cast<AMRLevelIdealMHD*>(m_amrlevels[1])->MedvedevOutput(m_cur_step);

}


void AMR_MHDAM::setInputFile(const char* a_inFile)
{
  m_inputFileName = a_inFile;
  std::ifstream infile(a_inFile);
  std::string buffer;

  while (!infile.eof())
  {
    getline(infile,buffer);
    m_inputFileContent+=buffer;
    m_inputFileContent+="\n";
  }
  infile.close();
}

void AMR_MHDAM::setSliceParameters(int a_sliceIntreval, int a_numSlices,
              std::vector<std::string>  a_slicePlanes,
              std::vector<Real> a_sliceValues)
{
  m_sliceInterval = a_sliceIntreval;
  m_numSlices     = a_numSlices;
  m_sliceValues   = a_sliceValues;

  int i;
  if (m_numSlices > 0)
  {
    CH_assert((m_numSlices == a_slicePlanes.size()) && (m_numSlices == a_sliceValues.size()) );
    m_slicePlanes.resize(m_numSlices);
    for (i=0; i<m_numSlices;i++)
    {
      char plane = a_slicePlanes[i][0];
      if ((plane == 'X') || (plane == 'x'))
        m_slicePlanes[i] = 'X';
      else
      if ((plane == 'Y') || (plane == 'y'))
        m_slicePlanes[i] = 'Y';
      else
      if ((plane == 'Z') || (plane == 'z'))
        m_slicePlanes[i] = 'Z';
      else MayDay::Error("Error in mhdam.slice_planes");
    }
  }
}

void AMR_MHDAM::writeSlices()
{
  for (int i=0; i<m_numSlices;i++)
    writeOneSlice(m_slicePlanes[i],m_sliceValues[i]);

}

void AMR_MHDAM::writeOneSlice(char plane, Real sValue) const
{
  if (CH_SPACEDIM != 3) return;

  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(m_amrlevels[0]);
  PhysProblem * phPr = Level0->getpatchMHDAM()->getPhysProblem();


  int iSlice = (int)floor(sValue+0.5);

  char file_name[100],buffer[1000];
  sprintf(file_name,"%ss%c%i_%06i.dat",m_plotfile_prefix.c_str(),plane,iSlice,m_cur_step);

#ifdef CH_MPI
  MPI_File tecplot_file;
  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_DELETE_ON_CLOSE, MPI_INFO_NULL, &tecplot_file);
  MPI_File_close(&tecplot_file);

  MPI_File_open(Chombo_MPI::comm,file_name,MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &tecplot_file);
  MPI_File_set_view(tecplot_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);
#else
  FILE*     tecplot_file=fopen(file_name,"w");
#endif

  sprintf(buffer,"TITLE = \"t=%.25g step=%i\"\n",m_cur_time,m_cur_step);  

#ifdef CH_MPI
  if (procID() == 0)
    MPI_File_write_shared(tecplot_file, buffer, strlen(buffer), MPI_CHAR, MPI_STATUS_IGNORE);
#else
  fprintf(tecplot_file,"%s",buffer);
#endif

  CoordinateSystemHandler * csh = static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->getpatchMHDAM()->getPhysProblem()->coordinateSystem();
  
  Real PI = 3.14159265358979323846;

  // write physics class per-level data
  for(int level = 0; level <= m_finest_level; ++level)
   {
     static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->writeSlice(tecplot_file, plane, sValue);
     if ((csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) && (plane == 'Y'))
       static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->writeSlice(tecplot_file, plane, fmod(sValue+PI,2*PI));
         
   }

#ifdef CH_MPI
  MPI_File_close(&tecplot_file);
#else
  fclose(tecplot_file);
#endif
}

void AMR_MHDAM::writeLines()
{
  for (int i=0; i<m_numLines;i++)
    writeOneLine(i,m_lineBgn[i],m_lineEnd[i]);
}


void AMR_MHDAM::writeOneLine(int iLine, RealVect a_Bgn, RealVect a_End) const
{

  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(m_amrlevels[0]);
  CoordinateSystemHandler * csh = static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->getpatchMHDAM()->getPhysProblem()->coordinateSystem();      
  EquationSystem  * eqSys = Level0->getpatchMHDAM()->getPhysProblem()->equationSystem();
  PhysProblem * phPr = Level0->getpatchMHDAM()->getPhysProblem();  
  
  RealVect crvBgn, crvEnd;
  char buffer[1000];
  
  csh->transCartesianCoordsToCurv(crvBgn, a_Bgn);
  csh->transCartesianCoordsToCurv(crvEnd, a_End);
  
  crvBgn = a_Bgn;
  crvEnd = a_End;
  
  IntVect iBgn, iEnd;
  csh->getClosestCell(iBgn, crvBgn, 0);
  csh->getClosestCell(iEnd, crvEnd, 0);
  
  const ProblemDomain& domain = m_amrlevels[0]->problemDomain();
    
  if (!domain.domainBox().contains(iBgn))
  {
    if (procID() == 0) pout() << "AMR_MHDAM::writeOneLine: Begining of the line (" << a_Bgn << ") must be inside of the problem domain" << endl;
    return;
  }
  if (!domain.domainBox().contains(iEnd))
  {
    if (procID() == 0) pout() << "AMR_MHDAM::writeOneLine: End of the line (" << a_End << ") must be inside of the problem domain" << endl;
    return;
  }
  
  int dirLine = -1,dir;
  for (dir = 0; dir < SpaceDim; dir++)
  {
    if (D_TERM(,(iBgn[(dir+1)%SpaceDim] == iEnd[(dir+1)%SpaceDim]), && (iBgn[(dir+2)%SpaceDim] == iEnd[(dir+2)%SpaceDim])))
    dirLine =  dir;
  }
  
  if (dirLine == -1)
  {
    if (procID() == 0) pout() << "AMR_MHDAM::writeOneLine: Line must coincide with the coordinate line" << endl;
    return;
  }
  
  IntVect iv_off (IntVect::Zero);
  iv_off[dirLine] = 1;
  
  IntVect ivLow,ivHi;
  D_TERM(
    ivLow[dirLine] = domain.domainBox().smallEnd()[dirLine];
    ivHi [dirLine] = domain.domainBox().bigEnd()[dirLine]; 
    ,    
    ivLow[(dirLine+1)%SpaceDim] = iBgn[(dirLine+1)%SpaceDim];
    ivHi [(dirLine+1)%SpaceDim] = iBgn[(dirLine+1)%SpaceDim];
    ,
    ivLow[(dirLine+2)%SpaceDim] = iBgn[(dirLine+2)%SpaceDim];
    ivHi [(dirLine+2)%SpaceDim] = iBgn[(dirLine+2)%SpaceDim];);  
        
  
  Box lineBox (ivLow, ivHi);  
    
  LevelData<FArrayBox>* lineData = new LevelData<FArrayBox>[m_finest_level+1];
    
  int level, i;  
  
  for(level = 0; level <= m_finest_level; ++level)
  {
    Vector<Box>  line_boxes;
    Vector<int>  procs;    
    
    if (level == 0)
    {
      line_boxes.push_back(lineBox);
      procs.push_back(0);            
    } else    
    {
      lineBox.refine(m_ref_ratios[level-1]);
      csh->getClosestCell(iBgn, crvBgn, level);      
      csh->getClosestCell(iEnd, crvEnd, level);
      
      CH_assert(D_TERM(,(iBgn[(dirLine+1)%SpaceDim] == iEnd[(dirLine+1)%SpaceDim]), && (iBgn[(dirLine+2)%SpaceDim] == iEnd[(dirLine+2)%SpaceDim])));
      
      D_TERM(  ,
        lineBox.setSmall((dirLine+1)%SpaceDim,iBgn[(dirLine+1)%SpaceDim]);
        lineBox.setBig  ((dirLine+1)%SpaceDim,iBgn[(dirLine+1)%SpaceDim]);,
        lineBox.setSmall((dirLine+2)%SpaceDim,iBgn[(dirLine+2)%SpaceDim]);
        lineBox.setBig  ((dirLine+2)%SpaceDim,iBgn[(dirLine+2)%SpaceDim]););              

#if CHOMBO_VERSION_MAJOR < 4                
      const Vector<Box> & level_boxes = m_amrlevels[level]->boxes();
#else       
      const Vector<Box>   level_boxes = m_amrlevels[level]->boxes().getVBox();
#endif
      
      for (i=0;i<level_boxes.size();i++)
      {
        Box b = level_boxes[i] & lineBox;
        if (b.isEmpty()) continue;
        line_boxes.push_back(b);
        procs.push_back(0);            
      }
    }
    DisjointBoxLayout strips(line_boxes, procs, m_amrlevels[level]->problemDomain());
    lineData[level].define(strips, eqSys->numStates());
    
    AMRLevelIdealMHD* myLevel=static_cast<AMRLevelIdealMHD*>(m_amrlevels[level]);
    const LevelData<FArrayBox>& U = myLevel->getStateNew();
    
    U.copyTo(U.interval(), lineData[level], U.interval());
    
  }
      
  if (procID() == 0)
  {
 
    char file_name[1000];
    sprintf(file_name,"line%i_%06i.dat",iLine,m_cur_step);
    FILE* tecplot_file=fopen(file_name,"w");
    fprintf(tecplot_file,"# Data along line (%.25g,%.25g",a_Bgn[0],a_Bgn[1]);
    if (SpaceDim > 2) fprintf(tecplot_file,",%.25g",a_Bgn[2]);
    fprintf(tecplot_file,") - (%.25g,%.25g",a_End[0],a_End[1]);
    if (SpaceDim > 2) fprintf(tecplot_file,",%.25g",a_End[2]);
    
    fprintf(tecplot_file,")\nTITLE = \"t=%.25g step=%i\"\n",m_cur_time,m_cur_step);        
        
    
    if (SpaceDim == 2) 
      strcpy(buffer,"VARIABLES=\"Distance\" \"X\" \"Y\"");
    else if (SpaceDim == 3) 
      strcpy(buffer,"VARIABLES=\"Distance\" \"X\" \"Y\" \"Z\"");
    else MayDay::Error("AMR_MHDAM::writeOneLine unsuported dimension");
    
    int numVars = eqSys->numPrimitives() +
                Level0->getpatchMHDAM()->getPhysProblem()->numPlotVars() +
                Level0->getpatchMHDAM()->getPhysProblem()->getSourceCalculator()->numPlotVars();
                
    Vector<string> varNames;
    varNames.reserve(numVars);        
                          
    varNames = eqSys->primitiveNames();
    for (i = 0; i < (int)varNames.size(); ++i)
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    varNames = Level0->getpatchMHDAM()->getPhysProblem()->plotNames();
    for (i = 0; i < (int)varNames.size(); ++i)
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    varNames = Level0->getpatchMHDAM()->getPhysProblem()->getSourceCalculator()->plotNames();
    for (i = 0; i < (int)varNames.size(); ++i)
    {
      strcat(buffer," \"");
      strcat(buffer,varNames[i].c_str());
      strcat(buffer,"\"");      
    }
    strcat(buffer,"\n");      
    fprintf(tecplot_file,"%s",buffer);
        

    fprintf(tecplot_file,"ZONE, DATAPACKING=POINT\n");		  
    std::string auxdata;
    phPr->auxDataTecplot(auxdata,m_cur_time,1);
    if (auxdata.size()>0)
    {
      fprintf(tecplot_file,"%s\n",auxdata.c_str());
    }
    
    Level0->writeOneLine(tecplot_file, m_amrlevels[0]->problemDomain().domainBox(), dirLine, lineData);
    
    fclose(tecplot_file);  
  }
     
  delete[] lineData;  
    
}

void AMR_MHDAM::write1DProbe()
{
  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(m_amrlevels[0]);
  PhysProblem * phPr = Level0->getpatchMHDAM()->getPhysProblem();  
  
  Interval dummyInt; int level, iProbe;
  
            
  for (level = 0; level <= m_finest_level; ++level) static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->fillGhostCellsCons(dummyInt);
  
  for (iProbe = 0; iProbe < m_pfi.numProbes(); ++iProbe)
  {
    const ProbeFilesInfo::ProbeInfo & probe = m_pfi.getProbe(iProbe);
    
    for (level = m_finest_level; level >= 0; --level)
    {
      AMRLevelIdealMHD* MHD_Level=static_cast<AMRLevelIdealMHD*>(m_amrlevels[level]);
      int res = MHD_Level->writeProbeData(probe);
      
      if (res == 1) break;
    }      
  }
  
}

bool AMR_MHDAM::needToWrite1DProbe(Real a_time)
{    
  if (m_pfi.numProbes() == 0) return false;
  
  Real time = a_time;
  if ((time - m_last1dprobe) > m_probeDt)
  {
    m_last1dprobe = time;
    return true;
  } else return false;
  
}


void AMR_MHDAM::setup1DProbeFiles()
{
  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(m_amrlevels[0]);
  PhysProblem * phPr = Level0->getpatchMHDAM()->getPhysProblem();  
  
  m_last1dprobe = -1.0;
  
  int numVars; std::string header;
  phPr->probeFilesHeader(numVars, header);
  m_probeDt = phPr->probeFilesDimensiolessDt(m_probeDt);
  
  m_pfi.setupProbeFiles(numVars, header);

/*    
  if (m_pfi.m_filenames.size() == 0) return;

  int ind;
  if (procID() == 0)
  for (ind=0; ind<m_pfi.m_filenames.size(); ind++)  
  {      
    struct stat stFileInfo;  
    int intStat; 
      
    FILE* plot_file;

    // Attempt to get the file attributes 
    intStat = stat(m_pfi.m_filenames[ind].c_str(),&stFileInfo); 
    plot_file=fopen(m_pfi.m_filenames[ind].c_str(),"a");
    if(intStat == 0)
    {             
      fprintf(plot_file,"\n===========\n");
    } else
    { 
      fprintf(plot_file,"%s\n",header.c_str());      
    }    
    fclose(plot_file);
  }
     
 
  
  for (ind=0; ind<m_pfi.m_filenames.size(); ind++)
  {
#ifdef CH_MPI    
    MPI_File f;
    MPI_File_open(Chombo_MPI::comm, const_cast<char*>(m_pfi.m_filenames[ind].c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &f);
    MPI_File_set_view(f, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);    
#else
    FILE* f = fopen(m_pfi.m_filenames[ind].c_str(),"a");    
#endif    
    m_pfi.m_files.push_back(f);
  }*/
  
}


int outputData(const Vector<LevelData<FArrayBox>* >& vectPhi,
               const Vector<LevelData<FArrayBox>* >& vectRhs,
               const Vector<DisjointBoxLayout>& vectGrids,
               const Vector<ProblemDomain>& vectDomain,
               const Vector<int>& vectRatio,
               Real dxCoarsest,
               int numlevels,
               int cur_step)
{
#ifdef CH_MPI
  MPI_Barrier(Chombo_MPI::comm);
#endif
  int ncomp = vectPhi[0]->nComp();
  

  string phiNameBase("phi");
  string rhsNameBase("rhs");

  Vector<string> vectName(2*ncomp);
  for (int comp=0; comp<ncomp; comp++)
    {
      char labelString[80];
      sprintf(labelString, "%s%d", phiNameBase.c_str(), comp);
      string phiName(labelString);
      vectName[comp] = phiName;
      sprintf(labelString, "%s%d", rhsNameBase.c_str(), comp);
      string rhsName(labelString);
      vectName[ncomp+comp] = rhsName;
    }
  vectName.push_back("procID");
  Box domain = vectDomain[0].domainBox();
  // Real dx = 1.;
  Real dt = 1.;
  Real time = 1.;
  Vector<LevelData<FArrayBox>* > vectPhiAndRHS(numlevels, NULL);
  for(int ilev = 0; ilev < numlevels; ilev++)
    {
      vectPhiAndRHS[ilev] = new LevelData<FArrayBox>(vectGrids[ilev], 2*ncomp+1);

      Interval phiInterval(0,ncomp-1);
      Interval rhsInterval(ncomp, 2*ncomp-1);
      vectPhi[ilev]->copyTo(vectPhi[ilev]->interval(),
                            *vectPhiAndRHS[ilev],
                            phiInterval);
      vectRhs[ilev]->copyTo(vectRhs[ilev]->interval(),
                            *vectPhiAndRHS[ilev],
                            rhsInterval);
	  DataIterator dit = vectPhiAndRHS[ilev]->dataIterator();
	  for(; dit.ok(); ++dit){
		vectPhiAndRHS[ilev]->operator[](dit()).setVal(procID(), 2*ncomp);
	  }
    }
#ifdef CH_USE_HDF5
  char iter_str[100];iter_str[0] = 0;
  sprintf(iter_str, "poissonOut%d.hdf5", cur_step );
  string filename(iter_str);
  WriteAMRHierarchyHDF5(filename,
                        vectGrids,
                        vectPhiAndRHS,
                        vectName,
                        domain,
                        dxCoarsest, dt, time,
                        vectRatio,
                        numlevels);
#endif

  for(int ilev = 0; ilev < numlevels; ilev++)
    delete vectPhiAndRHS[ilev];

  return 0;
}

void AMR_MHDAM::reinitializeLevelSet()
{
  for(int level = 0; level <= m_finest_level; ++level)
  {
    static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->reinitializeLevelSet();             
  }
}

void AMR_MHDAM::ApplyProjectionScheme()
{
/*  int numlevels = MIN(m_max_level-1,m_finest_level);
  int ilev;

  //pout() << "m_steps_since_regrid = ";
  //for(ilev = 0; ilev < m_steps_since_regrid.size(); ilev++)
  // {
  //  int a = m_steps_since_regrid[ilev];
  //  pout() << a << " ";
  //}
  //pout() << endl;

  if (m_PoissonInterval <=0 ) return;
  if (m_steps_since_Poisson < m_PoissonInterval)
  {
    pout() << m_steps_since_Poisson << " < " << m_PoissonInterval << " "<< endl;
    return;
  }

  for(ilev = 0; ilev <= numlevels; ilev++)
  if (m_steps_since_regrid[ilev] == 0)
  {
    pout() << "m_steps_since_regrid[ilev] == 0, level" << ilev << endl;
    return;
  }

  if (m_verbosity >= 3)
    pout() << "AMR_MHDAM::ApplyProjectionScheme" << endl;


  numlevels = m_finest_level + 1;

  int idir;
  AMRLevelIdealMHD* curLevel;
  DataIterator dit;

  PhysProblem* PhPr = static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->getpatchMHDAM()->getPhysProblem();

  DomainGhostBC domghostbc;

  // Set up boundary conditions
  SideIterator si;
  for(idir = 0; idir < SpaceDim; idir++)
  for (si.begin(); si.ok(); ++si)
  {
    Side::LoHiSide side = si();

    BoxGhostBC* bgbc = NULL;
    switch (PhPr->getBCFlags(idir,side))
    {
      case PhysProblem::BC_Periodic   : { bgbc = new PeriodicBC(idir,side);   break;}//bgbc = new AxisBC(idir,side);
      case PhysProblem::BC_Fixed      : { bgbc = new FixedBC(idir,side);      break;}
      case PhysProblem::BC_Continuous : { bgbc = new ContinuousBC(idir,side); break;}
      case PhysProblem::BC_Axis       : { bgbc = new AxisBC(idir,side);       break;}
      case PhysProblem::BC_Mixed      : { MayDay::Error("AMR_MHDAM::ApplyProjectionScheme: BC_Mixed is not supported yet"); break;}
      case PhysProblem::BC_Undefined  : { MayDay::Error("AMR_MHDAM::ApplyProjectionScheme: Undefined boundary condition");  break;}
      default: { MayDay::Error("AMR_MHDAM::ApplyProjectionScheme: Unknown boundary condition"); break;}
    }
    domghostbc.setBoxGhostBC(*bgbc);
    delete bgbc;
  }
  PoissonOp levelop;
  levelop.setDomainGhostBC(domghostbc);

  // Prepare initial data
  Vector<DisjointBoxLayout> vectGrids;
  Vector<ProblemDomain> vectDomain;
  Vector<int> vectRefRatio;
  Vector<Real> vectDx;

  vectGrids.resize(numlevels);
  vectDomain.resize(numlevels);
  vectDx.resize(numlevels);
  vectRefRatio.resize(numlevels);

  Vector<LevelData<FArrayBox>* > phi(numlevels, NULL);
  Vector<LevelData<FArrayBox>* > rhs(numlevels, NULL);
  for(ilev = 0; ilev < numlevels; ilev++)
  {
    curLevel = static_cast<AMRLevelIdealMHD*>(m_amrlevels[ilev]);

    const ProblemDomain& pd       = curLevel->problemDomain();
    DisjointBoxLayout grid = curLevel->getStateNew().disjointBoxLayout();
    
    Real dx = curLevel->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,ilev);

    vectDomain[ilev]   = pd;
    vectDx[ilev]       = dx;
    vectRefRatio[ilev] = curLevel->refRatio();
    DisjointBoxLayout& dbl = vectGrids[ilev];
    dbl.define(grid,pd);

    phi[ilev] = new LevelData<FArrayBox>(curLevel->getStateNew().disjointBoxLayout(), 1, IntVect::Unit);
    rhs[ilev] = new LevelData<FArrayBox>(curLevel->getStateNew().disjointBoxLayout(), 1, IntVect::Zero);
    LevelData<FArrayBox>& phifabs= *phi[ilev];
    for(dit = phifabs.dataIterator(); dit.ok(); ++dit)
    {
      phifabs[dit()].setVal(0.0);
    }

    LevelData<FArrayBox>& rhsfabs= *rhs[ilev];
    const LevelData<FArrayBox>& divBdata = curLevel->getdivB();
    for(dit = divBdata.dataIterator(); dit.ok(); ++dit)
    {
      FArrayBox& rfab = rhsfabs[dit()];
      rfab.copy(divBdata[dit()]);
      rfab.mult(1.0/dx);
    }
  }

  AMRSolver amrSolver(vectGrids,  vectDomain,
                      vectDx,     vectRefRatio,
                      numlevels, 0, &levelop);

  amrSolver.setVerbose(m_verbosity > 0);
  int maxiter = 20;
  amrSolver.setMaxIter(maxiter);
  int numvbot = 1;
  amrSolver.setNumVCyclesBottom(numvbot);

  amrSolver.setTolerance(1e-12);
  amrSolver.setMinIter(10);
  //amrSolver.setOperatorTolerance(1e-10) ;

  amrSolver.solveAMR(phi, rhs);

  //int solve_iterations = 5;
  //for(int i=0; i<solve_iterations; i++) {
  //  amrSolver.solveAMR(phi, rhs);
  //}

  if (m_verbosity >= 3)
  {
    Vector<Real> residualNorm = amrSolver.computeResidualNorm(phi, rhs, 0);
    for(ilev = 0; ilev < numlevels; ilev++)
    if (ilev < (int)residualNorm.size())
    {
      pout() << "Residual on level " << ilev << " is " << residualNorm[ilev] << endl;
    }
  }

  for(ilev = 0; ilev < numlevels; ilev++)
  {
    LevelData<FArrayBox>& phifabs= *phi[ilev];

    // Fill ghost cells, I'm not sure that it is really needed
    if (ilev > 0)
    {
      QuadCFInterp quadCFI(
        phi[ilev]->disjointBoxLayout(),&(phi[ilev-1]->disjointBoxLayout()),
        vectDx[ilev],vectRefRatio[ilev-1],1,vectDomain[ilev]);

      quadCFI.coarseFineInterp(phifabs, *phi[ilev-1]);
    }
    phifabs.exchange(phifabs.interval());

    curLevel = static_cast<AMRLevelIdealMHD*>(m_amrlevels[ilev]);
    curLevel->ApplyProjectionScheme(phifabs);
  }

  outputData(phi, rhs, vectGrids, vectDomain,
                       vectRefRatio, vectDx[0], numlevels, m_cur_step);

  for(ilev = 0; ilev < numlevels; ilev++)
  {
    delete phi[ilev];
    delete rhs[ilev];
  }

  m_steps_since_Poisson = 0;*/
}

void AMR_MHDAM::assignDt()
{
  //AMR::assignDt(); return;
  
  CH_assert(isDefined());

  if (m_verbosity >= 3)
    {
      pout() << "AMR_MHDAM::assignDt" << endl;
    }

  if(m_useSubcycling)
  {
    if (m_fixedDt > 0)
      {
        m_dt_base = m_fixedDt;
      }
    else
      {
        m_dt_base = m_dt_new[0];
                        
        bool use_input_step_factor = (m_input_step_factor[0]>0);

        // Multiply time step for each level by refinement factor, then min
        // across all levels.  Only go to finest_level_old because the higher
        // level dt's have not been set
        for (int level = 1; level <= m_finest_level_old; ++level)
          {
            int ref_factor = 1;

            for (int ilev = 0; ilev < level; ++ilev)
              {
                if (use_input_step_factor) 
                  ref_factor *= m_input_step_factor[ilev];
                else 
                  ref_factor *= m_ref_ratios[ilev];
              }

            Real dt_base_equiv = m_dt_new[level] * ref_factor;
                                    
            m_dt_base = Min(m_dt_base, dt_base_equiv);
          }

        // Check all the actual dt's (scaled by the refinement factor) and only
        // allow a growth of "m_maxDtGrow".
        for (int level = 0; level <= m_finest_level_old; ++level)
          {
            int ref_factor = 1;

            for (int ilev = 0; ilev < level; ++ilev)
              {
                ref_factor *= m_ref_ratios[ilev];
              }

            Real dt_base_equiv = m_dt_cur[level] * ref_factor * m_maxDtGrow;
            m_dt_base = Min(m_dt_base, dt_base_equiv);
          }
      }
          
    Real dt_level;
    Real prev_dt = m_dt_base;

    // refine base time step for all levels
    for (int level = 0; level <= m_max_level; ++level)
    {
      int ref_factor = 1;

      // reset reduction factors as there is no subcycling going on yet
      m_reduction_factor[level] = 1;
              
      ref_factor = (level == 0 ? 1 : m_ref_ratios[level-1]);
            
      if (m_dt_new[level]==0.0)
      {
        // This level just appeared or does not exist, m_dt_new[level] has not been defined yet
        dt_level = prev_dt / ref_factor;
        m_step_factor[level] = ref_factor;
      } else
      {
        if (m_dt_new[level] < prev_dt)
        {        
          int step_factor = (int) ceil(prev_dt/m_dt_new[level]);
          dt_level = prev_dt / step_factor;                
          m_step_factor[level] = step_factor;
        } else
        {
          dt_level = prev_dt;                
          m_step_factor[level] = 1;
        }
      }      
      m_amrlevels[level]->dt(dt_level);
      prev_dt  = dt_level;      
    }
      
    for (int level = m_finest_level+1; level <= m_max_level; ++level)
    {
      m_dt_new[level] = 0.0;      
    }
  }
  else
  {
    ///no subcycling.

    if (m_fixedDt > 0)
      {
        m_dt_base = m_fixedDt;
      }
    else
      {
        m_dt_base = m_dt_new[0];
        Real dt_base_equiv = m_dt_cur[0]*m_maxDtGrow;
        m_dt_base = Min(m_dt_base, dt_base_equiv);
      }
    //everybody gets the same time step if subcyling is turned off
    for (int level = 0; level <= m_max_level; ++level)
      {
        m_amrlevels[level]->dt(m_dt_base);
      }
  }
  
  //if ((m_verbosity >= 2) || ((m_verbosity <= -2)&&(procID() == 0)))
  //{
  //  pout() << "AMR_MHDAM::assignDt ";
  //  for (int level = 0; level <= m_max_level; ++level)
  //    {
  //      pout() << " level" << level << " dt = " << m_amrlevels[level]->dt() << ",";        
  //    }
  //  pout() << endl;
  //}
}

// This function performs a time step of level "a_level" and all finer
// levels.  It does subcycling in time as necessary (see below).  For this
// reason the number of time steps left, "a_stepsLeft", for this level is
// passed in.  If no additional subcycling occurs then this is returned.
// If additional subcycling occurs then this is used to compute the new
// number of time steps remaining (for this level) and this is returned.
int AMR_MHDAM::timeStep(int a_level, int a_stepsLeft, bool a_coarseTimeBoundary)
{
  //AMR::timeStep(a_level, a_stepsLeft, a_coarseTimeBoundary); return(a_stepsLeft);
  
  CH_assert(isDefined());
  CH_assert(isSetUp());

  if (m_verbosity >= 3)
    {
      pout() << "AMR_MHDAM::timeStep(" << a_level << ")" << endl;
    }

  if (a_level < m_max_level)
    {
      if (m_verbosity >= 4)
        {
          pout() << "Regrid (level " << a_level << ") needed - ";
        }

      // regrid if necessary
      if (needToRegrid(a_level,a_stepsLeft))
        {
          if (m_verbosity >= 4)
            {
              pout() << "yes" << endl;
            }

          regrid(a_level);
          
          if (abs(m_verbosity) >= 10)
            {                            
              writePlotFile();              
              //writeTecPlotFile();        
              writeSlices();              
            }
          
        }
      else
        {
          if (m_verbosity >= 4)
            {
              pout() << "no" << endl;
            }
        }
    }

  // If this wasn't just done by the next coarser level, check to see if
  // it is necessary to do additional subcycling in time.
  if ((!a_coarseTimeBoundary) && (m_fixedDt <= 0))
    {
      // The factor by which the current time step at the current level
      // has been divided (so far) for subcycling.
      int maxFactor = m_reduction_factor[a_level];

      // Compute the new subcycling factor for this level and all finer
      // levels and find the maximum
      for (int i = a_level; i <= m_max_level; i++)
        {
          int factor;
          Real dtCur = m_amrlevels[i]->dt();
          Real dtNew = m_dt_new[i];

          // The current factor for level "i"
          factor = m_reduction_factor[i];

          // While the current dt exceeds the new (max) dt by a tolerance
          // double the subcycling factor and half the current dt
          while (dtCur > m_dt_tolerance_factor*dtNew)
            {
              factor *= 2;
              dtCur *= 0.5;
            }

          if (factor > maxFactor)
            {
              maxFactor = factor;
            }
        }

      // More subcycling is necessary
      if (maxFactor > m_reduction_factor[a_level])
        {
          if (m_verbosity >= 3)
            {
              pout() << "  Subcycling --- maxFactor: " << maxFactor << endl;
            }

          // Adjust the number of time steps left for the current level
          a_stepsLeft = (a_stepsLeft+1)*maxFactor/m_reduction_factor[a_level] - 1;

          // Adjust the dt's on this and all finer levels
          for (int i = a_level; i <= m_max_level; i++)
            {
              int factor;

              factor = maxFactor/m_reduction_factor[i];
              m_amrlevels[i]->dt(m_amrlevels[i]->dt()/factor);

              if (m_verbosity >= 4)
                {
                  pout() << "    Level " << i << ": factor: " << factor
                         << " (" << m_reduction_factor[i] << "), dt: "
                         << m_amrlevels[i]->dt() << endl;
                }

              m_reduction_factor[i] = maxFactor;
            }
        }
    }

  // advance this level
  m_amrlevels[a_level]->advance();

  // get the new dt
  Real dt_level = m_amrlevels[a_level]->computeDt();

  // Save the current dt and the new (max) dt.
  m_dt_cur[a_level] = m_amrlevels[a_level]->dt();
  m_dt_new[a_level] = dt_level;

  // increment counter that gives the number of cells updates.
  long long numPts = 0;
  
#if CHOMBO_VERSION_MAJOR < 4  
  const Vector<Box >& levelBoxes = m_amrlevels[a_level]->boxes();
  for (int ll = 0;ll < levelBoxes.size();ll++)
    {
      numPts += levelBoxes[ll].numPts();
    }
#else
    int numptslev = m_amrlevels[a_level]->numPts();
    numPts = (long long)(numptslev);
    //above number is the number of BOXES.  no multiply by the 
    //number of points per box
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        numPts *= m_fixedGridSize;
      }
#endif  

  m_cell_updates[a_level] += numPts;

  if(a_level < m_max_level)
    {
      ++m_steps_since_regrid[a_level];
    }

  // advance the finer levels by subcycling
  if(a_level < m_finest_level)
    {
      int stepsLeft = m_step_factor[a_level+1];
      bool timeBoundary = true;

      while (stepsLeft > 0)
        {
          stepsLeft--;

          // Advance the finer level and take into account possible
          // subcycling by allowing for a change in "stepsLeft".
          //[NOTE: the if() test looks redundant with above, but it is not 
          //       because m_finest_level may change during a regrid();
          //       why you would regrid during a subcycle I don't know. <dbs>]
          if (a_level < m_finest_level)
            stepsLeft = timeStep(a_level+1,stepsLeft,timeBoundary);

          // The first time the next finer level time aligns with the current
          // level time.  After that this is not the case.
          //[NOTE: this if() test _is_ redundant. <dbs>]
          if (timeBoundary == true)
            {
              timeBoundary = false;
            }
        }
    }

  m_amrlevels[a_level]->postTimeStep();

  // Return the (possibly updated) number of time steps left on this level.
  return(a_stepsLeft);
}


void AMR_MHDAM::run(Real a_max_time, int a_max_step)
{
  CH_assert(isDefined());
  CH_assert(isSetUp());
  if(m_verbosity >= 3)
    {
      pout() << "AMR::coarseTimeStep:" << endl;
      pout() << "max_time = " << a_max_time << endl;
      pout() << "max_step = " << a_max_step << endl;
    }  

  char buf[20];
  timeval time1,time2;
  
  bool calcExecTime = ((m_verbosity >= 1)|| ((m_verbosity <= -1)&&(procID() == 0)) );

  if (calcExecTime)
  {
  gettimeofday(&time1,NULL);
  gettimeofday(&time2,NULL);
  }

  int level;
  
  Real old_dt_base = m_dt_base;
  old_dt_base = 0.0; // should be zero to determine correct 'conclude' time if we calculate nothing

  for( ; (m_cur_step < a_max_step) && (m_cur_time < a_max_time);
       ++m_cur_step, m_cur_time += old_dt_base)
    {
      old_dt_base = m_dt_base;
      for(level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->time(m_cur_time);
        }

      if((m_checkpoint_interval > 0) &&(m_lastcheck_step != m_cur_step) &&
         (m_restart_step != m_cur_step) &&
         (m_cur_step % m_checkpoint_interval == 0))
        {
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }
      if((m_plot_interval > 0) &&
         (m_restart_step != m_cur_step) &&
         (m_cur_step % m_plot_interval == 0) )
        {
          writePlotFile();
        }
      if((m_tecplot_interval > 0) &&
         (m_cur_step % m_tecplot_interval == 0) )
        {
          writeTecPlotFile();
        }

      if((m_sliceInterval > 0) &&
         (m_cur_step % m_sliceInterval == 0) )
        {
          writeSlices();
        }
      if((m_lineInterval > 0) &&
         (m_cur_step % m_lineInterval == 0) )
        {
          writeLines();
        }
      if ((m_ilsReinitInterval > 0) && 
          (m_cur_step % m_ilsReinitInterval == 0)) 
        {
          reinitializeLevelSet();
        }
      if (needToWrite1DProbe(m_cur_time)) 
        {
          write1DProbe();
        }

      //static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->writeAxisXData(m_cur_step, m_cur_time);
//      static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->writeLevel0SF(m_cur_step);

    if (m_pSourCalc != NULL)
    {
      if (m_pSourCalc->checkTime(m_cur_time, m_cur_step))
      {
        /*if((m_checkpoint_interval > 0) &&(m_lastcheck_step != m_cur_step) &&
         (m_restart_step != m_cur_step) )
        {
          writeCheckpointFile();
          m_lastcheck_step= m_cur_step;
        }
        if ((m_plot_interval > 0) && (m_restart_step != m_cur_step))
        {
          writePlotFile();
        }
        if ((m_tecplot_interval > 0) &&(m_restart_step != m_cur_step)) 
        {
          writeTecPlotFile();
        }*/
        m_pSourCalc->CalculateSources(getAMRLevels(), m_cur_time, m_cur_step);
      }
    }

      ApplyProjectionScheme();
      m_steps_since_Poisson++;

      level = 0;
      int stepsLeft = 0;
      bool timeBoundary = true;            

      (void)timeStep(level,stepsLeft,timeBoundary);
      
      Real curtime = m_cur_time + old_dt_base;
      
      if (fabs(curtime-a_max_time) < numeric_limits<Real>::epsilon())
      {
        m_cur_time = a_max_time;
      } else      
      if (curtime + m_dt_new[0] > a_max_time)
      {
        m_dt_new[0] = a_max_time - curtime;
      } else      
      if (curtime + 2.0*m_dt_new[0] >= a_max_time)
      {
        m_dt_new[0] = 0.5*(a_max_time - curtime);        
      }

      assignDt();

      if (calcExecTime)
      {
        gettimeofday(&time2,NULL);
          double stepTime = (double)(time2.tv_sec-time1.tv_sec) + 1e-6*((double)(time2.tv_usec-time1.tv_usec));
                    
          pout() << "coarse time step " << setw(3) << m_cur_step
                 << "  old time = " << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << m_cur_time
                 << "  old dt = "   << setw(12)
                 << setprecision(6)
                 << setiosflags(ios::showpoint)
                 << setiosflags(ios::scientific)
                 << old_dt_base;

          //pout().flush();
          //pout() << setiosflags(ios::fixed)
          sprintf(buf,"%.4g",stepTime);
          pout()  << " (" << buf << " sec) "
                 //(" << time2.tv_sec << ", " << time2.tv_usec << ")-"
                 //<<                           "(" << time1.tv_sec << ", " << time1.tv_usec << ") "
                 << endl;
        time1 = time2;
    }
  }
  
  if ((m_cur_step == a_max_step) || (m_cur_time >= a_max_time))
  {
    m_cur_time -= old_dt_base;
    
    for(level = 0; level <= m_max_level; ++level)
        {
          m_amrlevels[level]->time(m_cur_time);
        }
    
    if (m_tecplot_interval > 0) 
      {
        writeTecPlotFile();
      }

    if (m_sliceInterval > 0)       
      {
        writeSlices();
      }
      
    if (m_lineInterval > 0)        
      {
        writeLines();
      }
      
    // write physics class per-level data
    for(int level = 0; level <= m_finest_level; ++level)
    {
      PhysProblem* pPhysPr = static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->getpatchMHDAM()->getPhysProblem();
      pPhysPr->conclude(static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->getStateNew(),m_cur_time,m_cur_step, m_inputFileContent);            
      
      //static_cast<AMRLevelIdealMHD*>(m_amrlevels[level])->V1V2TSOutput(m_cur_step);
    }
    
    //static_cast<AMRLevelIdealMHD*>(m_amrlevels[0])->wrtieLevel0H5();
        
    
    
  }

}

