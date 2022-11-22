#include <sys/stat.h>

#include "ProbeFilesInfo.H"
#include "AMRLevelIdealMHD.H"
#include "MHDAMDefs.H"

ProbeFilesInfo::ProbeInfo::ProbeInfo(int a_numComp, std::string & a_fileName)
{
  int errcode,eclass,resultlen; char errstring[MPI_MAX_ERROR_STRING]; 
  
  m_fileName = a_fileName;
  
#ifdef CH_MPI      
  errcode = MPI_File_open(Chombo_MPI::comm, const_cast<char*>(m_fileName.c_str()), MPI_MODE_WRONLY | MPI_MODE_CREATE | MPI_MODE_APPEND, MPI_INFO_NULL, &m_file);
  if (errcode!=MPI_SUCCESS)
  {
    MPI_Error_class(errcode, &eclass);
    MPI_Error_string(errcode, errstring, &resultlen);
    if (procID() == 0) pout() << "ProbeFilesInfo::ProbeInfo::ProbeInfo, MPI_File_open, Error " <<  eclass << " " << errstring << endl;
    MPI_Abort(Chombo_MPI::comm, -1);
  }
  MPI_File_set_view(m_file, MPI_DISPLACEMENT_CURRENT, MPI_CHAR, MPI_CHAR, "external32", MPI_INFO_NULL);    
#else
  m_file = fopen(m_fileName.c_str(),"a");    
#endif
}

ProbeFilesInfo::ProbeInfo::~ProbeInfo()
{  
#ifdef CH_MPI      
  MPI_File_close(&m_file);  
#else
  if (m_file!=NULL) fclose(m_file);    
#endif
}

void ProbeFilesInfo::ProbeInfo::writeString(char * a_str) const
{
#ifdef CH_MPI          
    MPI_File_write_shared(m_file, a_str, strlen(a_str), MPI_CHAR, MPI_STATUS_IGNORE);
#else      
    if (m_file!=NULL) fprintf(m_file,"%s",a_str);      
#endif      
}


ProbeFilesInfo::ProbeInfoPoint::ProbeInfoPoint(int a_numComp, std::string & a_fileName, RealVect & a_uvw)
 : ProbeInfo(a_numComp, a_fileName), m_uvw(a_uvw)
{
}

ProbeFilesInfo::ProbeInfoTrajectory::ProbeInfoTrajectory(int a_numComp, std::string & a_fileName, std::string & a_trajFile)
  : ProbeInfo(a_numComp, a_fileName)
{
  Real time,lat,lon,vecHGI[3],vec[3],r; 
    
  FILE* V_file=fopen(a_trajFile.c_str(),"r");  
  if (V_file==NULL) 
  {
    std::string buf = a_trajFile + std::string(" not found\n");
    MayDay::Warning(buf.c_str());
  }
  
  if (procID() == 0)  pout() << "reading " << a_trajFile << endl;
  
  while (!feof(V_file))
  {    
    
    Real vdata[5];   int year; 
    fscanf(V_file,"%i %lf %lf %lf %lf", &year, &vdata[1], &vdata[2], &vdata[3], &vdata[4]);           
           
    int daysYear = (year%4 == 0 ? 366 : 365);
    time = year + (vdata[1]-1)/daysYear; 
    
    //if (time<m_startBCPhys) continue;
        
    
    r    = vdata[2];
    lat  = vdata[3]/180.0*d_PI;
    lon  = vdata[4]/180.0*d_PI;
           
    vecHGI[0] = cos(lat)*cos(lon);
    vecHGI[1] = cos(lat)*sin(lon);
    vecHGI[2] = sin(lat);
    
    vec[0]    = -0.9998430685*vecHGI[0]  + 0.01771548383*vecHGI[1];
    vec[1]    = -0.01771548383*vecHGI[0] - 0.9998430685*vecHGI[1];
    vec[2]    = vecHGI[2];    
    
    lon = acos(vec[0]/(sqrt(vec[0]*vec[0]+vec[1]*vec[1])));   
    if (vec[1]<0.0) lon = d_2PI-lon;
                           
    m_pos[time][0] = r;
    m_pos[time][1] = lon;
    m_pos[time][2] = d_PI_2-lat;
  }
  fclose(V_file);  
}


void ProbeFilesInfo::ProbeInfoTrajectory::get1Dposition(RealVect & a_uvw, Real a_time) const
{
  std::map<Real,RealVect>::const_iterator vpos_iter = m_pos.upper_bound(a_time);
  a_uvw[0] = vpos_iter->second[0];
  a_uvw[1] = vpos_iter->second[1];
  a_uvw[2] = vpos_iter->second[2];             
}

ProbeFilesInfo::ProbeFilesInfo()
{
  m_nProbes = 0;
  m_probes  = NULL;
}

ProbeFilesInfo::~ProbeFilesInfo()
{
  for (int i=0;i<m_nProbes;i++) delete m_probes[i];
  
  delete[] m_probes;  
}

void ProbeFilesInfo::setupProbeFiles(int a_num_vars, std::string& a_header)
{
  enum eProbeType {
                        PT_Undefined   = 0,
                        PT_1Dpoint     = 1,
                        PT_Trajectory  = 2
                        };
                        
  ParmParse parser("mhdam");      
  
  int i;
  int num_probes = 0;  
  parser.query("num_probes", num_probes);           
  
  if (num_probes <= 0) return;
  
  m_nProbes = num_probes;
  m_probes  = new ProbeInfo*[m_nProbes];
  
  
  char buf[20];
  
  for (i = 0;i < num_probes; i++)
  {
    sprintf(buf,"probe%i",i);
    ParmParse probe_parser(buf);      
    
    int type = PT_Undefined;
    probe_parser.query("type", type);     

    std::string outputFile;
    probe_parser.query("filename", outputFile);     
    
    if (procID() == 0)  
    {      
      struct stat stFileInfo;  
      int intStat; 
      
      pout() << "creating probe file " << outputFile << std::endl;
        
      FILE* plot_file;

      // Attempt to get the file attributes 
      intStat   = stat(outputFile.c_str(),&stFileInfo); 
      plot_file = fopen(outputFile.c_str(),"a");
      if(intStat == 0)
      {             
        fprintf(plot_file,"\n===========\n");
      } else
      { 
        fprintf(plot_file,"%s\n",a_header.c_str());      
      }    
      fclose(plot_file);
    }
    
#ifdef CH_MPI  
    MPI_Barrier(Chombo_MPI::comm);
#endif
    
    if (type == PT_1Dpoint)
    {
      Vector<Real> coords;
      probe_parser.queryarr("coords", coords, 0, SpaceDim);
      
      RealVect rv(coords);
      m_probes[i] = new ProbeInfoPoint(a_num_vars, outputFile, rv);
    } else
    if (type == PT_Trajectory)
    {
      std::string coordsFile;
      probe_parser.query("coords", coordsFile);           
            
      m_probes[i] = new ProbeInfoTrajectory(a_num_vars, outputFile, coordsFile);
    } else MayDay::Error("Unknown probe type");   
  }
    
}
