#ifndef _PROBEFILESINFO_H_
#define _PROBEFILESINFO_H_

#include "Vector.H"
#include "CH_HDF5.H"
#include "RealVect.H"


class ProbeFilesInfo
{  
  public:
  
  // Not classical object orientd approach
  class ProbeInfo
  {    
    public:    
      ProbeInfo(int a_numComp, std::string & a_fileName);
      virtual ~ProbeInfo();
            
      virtual void get1Dposition(RealVect & a_uvw, Real a_time) const = 0;
      
      virtual void writeString(char * a_str) const;
      
    protected:

      std::string m_fileName;
  #ifdef CH_MPI  
      mutable MPI_File m_file;
  #else
      mutable FILE*    m_file;
  #endif
      
      // Reserved for future use
      //void incrementBuf(AMRLevelIdealMHD & a_level0, Real a_time);  
      //Vector<Real> m_buf;      
  };
  
  class ProbeInfoPoint : public ProbeInfo
  {
    public:
      ProbeInfoPoint(int a_numComp, std::string & a_fileName, RealVect & a_uvw);
      virtual void get1Dposition(RealVect & a_uvw, Real a_time) const {a_uvw = m_uvw;}      
      
    protected:
      RealVect  m_uvw;        
  };
  
  class ProbeInfoTrajectory : public ProbeInfo
  {
    public:
      ProbeInfoTrajectory(int a_numComp, std::string & a_fileName, std::string & a_trajFile);
      virtual void get1Dposition(RealVect & a_uvw, Real a_time) const;
      
    protected:    
      std::map<Real,RealVect> m_pos;
    
  };
    
    ProbeFilesInfo();
    virtual ~ProbeFilesInfo();
    
    inline int numProbes() {return m_nProbes;}
    
    const ProbeInfo & getProbe(int a_iProbe)
    {
      CH_assert((a_iProbe>=0) && (a_iProbe<m_nProbes));
      return *(m_probes[a_iProbe]);
    }
    
    void setupProbeFiles(int a_num_vars, std::string& a_header);
    
  protected:
  
    int m_nProbes;
    ProbeInfo ** m_probes;
};


#endif