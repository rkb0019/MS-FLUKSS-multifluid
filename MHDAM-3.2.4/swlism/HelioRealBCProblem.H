#ifndef  HelioRealBCProblem_H_
#define  HelioRealBCProblem_H_

#include "HeliosphericProblem.H"

/// Base class for 3D heliospheric problem
/**
   Parameters:
 */
class HelioRealBCProblem : public HeliosphericProblem
{  
public:
  struct GMIRData
  {
    int nPoints;
    Real * time;
    Real * N;
    Real * V;
    Real * B;
    Real * T;    
  } m_GMIRData;
  
  struct ThetaDistrib
  {
    Real time;        
    Real * N;
    Real * V;    
    Real * T;
    ThetaDistrib(){N=NULL;V=NULL;T=NULL;}
    ~ThetaDistrib()
    {
      if (N!=NULL) delete[] N;
      if (V!=NULL) delete[] V;
      if (T!=NULL) delete[] T;      
    }
        
  } * m_ThetaDistrib;
  
  int m_nThetaDistrib;
  
  struct OMNIData
  {
    int nPoints;
    Real * time;
    Real * N;
    Real * V;
    Real * Br;
    Real * Bp;
    Real * Bt;
    Real * B;
    Real * Tilt;
    Real * T;
    Real * Z2;
  } m_OMNIData;
  
  
  struct WSOData
  {
    int nPoints;
    Real * time;    
    Real * Tilt;    
    WSOData(){nPoints=0;time=NULL;Tilt=NULL;}
  } m_WSOData;
  
  //std::map<Real,RealVect> m_V1pos;
  //std::map<Real,RealVect> m_V2pos;
  //std::map<Real,RealVect> m_Ulyssespos;
  //std::map<Real,RealVect> m_Earthpos;
  
  
  
  /// Null constructor
  /**
   */
  HelioRealBCProblem();

  /// Destructor
  /**
   */
  ~HelioRealBCProblem();
  


  /// Factory method - this object is its own factory
  /**
     Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
     its define() must be called before it is used) and m_isFortranCommonSet
     set to value of m_isFortranCommonset in the current (factory) object.
   */
  virtual PhysProblem *new_PhysProblem();
  
  /// Copy method 
  /**
     Copy all data from a_PP to this instance.
   */
  virtual void copy_PhysProblem(const PhysProblem* a_PP);
  
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );
  
  /// Creates "DATASETAUXDATA" fields for tecplot
  /**
  */
  virtual void auxDataTecplot(std::string & a_str,
                              Real          a_time,
                              int           a_type); 


  /// Fill ghost cells
  /**
   */
  virtual void fillGhostCells(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );

  /// Set boundary fluxes
  /**
   */
  virtual void fluxBC(       FArrayBox&      a_F,
                             FArrayBox&      a_Bn,
                       const FArrayBox&      a_W,
                       const FArrayBox&      a_Wextrap,
                       const int&            a_dir,
                       const Side::LoHiSide& a_side,
                       const Real&           a_time);
                       
  virtual void define(const ProblemDomain& a_domain,                                  
                      const int            a_level);


  /// Set up initial conditions
  /**
   */
  virtual void initialize(LevelData<FArrayBox>& a_U);
  
  /// Problem specific postprocessing
  /**
   */
  virtual void postprocessing(       FArrayBox & a_U,
                               const FArrayBox & a_W,
                               const Real      & a_dt,
                               const Real      & a_time,                               
                               const Box       & a_box);

  

  /// Creates tagged cells for dynamic mesh refinement
  /**
    Problem specific cells tagging
   */
  virtual void tagCells(const FArrayBox&  a_U,
                        const Box&        a_box,
                              IntVectSet& a_tags);
                              
  void tagCellsCIR(const FArrayBox&  a_U,
                    const Box&        a_box,
                          IntVectSet& a_tags);
                          
  void tagCellsCME(const FArrayBox&  a_U,
                    const Box&        a_box,
                          IntVectSet& a_tags);

  /// Check geometrical/problem limitations for grid adaptation
  /**
   */  
  virtual void lockedCellsRegrid( BaseFab<int> & a_flag,
                            const FArrayBox&  a_U,
                            const Box&     a_box);

  /// Return boundary condition flags for all boundaries
  /**
   */
  virtual void getBCFlags( eBoundaryConditions leftBC,
                           eBoundaryConditions rightBC,
                           eBoundaryConditions bottomBC,
                           eBoundaryConditions topBC,
                           eBoundaryConditions frontBC,
                           eBoundaryConditions behindBC );

  /// Return the boundary condition flag for given boundary
  /**
   */
  virtual eBoundaryConditions getBCFlags( int a_dir, Side::LoHiSide a_sd );
  
                                
  
  /// Modifies curvilinear mesh (applies to curvilinear meshes only)
  /**
   */                               
  virtual void defineMesh(const ProblemDomain & a_prob_domain,
                          const Vector<Real>  & a_domainBox);
  
  /// For a tracking surface 'a_s' returns first component of the velocity field for the level set method
  /**
   */                               
  virtual int lsIndexField(int a_s);
  
  /// Converts dimensionless time to problem specific time
  /**
   */                               
  virtual Real getPhysTime(Real a_time); 
  
  /// Converts 'probe_dt' value in inputs file to dimensionless units
  /**
   */
  virtual Real probeFilesDimensiolessDt(Real a_dt);
  
    
 
  /// Returns curvilinear coordinates for 1D probes
  /**
   */                        
  virtual void probeFilesHeader(int & a_numVars, std::string& a_header);

                             
  //void readVoyagersPos();
  //void readUlyssesPos();
  //void readEarthPos();

protected:
                                  
  FArrayBox m_OdstrcilData; // Odstrcil's initial data  
  void readOdstrcilData();
  
  
  Real m_startBC;     // Dimemsionless time when RealBC should be started (usually time in a checkpoint the probelm is started from)
  Real m_startBCPhys; // Physical time when BC are initiated
       
  void readGMIRData();  
  void readThetaDistrib();
  
  std::string m_OMNIDataFile;
  void readOMNIData();
  void readOMNIDataTurb();
  
  void readWSOData();

  // Data for HPBD_OMNIPROBLEM
  bool m_iOldOMNIFile;
  
  // Data for HPBC_CIRPROBLEM and HPBC_HDF5INPUT problems
  Real m_dataR0;
  std::string m_h5File;
  hid_t       m_h5fd;
  int m_nDataSets;  
  int m_iData1;
  int m_iData2;
  FArrayBox   m_h5Data1;
  FArrayBox   m_h5Data2;
  Real m_time1;
  Real m_time2;
  Real m_firstDataSetTime;
  Real m_lastDataSetTime;
  bool m_bDataRotatingFrame;
  
  void readCIRData();
  void readH5DataSet(hid_t a_h5f, int a_iDataSet, FArrayBox & a_data, Real & a_time);
  void prepareHDF5InputProblem();
    
  
  
  void fillGhostCellsOdstrcil(       FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsGMIR(           FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsSuess(          FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsOMNI(           FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsWSO(            FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsCIR(            FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
  void fillGhostCellsH5Input(        FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time );
                               
                               
                               
                                

private:
                                         // Disallowed for all the usual reasons
  void operator=(const HelioRealBCProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
                                         // Disallowed for all the usual reasons
  HelioRealBCProblem(const HelioRealBCProblem& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif /*  HelioRealBCProblem_H_ */


