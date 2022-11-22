#ifndef _KINETICSOURCES2D_H
#define _KINETICSOURCES2D_H

#include "SourceCalculator.H"

#include <queue>

class KineticSources2D : public SourceCalculator
{
public:

  /// Constructor
  /**
   */
  KineticSources2D();

  /// Destructor
  /**
   */
  virtual ~KineticSources2D();

  /// Factory method - this object is its own factory
  /**
     Return a point to a new KineticSources object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual SourceCalculator* new_SourceCalculator( void );

  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );
  
  // External source calculator may require some initialization before it 
  // can be used. Put initialization statements here.
  //
  // This method is called only once before time loop starts.
  virtual void initExternalSC(Vector<AMRLevel*> a_levels);
  
  // number of variables that are needed for external source calculator
  virtual int numVarsForSC() ;
  
  /// Number of calculated source terms
  /**
   */
  virtual int numSourceTerms( void );
  
  // Using conservative variables a_U, the method prepares data (stored in a_Buf)
  // that will be used in an external source calculator.
  virtual void PrepareDataForSC(const FArrayBox& a_U, FArrayBox& a_Buf);
    
    
  // Modification of source terms calculated by external source calculator.
  virtual void ModifySourceTerms(FArrayBox& a_SourceTerms, const Box& a_box);
  
      

  /// Check time for source terms calculation
  /**
    dTime is current time, iStep - current step
    Returns "true" if source terms must be calculated
   */
  virtual bool checkTime( Real dTime, int iStep );

  /// Add external source terms
  /**
   */       
  virtual void addExternalSources(       FArrayBox & a_U,
                                   const FArrayBox & a_S,
                                   const FArrayBox & a_W,
                                   const Real      & a_dt,
                                   const Real      & a_dx,
                                   const Box       & a_box);
                                   
                                   
  /// Write information for restarting source calculator
  virtual void writeCheckpointFile(const char* chkFileName) const;

  virtual void CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep);

  /// Number additional variables for writing to plot file
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */                                
  virtual int numPlotVars();
  
  /// Names of the additional variables for writing to plot file
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual Vector<std::string> plotNames();  
  
  
protected:
  Real m_totalNeutrals;         // it will be rounded to 64-bit integer
  
  Real m_runtimeYears;          // Kinetic run time (in years)
  Real m_initialRuntimeYears;
  Real m_timeFactor;            // Time in years = m_timeFactor * TimeMHDAM  
  Real m_timeIntervalYears;     // in years
  Real m_timeIntervalMHDAM;     // in dimensionless units    
  int  m_stepInterval;       
  int  m_photoionize;
        
  std::string m_restartFile;   
  Real m_nextTimeMHDAM;   
  
  // SW-Lism parameters
  Real  m_sunXC;
  Real  m_sunYC;
  Real  m_sunZC;
  Real  m_lismN;
  Real  m_lismV;
  Real  m_lismT;
  Real  m_netN;
  // Node coordinates of the Sun on base level (lower corner of corresponding cell)
  // This is (1,1) coordinates in Jacob's code on all levels.
  int m_sunIC;  
  int m_sunJC;  
  int m_sunKC;  
  
  int m_ngrids;  
  Vector < IntVect > m_gridsDim;        // Kinetic grid dimensions (stores Jacob's nz,nx)
  Vector < RealVect> m_gridsBoundaries; // in AU
  Vector < Box >     m_gridsBoxes;      // Box in MHDAM mesh that coincides with kinetic mesh
  
  bool m_firstCall;
  
  
  // Send data to zero processor
  void sendDataToZeroProc(Vector<AMRLevel*> a_levels);
  
  // Gather data on zero processor. Gathered data are stored in a_Data
  void GatherDataOnZeroProc(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_Data);
  
  // Broadcast source terms after they have been computed
  void BroadcastSourceTerms(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_Data);
    
  // Fill field AMRLevelIdealMHD.m_SCData on all levels with correct data
  void PrepareSourceTermsOnLevels(Vector<AMRLevel*> a_levels, Vector<FArrayBox*>& a_SCData, Vector< BaseFab<int>* > Region);

  // Contains file names of last written m_maxChkFiles kinetic checkpoints
  std::queue<std::string> m_lastKineticFiles; // Don't forget to remove #include <queue>
  
  std::string m_chk_prefix;
      
private:
  int m_cur_step;
  

  // Disallowed for all the usual reasons
  void operator=(const KineticSources2D& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  KineticSources2D(const KineticSources2D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
