#ifndef _KINETICSOURCES3D_H
#define _KINETICSOURCES3D_H

#include "SourceCalculator.H"

#include <queue>

class KineticSources3D : public SourceCalculator
{
public:

  // This should be bigger class similar to Box
  struct RealBox
  {
    RealVect m_smallEnd;
    RealVect m_bigEnd;    
  };

  /// Constructor
  /**
   */
  KineticSources3D();

  /// Destructor
  /**
   */
  virtual ~KineticSources3D();

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
  virtual int numVarsForSC( void ) ;
  
  /// Number of calculated source terms
  /**
   */
  virtual int numSourceTerms( void );
      

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
  Real m_totalNeutrals;         // it will be rounded to 64-bit integer in fortran code
  
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
  RealVect  m_sunXYZ;
  IntVect   m_sunIJK;  
  Real  m_lismN;
  Real  m_lismV;
  Real  m_lismUX;
  Real  m_lismUY;
  Real  m_lismUZ;
  Real  m_lismT;
  Real  m_netN;
  // Node coordinates of the Sun on base level (lower corner of corresponding cell)
  // This is (1,1) coordinates in Jacob's code on all levels.
  
  
  Vector < RealBox >  m_gridsBoundaries; // in AU, in my notation (Sun is not at (0,0,0))
  Vector < Box >      m_kineticBoxes;       // Box in MHDAM mesh that coincides with kinetic mesh
  Vector <LevelData<FArrayBox>* > m_plasmaDataLD;
  Vector <FArrayBox*>     m_plasmaDataFAB;
  Vector <FArrayBox*>     m_sourceTerms;
  Vector <FArrayBox*>     m_neutralData;  
  
    
  int m_ngrids;    
    
  bool m_firstCall;
    
  
  // Contains file names of last written m_maxChkFiles kinetic checkpoints
  std::queue<std::string> m_lastKineticFiles; // Don't forget to remove #include <queue>
  
  std::string m_chk_prefix;
      
private:
  int m_cur_step;
  

  // Disallowed for all the usual reasons
  void operator=(const KineticSources3D& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  KineticSources3D(const KineticSources3D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
