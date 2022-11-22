#ifndef _KINETICSOURCES1D_H
#define _KINETICSOURCES1D_H

#include "KineticSources.H"

class KineticSources1D : public KineticSources
{
public:

  /// Constructor
  /**
   */
  KineticSources1D();

  /// Destructor
  /**
   */
  virtual ~KineticSources1D();

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
  
  
  // Calls external source calculator. Source terms are returned in a_ST
  virtual void CallExternalSC(const FArrayBox& Level0Data, Real a_dt, Real a_dx, FArrayBox& a_SourceTerms);
  
  
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

  
  
protected:
  int m_iTotalNeutrals;  
  Real m_dKineticRuntimeYears;        // in years
  Real m_dTimeFactor;                 // Time in years = m_dTimeFactor * TimeMHDAM
  
  Real m_dKineticTimeIntervalYears;   // in years
  Real m_dKineticTimeIntervalMHDAM;   // in dimensionless units
  Real m_dNextKineticTimeMHDAM;   
  
  int m_iKineticStepInterval;       
  int m_iKineticRuntimeConsistent;
        
  std::string m_sKineticRestartFile;
  
  Real  m_dTemperatureL;  
  Real  m_dNetnumL;   
  Real  m_dVelxL;
  
  Real m_dExtraLength;
  int  m_iNumExtraCells;
    
  
private:
  int m_cur_step;
  

  // Disallowed for all the usual reasons
  void operator=(const KineticSources1D& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  KineticSources1D(const KineticSources1D& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
