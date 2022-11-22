#ifndef _KINETICSOURCES_H
#define _KINETICSOURCES_H

#include "SourceCalculator.H"

class AMRLevelIdealMHD;

class KineticSources : public SourceCalculator
{
public:

  /// Constructor
  /**
   */
  KineticSources();

  /// Destructor
  /**
   */
  virtual ~KineticSources();

  /// Factory method - this object is its own factory
  /**
     Return a point to a new KineticSources object with m_isDefined = false (i.e.,
     its define() must be called before it is used).
   */
  virtual SourceCalculator* new_SourceCalculator( void );

  // number of variables that are needed for external source calculator
  virtual int numVarsForSC() = 0;
  
  // Using conservative variables a_U, the method prepares data (stored in a_Buf)
  // that will be used in an external source calculator.
  virtual void PrepareDataForSC(const FArrayBox& a_U, FArrayBox& a_Buf) = 0;
  
  
  // Calls external source calculator. Source terms are returned in a_ST
  virtual void CallExternalSC(const FArrayBox& Level0Data, Real a_dt, Real a_dx, FArrayBox& a_SourceTerms) = 0;  
    
  
  // Modification of source terms calculated by external source calculator.
  virtual void ModifySourceTerms(FArrayBox& a_SourceTerms, const Box& a_box) = 0;
      
                                  
  /// Calculate source terms
  /**
    During one time step this method must be called only once on each processor
   */           
  virtual void CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep);
  
protected:
  
  // Send data from Level "0" to zero processor
  void sendDataToZeroProcessor(const AMRLevelIdealMHD* aLevel0);
  
  // Gather data on zero processor. Gathered data are stored in Level0Data
  void GatherDataOnZeroProcessor(const AMRLevelIdealMHD* aLevel0, FArrayBox& Level0Data);
  
  // Broadcast source terms after they have been computed
  void BroadcastSourceTerms(AMRLevelIdealMHD* aLevel0, FArrayBox& a_SourceTerms);
  
  /// Interpolate source terms from zero level to higher levels.
  /**
      It is usually called from CalculateSources. This method fills m_SCData field
      of AMRLevelIdealMHD class.
   */
  virtual void InterpolateSourceTerms(Vector<AMRLevel*> a_levels);
    
  
  
private:


  // Disallowed for all the usual reasons
  void operator=(const KineticSources& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  KineticSources(const KineticSources& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
