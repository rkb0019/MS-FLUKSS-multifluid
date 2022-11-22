#ifndef _SOURCECALCULATOR_H
#define _SOURCECALCULATOR_H

#include "FArrayBox.H"
#include "ParmParse.H"
#include "REAL.H"
#include "LevelData.H"
#include "Vector.H"
#include "AMRLevel.H"

#include "PhysProblem.H"

class SourceCalculator
{
public:

  /// Constructor
  /**
   */
  SourceCalculator();

  /// Destructor
  /**
   */
  virtual ~SourceCalculator();

  /// Factory method - this object is its own factory
  /**
     Return a point to a new SourceCalculator object with m_isDefined = false (i.e.,
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
  virtual void initExternalSC(Vector<AMRLevel*>  a_levels);
  
  // Initial values of source terms 
  // By default, sets the values to zero.
  virtual void initialize(LevelData<FArrayBox>& a_S);

  /// Number of calculated source terms
  /**
   */
  virtual int numSourceTerms( void ){ return 0; }
  
  /// Check time for source terms calculation
  /**
    dTime is current time, iStep - current step
    Returns "true" if source terms must be calculated
   */
  virtual bool checkTime(Real dTime, int iStep);


  /// Add external source terms
  /**
   */
  virtual void addExternalSources(       FArrayBox & a_U,
                                   const FArrayBox & a_S,
                                   const FArrayBox & a_W,
                                   const Real      & a_dt,
                                   const Real      & a_dx,
                                   const Box       & a_box);
                                   
  /// Calculate source terms
  /**
    During one time step this method must be called only once on each processor
    The main method of the class.
   */           
  virtual void CalculateSources(Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep);
  
  // Write information for restarting source calculator
  /**
   */
  virtual void writeCheckpointFile(const char* chkFileName) const;  
  
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

  int m_verbosity;
  
  
  

private:
  // Disallowed for all the usual reasons
  void operator=(const SourceCalculator& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  SourceCalculator(const SourceCalculator& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
