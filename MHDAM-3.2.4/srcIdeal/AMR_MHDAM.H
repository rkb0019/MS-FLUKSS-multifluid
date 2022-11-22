#ifndef _AMR_MHDAM_H_
#define _AMR_MHDAM_H_

#include <queue>

#include "AMR.H"
#include "RealVect.H"
#include "ProbeFilesInfo.H"


class SourceCalculator;


class AMR_MHDAM : public AMR
{
public:
  ///
  /**
     Default constructor.  Use must subsequently call a define() function
     and a setUp function.
  */
  AMR_MHDAM();

  ///
  /**
     Destructor.
  */
  virtual ~AMR_MHDAM();

  ///
  /**
     Defines this object.   User must call a setup function
     before running.
     Arguments:

     a_prob_domain (not modified): problem domain on the base level.

     a_ref_ratios (not modified): refinement ratios.  There
     must be at least a_max_level+1 elements,  or an error will
     result. Element   zero is the base level.

     a_max_level (not modified): the maximum level allowed, where the
     base level is zero.  There will be a total of a_max_level+1
     levels, since level zero and level a_max_level will both
     exist.

     a_amrlevel_ptr (not modified): pointer to a physics class factory
     object.
     The object it points to is used to construct the collection of
     AMRLevels in this AMR as objects of the physics class
     type.  It can be destructed after this define()
     returns.
  */
  void define(int a_max_level,
              const Vector<int>& a_ref_ratios,
              const Box& a_prob_domain,
              const AMRLevelFactory* const a_amrLevelFact,
              SourceCalculator* a_SourCalc = NULL);


  ///
  /**
     Defines this object.   User must call a setup function
     before running.

     a_prob_domain (not modified): problem domain on the base level.

     a_ref_ratios (not modified): refinement ratios.  There
     must be at least a_max_level+1 elements,  or an error will
     result. Element   zero is the base level.

     a_max_level (not modified): the maximum level allowed, where the
     base level is zero.  There will be a total of a_max_level+1
     levels, since level zero and level a_max_level will both
     exist.

     a_amrlevel_ptr (not modified): pointer to a physics class factory
     object.
     The object it points to is used to construct the collection of
     AMRLevels in this AMR as objects of the physics class
     type.  It can be destructed after this define()
     returns.
  */
  void define(int a_max_level,
              const Vector<int>& a_ref_ratios,
              const ProblemDomain& a_prob_domain,
              const AMRLevelFactory* const a_amrLevelFact,
              SourceCalculator* a_SourCalc = NULL);

#ifdef CH_USE_HDF5
  /// Setup for restart
  /**
   */
  void setupForRestart(HDF5Handle& a_handle,
                       int  a_fixedGridSize);
#endif

  /// Go baby, go
  /**
   */
  virtual void run(Real a_max_time, int a_max_step);

  /// ????
  /**
   */
  int getCurrentStep() const;
      

  /// Clear divB by applying projection scheme (solving Poisson equation)
  /**
   */
  void ApplyProjectionScheme();

  /// set input file
  /**
   */
  void setInputFile(const char* a_inFile);

  /// set slice parameters
  /**
   */
  void setSliceParameters(int a_sliceIntreval, int a_numSlices,
              std::vector<std::string>  m_slicePlanes,
              std::vector<Real> m_sliceValues);
              
  void setStepFactor(const Vector<int>& step_factor);
  
              
  // computes maximum stable time step given the maximum stable time
  // step on the individual levels.
  virtual void assignDt();
  
  /// Input parameters
  /**
   */
  void input( ParmParse & a_parser, int a_verbosity );
  
  // Creates probe files and allocates m_pfi
  void setup1DProbeFiles();


protected:


  // Source calculator
  SourceCalculator* m_pSourCalc;

  // Tecplot output interval
  int  m_tecplot_interval;

  // Leave only last m_maxChkFiles checkpoints
  int m_maxChkFiles;
  
  // True if map files should be written for checkpoints
  bool m_outputChkMapFiles;

  // Contains file names of last written m_maxChkFiles checkpoints
  // "mutable" means that we can change this field in const methods
  mutable std::queue<std::string> m_lastChkFiles;

  // Contains file names of last written m_maxChkFiles kinetic checkpoints
  mutable std::queue<std::string> m_lastKineticFiles;

  // Store here an input file name
  std::string m_inputFileName;

  // Store here an input file content to write it to checkpoint
  std::string m_inputFileContent;

  // Slice parameters
  int m_sliceInterval;
  int m_numSlices;
  std::vector<int>  m_slicePlanes;
  std::vector<Real> m_sliceValues;
  
  int m_numLines;
  int m_lineInterval;
  std::vector<RealVect> m_lineBgn;
  std::vector<RealVect> m_lineEnd;
  
  // LevelSet renitialization interval
  int m_ilsReinitInterval;

  // How often (in time steps) Poisson solver is used
  int m_PoissonInterval;
  int m_steps_since_Poisson;
  
  // Determines subcycling in time on finer levels
  Vector<int> m_step_factor;
  
  // Subcycling in time on finer levels as it appears in input file.
  Vector<int> m_input_step_factor;
  
  // File info for probe files
  ProbeFilesInfo m_pfi;
  
  Real m_probeDt;    // How often output probes (in dimensionless)    
  Real m_last1dprobe;// Time when last probe was written
    
  
  // advance by dt on this level and all finer levels and return the
  // number of steps left - given the number of steps left on entry.
  virtual int timeStep(int a_level, int a_stepsLeft, bool a_coarseTimeBoundary);
    
  // Write plot file in hdf5 format
  virtual void writePlotFile() const;

  // Write checkpoint plot file
  virtual void writeCheckpointFile()  const;
  
  // Write map file for cruvilinear meshes 
  void writeMappedGeometry(const std::string& a_fileRoot, IntVect a_outputGhost = IntVect::Zero) const;

  /// Writes data to tecplot format
  void writeTecPlotFile();

  /// writes slices to disk
  void writeSlices();

  /// Creates slice
  /**
    plane = 'X', x = sValue
    plane = 'Y', y = sValue
    plane = 'Z', z = sValue
   */
  void writeOneSlice(char plane, Real sValue) const;
  
  /// extract data along lines
  void writeLines();
  
  /// Extract data along one line
  void writeOneLine(int iLine, RealVect a_Bgn, RealVect a_End) const;
    
  void write1DProbe();
    
  
  void reinitializeLevelSet();
  
  /// Returns true if 1D probe should be written
  /**
   */
  bool needToWrite1DProbe(Real a_time);

};
#endif

