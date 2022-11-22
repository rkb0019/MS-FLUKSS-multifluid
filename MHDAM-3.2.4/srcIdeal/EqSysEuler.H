#ifndef _EQSYSEULER_H_
#define _EQSYSEULER_H_

#include "EquationSystem.H"

///
/**
   The base class EquationSystem descrbes the system of hyperbolic equations 
   solved in the package
 */
  
class EqSysEuler : public EquationSystem
{
public:

  /// Constructor
  /**
   */
  EqSysEuler(int a_ts);

  /// Destructor
  /**
   */
  virtual ~EqSysEuler();
  
  
  /// Number of state variables
  /**
     Return the number of state variables.
   */
  virtual int numStates();
  
  /// Names of the conserved variables
  /**
     Return the names of the conserved variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numConserved()-1.
   */
  virtual Vector<string> stateNames();
  
  /// Number of primitive variables
  /**
     Return the number of primitive variables.  This may be greater than the
     number of conserved variables if derived/redundant quantities are also
     stored for convenience.
   */
  virtual int numPrimitives();
  
  /// Names of the primitive variables
  /**
     Return the names of the primitive variables.  A default implementation is
     provided that puts in generic names (i.e., "variable#" which "#" ranges
     for 0 to numPrimitive()-1.
   */
  virtual Vector<string> primitiveNames();

  /// Number of flux variables
  /**
     Return the  number of flux variables.  This can be greater than the number
     of conserved variables if addition fluxes/face-centered quantities are
     computed.
   */
  virtual int numFluxes();

  /// Number of primitive variables for which slopes are computed
  /**
     Return the number of primitive variables for which slopes are computed.
     Only slopes corresponding to primitive variables in the interval 0 to
     numSlopes() - 1 are computed and only primitive variables in that interval
     are updated using the slopes.
   */
  virtual int numSlopes();


  /// Compute the primitive variables from the conserved variables within a_box
  /**
   */
  virtual void stateToPrim(      FArrayBox& a_W,
                          const FArrayBox& a_U,
                          const Box&       a_box);
                          
  /// Compute the conserved variables from the primitive variables within a_box
  /**
   */
  virtual void primToState(      FArrayBox& a_U,
                          const FArrayBox& a_W,
                          const Box&       a_box);


  /// Transform a_dWLeft and a_dWRight from primitive to characteristic variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the primitive variables.
     On output, it contains the increments in the characteristic variables.
  */
  virtual void charAnalysis(       FArrayBox & a_dWLeft,
                                   FArrayBox & a_dWRight,
                             const FArrayBox & a_W,
                             const int &       a_dir,
                             const Box &       a_box);

  /// Transform a_dWLeft and a_dWRight from characteristic to primitive variables
  /**
     On input, a_dWLeft and a_dWRight contain the increments of the characteristic variables.
     On output, it contains the increments in the primitive variables.
  */
  virtual void charSynthesis(       FArrayBox & a_dWLeft,
                                    FArrayBox & a_dWRight,
                              const FArrayBox & a_W,
                              const int &       a_dir,
                              const Box &       a_box);
                              
  virtual void vectorVars(Vector<int> & a_vars) const;
  
  /// Returns first components of velocities in the array of conservative/primitive variables
  /**
   */
  virtual void velocityVars(Vector<int> & a_vars) const;

protected:
    
  
private:
  // Disallowed for all the usual reasons
  void operator=(const EqSysEuler& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  EqSysEuler(const EqSysEuler& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
