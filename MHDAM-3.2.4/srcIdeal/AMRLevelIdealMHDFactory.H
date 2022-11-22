#ifndef _AMR_LEVEL_IDEAL_MHD_FACTORY_H_
#define _AMR_LEVEL_IDEAL_MHD_FACTORY_H_

#include "AMRLevelFactory.H"
#include "AMRLevelIdealMHD.H"

/// AMR MHDAM factory
/**
 */
class AMRLevelIdealMHDFactory : public AMRLevelFactory
{
  
public:
  /// Null constructor
  /**
   */
  AMRLevelIdealMHDFactory();

  /// Virtual constructor
  /**
   */
  virtual AMRLevel* new_amrlevel() const;

  /// Destructor
  /**
   */
  virtual ~AMRLevelIdealMHDFactory();

  /// CFL number
  /**
   */
  virtual void CFL(Real a_cfl);

  /// Physical dimension of the longest side of the domain
  /**
   */
  virtual void domainLength(RealVect a_domainLength);

  /// ????
  /**
   */
  virtual void verbosity(const int& verbosity);


  /// Tag buffer size
  /**
   */
  void tagBufferSize(int a_tagBufferSize);

  /// Initial dt multiplier
  /**
   */
  void initialDtMultiplier(Real a_initialDtMultiplier);

  /// Set the gradient outputs flag
  /**
   */
  void output_density_gradient( bool a_output_density_gradient );
  void output_B_gradient( bool a_output_B_gradient );

  /// PatchMHDAM object (used as a factory)
  /**
   */
  virtual void patchMHDAM(RefCountedPtr<PatchMHDAM> a_patchMHDAM);

  /// Check that everything is defined
  /**
   */
  bool isDefined() const;
  
  /// Read parameters
  /**
   */
  void input( ParmParse & parser, int verbosity );


protected:
  // Some default values
  void setDefaultValues();

  int m_verbosity;

  // CFL number
  Real m_cfl;
  bool m_cflSet;

  // Physical dimension of the longest side of the domain
  RealVect m_domainLength;
  bool     m_domainLengthSet;

  // Tag buffer size
  int  m_tagBufferSize;
  bool m_tagBufferSizeSet;

  // Initial dt multiplier
  Real m_initialDtMultiplier;
  bool m_initialDtMultiplierSet;

  // Patch integrator (used as a factory)
  RefCountedPtr<PatchMHDAM> m_patchMHDAM;
  bool                      m_patchMHDAMSet;

  // true if we need to output gradients in tecplot format
  bool m_output_density_gradient;
  bool m_output_B_gradient;
  bool m_output_divB;
  
  // true if vectors are written in cartesian coordinates
  bool m_output_vecCS;

private:
  // Disallowed for all the usual reasons
  void operator=(const AMRLevelIdealMHDFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }

  // Disallowed for all the usual reasons
  AMRLevelIdealMHDFactory(const AMRLevelIdealMHDFactory& a_input)
  {
    MayDay::Error("invalid operator");
  }
};

#endif
