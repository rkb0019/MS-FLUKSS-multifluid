
#include "EqSysEuler.H"

#include "PatchEulerF_F.H"
#include "LGintegrator.H"


EqSysEuler::EqSysEuler(int a_ts)
: EquationSystem()
{    
  
  int nPrim   = WNUM_E + a_ts; 
  Interval lvlsStateInterval, lvlsPrimInterval;
  
  Interval consStateInterval(0,UNUM_E-1);
  
  if (a_ts > 0)
  {
    lvlsStateInterval.define(consStateInterval.end()+1, consStateInterval.end()+a_ts);
    lvlsPrimInterval = lvlsStateInterval;
  }
  
  define(consStateInterval,    // conservative variables interval,
         nPrim,                // number of primitive variables,     
         a_ts, 
         lvlsStateInterval,                
         lvlsPrimInterval);     
}

EqSysEuler::~EqSysEuler()
{  
  
}


// Number of conserved variables
int EqSysEuler::numStates()
{
  return UNUM_E;
}

// Names of the conserved variables
Vector<string> EqSysEuler::stateNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-momentum");
  retval.push_back("y-momentum");
  retval.push_back("z-momentum");
  retval.push_back("energy-density");
  return retval;
}

Vector<string> EqSysEuler::primitiveNames()
{
  Vector<string> retval;

  retval.push_back("density");
  retval.push_back("x-vel");
  retval.push_back("y-vel");
  retval.push_back("z-vel");
  retval.push_back("p");

  return retval;
}


// Number of flux variables
int EqSysEuler::numFluxes()
{
  // In some computations there may be more fluxes than conserved variables
  return UFLU_E;
}

// Number of primitive variables
int EqSysEuler::numPrimitives()
{
  // This doesn't equal the number of conserved variables because
  // auxiliary/redundant variable may be computed and stored
  return WNUM_E;
}

//  Number of primitive variables for which slopes are computed
int EqSysEuler::numSlopes()
{
  // This may be less than the number of primitive variables for the
  // reason given in numPrimitives() comments
  return WSLO_E;
}

// Compute the primitive variables from the conserved variables
void EqSysEuler::stateToPrim(       FArrayBox& a_W,
                             const FArrayBox& a_U,
                             const Box&       a_box)
{  
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  int startRho  = URHO;
  int startRhoW = WRHO;

  FORT_CONSTOPRIM_E( CHF_FRA(a_W),
                     CHF_CONST_FRA(a_U),
                     CHF_CONST_INT(startRho),
                     CHF_CONST_INT(startRhoW),
                     CHF_BOX(a_box));
}

// Compute the conserved variables from the primitive variables
void EqSysEuler::primToState(       FArrayBox& a_U,
                             const FArrayBox& a_W,
                             const Box&       a_box)
{  
  CH_assert(a_U.box().contains(a_box));
  CH_assert(a_W.box().contains(a_box));

  int startRho  = URHO;
  int startRhoW = WRHO;

  FORT_PRIMTOCONS_E( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_INT(startRho),
                     CHF_CONST_INT(startRhoW),
                     CHF_BOX(a_box));
}

  // Transform a_dWLeft and a_dWRight from primitive to characteristic variables
void EqSysEuler::charAnalysis(       FArrayBox & a_dWLeft,
                                     FArrayBox & a_dWRight,
                               const FArrayBox & a_W,
                               const int &       a_dir,
                               const Box &       a_box)
{
  FORT_CHARANALYSIS_E( CHF_FRA(a_dWLeft),
                       CHF_FRA(a_dWRight),
                       CHF_CONST_FRA(a_W),
                       CHF_CONST_INT(a_dir),
                       CHF_BOX(a_box));
}

  // Transform a_dWLeft and a_dWRight from characteristic to primitive variables
void EqSysEuler::charSynthesis(       FArrayBox & a_dWLeft,
                                      FArrayBox & a_dWRight,
                                const FArrayBox & a_W,
                                const int &       a_dir,
                                const Box &       a_box)
{
  FORT_CHARSYNTHESIS_E( CHF_FRA(a_dWLeft),
                        CHF_FRA(a_dWRight),
                        CHF_CONST_FRA(a_W),
                        CHF_CONST_INT(a_dir),
                        CHF_BOX(a_box));
}

void EqSysEuler::vectorVars(Vector<int> & a_vars) const
{
  a_vars.clear();
  a_vars.push_back(UMOMX);
}


void EqSysEuler::velocityVars(Vector<int> & a_vars) const
{
  a_vars.clear();
  a_vars.push_back(UMOMX);
}
