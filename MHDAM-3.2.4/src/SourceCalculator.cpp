
#include "SourceCalculator.H"



//                                                                    Costructor
SourceCalculator::SourceCalculator()
{
  m_verbosity = 0;
}
//                                                                    Destructor
SourceCalculator::~SourceCalculator()
{
}
//                               Factory method - this object is its own factory
SourceCalculator * SourceCalculator :: new_SourceCalculator( void )
{
  SourceCalculator* retval = new SourceCalculator();

  return retval;
}
//                              Set time parameters for source terms calculation
//                                                              Input parameters
void SourceCalculator :: input( ParmParse & parser, int verbosity )
{
  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Source term calculator parameters:"      << endl;
  }
  
  m_verbosity = verbosity;
}

// Initialization of external source calculator.
void SourceCalculator::initExternalSC(Vector<AMRLevel*> a_levels)
{
}

// Initial values of source terms 
// By default, set the values to zero.
void SourceCalculator::initialize(LevelData<FArrayBox>& a_S)
{
  if (!a_S.isDefined()) return;

  for(DataIterator dit = a_S.dataIterator(); dit.ok(); ++dit)
  {
    a_S[dit()].setVal(0.0);
  }
}



// Check time for source terms calculation  
bool SourceCalculator::checkTime(Real dTime, int iStep)
{
  return false;
}

//                                                     Add external source terms
void SourceCalculator :: addExternalSources(       FArrayBox & a_U,
                                             const FArrayBox & a_S,
                                             const FArrayBox & a_W,
                                             const Real      & a_dt,
                                             const Real      & a_dx,                                             
                                             const Box       & a_box)
{
}

void SourceCalculator::CalculateSources(
                                Vector<AMRLevel*> a_levels,
                                Real a_time,
                                int  a_curStep)
{   
}

// Write information for restarting source calculator
void SourceCalculator::writeCheckpointFile(const char* chkFileName) const
{
}

// Number additional variables for writing to plot file
int SourceCalculator::numPlotVars()
{
  return 0;
}

//               Generate default names for the primitive variables, "variable#"
Vector<std::string> SourceCalculator::plotNames()
{
  Vector<std::string> retval;

  int cnum = numPlotVars();

  for (int ivar = 0; ivar < cnum; ivar++)
  {
    char varNameChar[80];
    sprintf(varNameChar,"variable%d",ivar);
    retval.push_back(std::string(varNameChar));
  }

  return retval; 
}

