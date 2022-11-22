#include "AMRLevel.H"
#include "AMRLevelIdealMHD.H"
#include "AMRLevelIdealMHDFactory.H"

AMRLevelIdealMHDFactory::AMRLevelIdealMHDFactory()
{
  setDefaultValues();
}

// Virtual constructor
AMRLevel* AMRLevelIdealMHDFactory::new_amrlevel() const
{
  // Make sure everything is defined
  CH_assert(isDefined());

  // Create a new AMRLevelIdealMHD
  AMRLevelIdealMHD* amrMHDPtr = new AMRLevelIdealMHD();

  // Set up new object
  amrMHDPtr->CFL(m_cfl);
  amrMHDPtr->domainLength(m_domainLength);  
  amrMHDPtr->tagBufferSize(m_tagBufferSize);
  amrMHDPtr->initialDtMultiplier(m_initialDtMultiplier);
  amrMHDPtr->patchMHDAM(m_patchMHDAM);
  amrMHDPtr->verbosity(m_verbosity);  
  
  
  amrMHDPtr->output_density_gradient( m_output_density_gradient );
  amrMHDPtr->output_B_gradient( m_output_B_gradient );
  
  amrMHDPtr->m_output_vecCS = m_output_vecCS;
  amrMHDPtr->m_output_divB  = m_output_divB;

  // Return it
  return (static_cast <AMRLevel*> (amrMHDPtr));
}

AMRLevelIdealMHDFactory::~AMRLevelIdealMHDFactory()
{
}

void AMRLevelIdealMHDFactory::input( ParmParse & parser, int verbosity )
{
  int density_gradient = 0;
  parser.query ("output_density_gradient",density_gradient);

  int B_gradient = 0;
  parser.query ("output_B_gradient",B_gradient);
  
  output_density_gradient( density_gradient );
  output_B_gradient( B_gradient );  
  
  int output_vecCS = 0;
  parser.query ("output_vec_cs",output_vecCS);
  m_output_vecCS = (output_vecCS == 1);
  
  int output_divb = 0;
  parser.query ("output_divb",output_divb);
  m_output_divB = (output_divb == 1);


}

// CFL number
void AMRLevelIdealMHDFactory::CFL(Real a_cfl)
{
  m_cfl    = a_cfl;
  m_cflSet = true;
}

void AMRLevelIdealMHDFactory::verbosity(const int& a_verbosity)
{
  m_verbosity = a_verbosity;
}

// Physical dimension of the longest side of the domain
void AMRLevelIdealMHDFactory::domainLength(RealVect a_domainLength)
{
  m_domainLength    = a_domainLength;
  m_domainLengthSet = true;
}


// Tag buffer size
void AMRLevelIdealMHDFactory::tagBufferSize(int a_tagBufferSize)
{
  m_tagBufferSize    = a_tagBufferSize;
  m_tagBufferSizeSet = true;
}

// Initial dt multiplier
void AMRLevelIdealMHDFactory::initialDtMultiplier(Real a_initialDtMultiplier)
{
  m_initialDtMultiplier    = a_initialDtMultiplier;
  m_initialDtMultiplierSet = true;
}

// Set the density gradient output flag
void AMRLevelIdealMHDFactory::output_density_gradient( bool a_output_density_gradient )
{	m_output_density_gradient=a_output_density_gradient;}

// Set the B gradient output flag
void AMRLevelIdealMHDFactory::output_B_gradient( bool a_output_B_gradient )
{	m_output_B_gradient=a_output_B_gradient;}

// PatchMHDAM object (used as a factory)
void AMRLevelIdealMHDFactory::patchMHDAM(RefCountedPtr<PatchMHDAM> a_patchMHDAM)
{
  m_patchMHDAM    = a_patchMHDAM;
  m_patchMHDAMSet = true;
}

// Check that everything is defined
bool AMRLevelIdealMHDFactory::isDefined() const
{
  return (m_cflSet &&
          m_domainLengthSet &&          
          m_tagBufferSizeSet &&
          m_initialDtMultiplierSet &&
          m_patchMHDAMSet);
}

// Some default values
void AMRLevelIdealMHDFactory::setDefaultValues()
{
  CFL                (0.8);
  domainLength       (RealVect(D_DECL(1.0,1.0,1.0)));
  
  tagBufferSize      (2  );
  initialDtMultiplier(0.1);

  m_patchMHDAMSet        = false;
  m_verbosity            = 0;

  m_output_density_gradient = false;
  m_output_B_gradient       = false;
  
  m_output_vecCS            = false;
}
