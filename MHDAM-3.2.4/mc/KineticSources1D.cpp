#include "KineticSources1D.H"
#include "KineticSources1DF_F.H"
#include "AMRLevelIdealMHD.H"
#include "PatchIdealMHDF_F.H"
#include "LGintegrator.H"
#include "KS1DIntergrator.H"
#include "EosCommon.H"
#include "FORT_PROTO.H"


#include "DebugF_F.H"
#include "TecplotIO.H"

#include <list>
  
extern "C"
{

// Prototype for Fortran procedure mc_init 
#define FORT_MC_INIT_1D FORTRAN_NAME( RUN_MC_INIT_1D ,run_mc_init_1d )
void FORT_MC_INIT_1D
(
    int* const p_nz_in, 
    Real* const p_zmin_in, 
    Real* const p_zmax_in, 
    int* const restart_in, 
    char* const loadfile, 
    int* const total_neutrals,
    Real* const lism_nh_in, 
    Real* const lism_vh_in, 
    Real* const lism_th_in,
    int* const verbosity_in
);

// Prototype for Fortran procedure mc_neutrals
#define FORT_MC_NEUTRALS_1D FORTRAN_NAME( RUN_MC_NEUTRALS_1D, run_mc_neutrals_1d)
void FORT_MC_NEUTRALS_1D(
    Real* const run_time, 
    Real* const min_nchex, 
    Real* const p_dens_in, 
    Real* const p_ux_in, 
    Real* const p_uy_in, 
    Real* const p_uz_in, 
    Real* const p_temp_in,
    Real* const p_region_in,
    Real* const src_momx,
    Real* const src_momy,
    Real* const src_momz,
    Real* const src_energy,
    int* const  n_elements,
    Real* const grid_dens1, 
    Real* const grid_ux1, 
    Real* const grid_uy1, 
    Real* const grid_uz1, 
    Real* const grid_temp1
);

#define FORT_MC_OUTPUT_RAW_1D FORTRAN_NAME( RUN_MC_OUTPUT_RAW_1D ,run_mc_output_raw_1d )
void FORT_MC_OUTPUT_RAW_1D
(    
  const char* const chkfile     
);


}



                                                                   // Costructor
KineticSources1D::KineticSources1D()
{
  m_iTotalNeutrals = 0;
    
  m_dKineticRuntimeYears = -1.0;        
  m_dTimeFactor = 1.0;                 
  
  m_dKineticTimeIntervalYears = -1.0;   
  m_dKineticTimeIntervalMHDAM = -1.0;   
  m_dNextKineticTimeMHDAM = 0.0;   
  
  m_iKineticStepInterval = -1;   
  m_iKineticRuntimeConsistent = -1;       
}
                                                                   // Destructor
KineticSources1D::~KineticSources1D()
{
}
                              // Factory method - this object is its own factory
SourceCalculator * KineticSources1D :: new_SourceCalculator( void )
{
  SourceCalculator* retval = new KineticSources1D();

  return retval;
}

//                                             Number of calculated source terms
//                                                    and neutral parameters
int KineticSources1D :: numSourceTerms( void )
{
  return 7;
}
                                                             // Input parameters
void KineticSources1D :: input( ParmParse & parser, int verbosity )
{   
  m_dTemperatureL = 5000.0;  
  m_dNetnumL      = 1.0;   
  m_dVelxL        = 25000.0;
  
  parser.query( "velxL",        m_dVelxL );    
  parser.query( "temperatureL", m_dTemperatureL );
  parser.query( "netnumL",      m_dNetnumL      );    
    
  
  m_dKineticRuntimeYears = 1.0;
  parser.query("kinetic_runtime",  m_dKineticRuntimeYears);    
  
  m_iKineticRuntimeConsistent = -1;
  parser.query("kinetic_runtime_consistent",  m_iKineticRuntimeConsistent);    
  if (m_iKineticRuntimeConsistent != 1) m_iKineticRuntimeConsistent = -1;  
  
  
  m_dKineticTimeIntervalYears = -1.0;     
  parser.query("kinetic_time_interval",  m_dKineticTimeIntervalYears);    
  m_dTimeFactor = (eos_AU/m_dVelxL)/(60.0*60.0*24.0*365.0);   
  m_dKineticTimeIntervalMHDAM = m_dKineticTimeIntervalYears/m_dTimeFactor;     // Transform run time to our units.
  m_dNextKineticTimeMHDAM = 0.0;
  
  m_iKineticStepInterval = -1;
  parser.query("kinetic_step_interval",  m_iKineticStepInterval);        
    
  m_iTotalNeutrals = 20000;
  parser.query( "total_neutrals",  m_iTotalNeutrals);
  
  if (parser.contains("kinetic_restart_file"))  parser.query("kinetic_restart_file", m_sKineticRestartFile);    
  
  m_verbosity = verbosity;
  
  m_dExtraLength = 600.0;
  
  //m_verbosity = 0;
 
                                                             // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << endl << "Source term calculator parameters:"      << endl;        
    pout() << "    kinetic_runtime       = " << m_dKineticRuntimeYears  << endl;
    pout() << "    runtime consistent    : " << (m_iKineticRuntimeConsistent == 1 ? "yes" : "no")  << endl;
    pout() << "    total_neutrals        = " << m_iTotalNeutrals  << endl;
    pout() << "    kinetic_restart_file  = " << m_sKineticRestartFile  << endl;
    pout() << "    kinetic_time_interval = " << m_dKineticTimeIntervalYears  << endl;
    pout() << "    kinetic_step_interval = " << m_iKineticStepInterval  << endl;    
  }
  
}


// Initialization of external source calculator.
void KineticSources1D :: initExternalSC(Vector<AMRLevel*> a_levels)
{
  Vector<AMRLevel*>& amrlevels = a_levels;
  AMRLevelIdealMHD* Level0=static_cast<AMRLevelIdealMHD*>(amrlevels[0]);  
  const Box& Level0Box = amrlevels[0]->problemDomain().domainBox();
  
  Real AU2Meters = 0.01*eos_AU;
  
  Real dx = Level0->getpatchMHDAM()->getPhysProblem()->coordinateSystem()->dx(0,0);
  
  m_iNumExtraCells = (int) (m_dExtraLength/dx);
  pout() << "    m_iNumExtraCells = " << m_iNumExtraCells  << endl;    
    
  int p_nz = Level0Box.size(0)+m_iNumExtraCells;
  Real p_zmin = 0.0;                             // AU2Meters*(0.5*Level0->dx());
  Real p_zmax = AU2Meters*( p_nz*dx ); // AU2Meters*( (0.5+p_nz)*Level0->dx() );  
  int  total_neutrals = m_iTotalNeutrals;
  
  int  restart = (m_sKineticRestartFile.size() > 0 ? 1 : 0);
  char loadfile[60]={0}; loadfile[0] = 0;
  if (restart==1) m_sKineticRestartFile.copy(loadfile,m_sKineticRestartFile.size());
  
  
    
  FORT_MC_INIT_1D(&p_nz, &p_zmin, &p_zmax, &restart,   loadfile,   &total_neutrals, &m_dNetnumL, &m_dVelxL, &m_dTemperatureL,&m_verbosity);
}


//                                       Check time for source terms calculation
bool KineticSources1D :: checkTime( Real dTime, int iStep )
{
  //return true;
  
  bool bRetCode  = false;
  m_cur_step = iStep;
  
      
  if( m_dKineticTimeIntervalMHDAM > 0.0 )
  {
    if( dTime >= m_dNextKineticTimeMHDAM )
    {
      m_dNextKineticTimeMHDAM += m_dKineticTimeIntervalMHDAM;
      bRetCode     = true;
    }
    if( dTime >= m_dNextKineticTimeMHDAM )
    {
      m_dNextKineticTimeMHDAM = dTime + m_dKineticTimeIntervalMHDAM;
      bRetCode     = true;
    }
  }

  if( m_iKineticStepInterval > 0 )
  {
    if( iStep - (iStep/m_iKineticStepInterval)*m_iKineticStepInterval == 0 )
    {
       bRetCode    = true;
    }
 }

  return bRetCode;
}

// number of variables that are needed for external source calculator
int KineticSources1D::numVarsForSC()
{
  return 4;
}

// Using conservative variables a_U, the method prepares data (stored in a_Buf)
// that will be used in an external source calculator.
void KineticSources1D::PrepareDataForSC(const FArrayBox& a_U, FArrayBox& a_Buf) 
{
  const Box& b = a_Buf.box();
  FORT_PREPAREDATAFORKINETICSC1D(
      CHF_CONST_FRA(a_U),
      CHF_FRA(a_Buf),
      CHF_BOX(b));
}

// Calls external source calculator. Source terms are returned in a_SourceTerms
void KineticSources1D::CallExternalSC(const FArrayBox& Level0Data, Real a_dt, Real a_dx, FArrayBox& a_SourceTerms) 
{
  const Box& Level0Box = Level0Data.box();
  Real min_nchex;  
      
  int arraysize = Level0Box.size(0)+m_iNumExtraCells;
  
  Real* p_dens_in = new Real [arraysize];
  Real* p_ux_in   = new Real [arraysize];
  Real* p_uy_in   = new Real [arraysize];
  Real* p_uz_in   = new Real [arraysize];
  Real* p_temp_in = new Real [arraysize];
  
  Real* p_region_in = new Real [arraysize];  
  
  Real* src_momx   = new Real [arraysize];
  Real* src_momy   = new Real [arraysize];
  Real* src_momz   = new Real [arraysize];
  Real* src_energy = new Real [arraysize];  
  
  Real* grid_dens1 = new Real [arraysize];
  Real* grid_ux1   = new Real [arraysize];
  Real* grid_uy1   = new Real [arraysize];
  Real* grid_uz1   = new Real [arraysize];
  Real* grid_temp1 = new Real [arraysize];
  
  IntVect aux_IntVect;int i,j;
  IntVect LCorner=Level0Box.smallEnd();
  IntVect UCorner=Level0Box.bigEnd();
  
  int curProc = procID();  
    
  if (curProc==0)
  {    
            
    for (i=LCorner[0];i<=UCorner[0];i++)
    { 
      aux_IntVect[0]=i;aux_IntVect[1]=1;		            
            
      p_dens_in[i]   = Level0Data.get(aux_IntVect,SWRHO);
      p_ux_in[i]     = Level0Data.get(aux_IntVect,SWVELY);
      p_uy_in[i]     = 0.0;
      p_uz_in[i]     = Level0Data.get(aux_IntVect,SWVELX);
      p_temp_in[i]   = Level0Data.get(aux_IntVect,SWTEMP);
      p_region_in[i] = 1.1;          
    }

    aux_IntVect[0]=UCorner[0];aux_IntVect[1]=1;
    for (i=UCorner[0];i<arraysize;i++)
    { 
      p_dens_in[i]   = Level0Data.get(aux_IntVect,SWRHO);
      p_ux_in[i]     = Level0Data.get(aux_IntVect,SWVELY);
      p_uy_in[i]     = 0.0;
      p_uz_in[i]     = Level0Data.get(aux_IntVect,SWVELX);
      p_temp_in[i]   = Level0Data.get(aux_IntVect,SWTEMP);
      p_region_in[i] = 1.1;          
    }      
  }
  
  Real TimeStep = m_dKineticRuntimeYears;
  if (m_iKineticRuntimeConsistent == 1) 
  {
    //pout() << "dt: " << a_dt << ", m_dTimeFactor: " << m_dTimeFactor;
    TimeStep = a_dt*m_dTimeFactor;
  }
  
  
  int n_elements = arraysize;
  
  if ( m_verbosity >= 2 )
  {
    pout() <<  "Call MC_NEUTRALS (runtime = "<< TimeStep << " years)"  << endl;        
  }
      
  FORT_MC_NEUTRALS_1D(&TimeStep, &min_nchex,
          p_dens_in, p_ux_in, p_uy_in, p_uz_in, p_temp_in,
          p_region_in,  
          src_momx, src_momy, src_momz, src_energy,
          &n_elements,
          grid_dens1, grid_ux1, grid_uy1, grid_uz1, grid_temp1);
  
  if ( m_verbosity >= 2 )
  {
    pout() <<  "End MC_NEUTRALS"  << endl << endl;        
  }    
  
  if ((curProc==0)&&(m_verbosity >= 4))
  {
    pout() <<  "Source terms: (momx, momy, momz, energy)"  << endl;        
    char buf[50];
    for (i=LCorner[0];i<=UCorner[0];i++) 
    {
      pout() << i << " ";
      sprintf(buf,"%.6e",src_momx[i]); pout() << buf << " " ;
      sprintf(buf,"%.6e",src_momy[i]); pout() << buf << " " ;
      sprintf(buf,"%.6e",src_momz[i]); pout() << buf << " " ;
      sprintf(buf,"%.6e",src_energy[i]);   pout() << buf << endl;        
    }
  }

  
  if (curProc==0)
  {                
    for (i=LCorner[0];i<=UCorner[0];i++)
    for (j=LCorner[1];j<=UCorner[1];j++)
    { 
      aux_IntVect[0]=i;aux_IntVect[1]=j;		            
      
      a_SourceTerms.set(aux_IntVect,SMOMX, src_momz[i-LCorner[0]]);
      a_SourceTerms.set(aux_IntVect,SMOMY, src_momx[i-LCorner[0]]);
      a_SourceTerms.set(aux_IntVect,SENG,  src_energy[i-LCorner[0]]);      
      a_SourceTerms.set(aux_IntVect,NRHO,  grid_dens1[i-LCorner[0]]);      
      a_SourceTerms.set(aux_IntVect,NVELX, grid_uz1[i-LCorner[0]]);      
      a_SourceTerms.set(aux_IntVect,NVELY, grid_ux1[i-LCorner[0]]);      
      a_SourceTerms.set(aux_IntVect,NTEMP, grid_temp1[i-LCorner[0]]);      
    }      
    
  //#ifndef NDEBUG
    if ( m_verbosity >= 1 )
    {
      char file_name[50];  
          
      FArrayBox tmpST(a_SourceTerms.box(),a_SourceTerms.nComp()); 
      tmpST.copy(a_SourceTerms);
      ModifySourceTerms(tmpST, a_SourceTerms.box());            
      
      
    
      if (m_cur_step%100 == 0)
      {                
        sprintf(file_name,"mhdk_st%06i.dat",m_cur_step);
        FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"source terms\"\n VARIABLES=\"X\" \"Y\" \"SMOMX\" \"SMOMY\" \"SENG\" \"NRHO\" \"NVELX\" \"NVELY\" \"NTEMP\" ");
        WriteFArrayBoxToTecplotFile(tfile, tmpST, a_SourceTerms.box(), Interval(0,numSourceTerms()-1), a_dx);
        CloseTecplotFile(tfile);
      
        //sprintf(file_name,"mhdk_st_orig%04i.dat",m_cur_step);
        //FILE* tfile1 = OpenTecplotFile(file_name,"TITLE = \"source terms\"\n VARIABLES=\"X\" \"Y\" \"SMOM\" \"SENG\" \"NRHO\" \"NVEL\" \"NTEMP\" ");        
        //WriteFArrayBoxToTecplotFile(tfile1, a_SourceTerms, a_SourceTerms.box(), Interval(0,4), a_dx);    
        //CloseTecplotFile(tfile1); 
      }
     }
     //#endif
  }
  delete[] p_dens_in;
  delete[] p_ux_in;
  delete[] p_uy_in;
  delete[] p_uz_in;
  delete[] p_temp_in;
  delete[] p_region_in;  
  delete[] src_momx;
  delete[] src_momy;
  delete[] src_momz;
  delete[] src_energy;  
  delete[] grid_dens1;
  delete[] grid_ux1;
  delete[] grid_uy1;
  delete[] grid_uz1;
  delete[] grid_temp1;
}

// Modification of source terms calculated by external source calculator.
void KineticSources1D::ModifySourceTerms(FArrayBox& a_SourceTerms, const Box& a_box)
{  
  FORT_MODIFYSOURCEAFTERKINETICSC1D(
        CHF_FRA(a_SourceTerms),
        CHF_BOX(a_box));            

}


//                                                     Add external source terms
void KineticSources1D :: addExternalSources(       FArrayBox & a_U,
                                           const FArrayBox & a_S,
                                           const FArrayBox & a_W,
                                           const Real      & a_dt,
                                           const Real      & a_dx,
                                           const Box       & a_box)
{
  FORT_ADDKINETICSOURCE1D( CHF_FRA(a_U),
                         CHF_CONST_FRA(a_S),
                         CHF_CONST_REAL(a_dt),
                         CHF_BOX(a_box) );
}

void KineticSources1D::writeCheckpointFile(const char* chkFileName) const
{
  //return;
  FORT_MC_OUTPUT_RAW_1D(chkFileName);
}
