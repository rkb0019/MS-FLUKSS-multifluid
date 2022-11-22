#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <limits>

using std::ifstream;
using std::ios;

#include "DebugOut.H"
#include "DebugF_F.H"

#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "SWLISMProblem.H"
#include "SWLISMF_F.H"
#include "ChargeExchange2F_F.H"
#include "bedfordcxF_F.H"
#include "bedford4FPIREGF_F.H"

#include "PatchMHDAMF_F.H"
#include "PatchIdealMHDF_F.H"
#include "RiemannSolver.H"
#include "LGintegrator.H"
#include "MHDAMDefs.H"
#include "TecplotIO.H"
#include "EqSysMHDMF.H"
#include "EosCommon.H"
#include "DednerF_F.H"
#include "SWLISMTurbF_F.H"

// Null constructor
SWLISMProblem::SWLISMProblem()
  : MultiFluidProblem()
{
  m_physModel          = PP_Undefined;  
  m_isFortranCommonSet = false;
    
  m_stop_ref           = -1.0;
  m_max_levelTS        = 10000;  
  m_region_tracer      = false;
  m_subproblem         = 0;  
  m_dNextTurbLISMTime  = 0.0;
  
  m_photoionize = false;  
  m_const_H     = false;

}

SWLISMProblem::SWLISMProblem(ePhysicalModel ePM)
  : MultiFluidProblem()
{
  m_physModel          = ePM;  
  m_isFortranCommonSet = false;
  m_stop_ref           = -1.0;
  m_max_levelTS        = 10000;  
  m_region_tracer      = false;
  m_subproblem         = 0;  
  m_dNextTurbLISMTime  = 0.0;
  
  m_photoionize = false;  
  m_const_H     = false;

}


// Input parameters
void SWLISMProblem::input( ParmParse & parser, int verbosity )
{  
  // Read problem specific parameters here
  
  parser.query( "gamma",   m_gamma );
  parser.query( "lismN",   m_lismN );
  parser.query( "lismV",   m_lismV );
  parser.query( "lismT",   m_lismT );
  parser.query( "lismB",   m_lismB );
  parser.query( "sunN",    m_sunN  );
  parser.query( "sunV",    m_sunV  );
  parser.query( "sunT",    m_sunT  );
  parser.query( "initR",   m_initR );
  parser.query( "subproblem",     m_subproblem     );
  
  std::string problemString;
  parser.query("problem",problemString);              
  if (problemString == "swlismKinetic") m_subproblem = SW_KINETIC;

  
  Real rho        = m_lismN*eos_mp;
  Real pref       = rho*m_lismV*m_lismV;
  Real p          = 2.0*eos_k*m_lismN*m_lismT;
  Real lismP      = p/pref;
  m_lismM         = 1.0/sqrt( m_gamma*lismP );
  
  if (m_physModel == PP_Undefined)
  {
    int iFluids  = 1;
    
    parser.query( "fluids",  iFluids   );

    switch( iFluids ){
    case 2  : m_physModel  = PP_2FluidPM; break;
    case 3  : m_physModel  = PP_3FluidPM; break;
    case 4  : m_physModel  = PP_4FluidPM; break;
    default : m_physModel  = PP_MHDPM;}
  }
  
  if (m_subproblem == SW_TURB)
  {
    Real timeStep;  
    m_dTimeFactor = (eos_AU/m_lismV)/(60.0*60.0*24.0*365.0);   
        
    parser.query( "TurbLISMTimeStep", timeStep);  
    m_dTurbLISMTimeStep = timeStep/m_dTimeFactor;
    
    parser.query( "TurbSWTimeStep",   timeStep );  
    m_dTurbSWTimeStep   = timeStep/m_dTimeFactor;
    
    parser.query( "LISMDeviation", m_dLISMDeviation );  
    parser.query( "SWDeviation",   m_dSWDeviation );  
    
    parser.query( "R0Turb",   m_R0Turb );    
    
    m_dNextTurbLISMTime = 0.0;
  }
  
  m_TMLIM = 50000.0; // Temperature separating regions 1 and 2
  parser.query( "TMLIM", m_TMLIM);
    
  int nRegionTracer = -1;
  parser.query("region_tracer",nRegionTracer);  
  m_region_tracer = (nRegionTracer > 0);
  
  int photoionize = 0;
  parser.query( "photoionize",   photoionize   );
  m_photoionize = (photoionize==1);
  
  int const_H = 0;
  parser.query( "const_H",   const_H);
  m_const_H = (const_H==1);
  
  int output_region = -1;
  parser.query( "output_region", output_region);
  m_output_region = output_region > 0;
      
      
  D_TERM(parser.query( "XC",      m_sunXYZ[0]    );,
         parser.query( "YC",      m_sunXYZ[1]    );,
         parser.query( "ZC",      m_sunXYZ[2]    ););
         
  parser.query( "R0",      m_R0    );
    
  parser.query( "stop_ref",  m_stop_ref);
  parser.query( "max_levelTS", m_max_levelTS);
  
  int i,j;
  
  // We need to sligtly modify m_sunXC, m_sunYC, m_sunZC in order 
  // to Sun position lies strictly at grid node
  
  Real domainLength = 1.0;
  parser.query("domain_length",domainLength);
  
  // Set the resolution of the coarsest level
  vector<int> numCells(SpaceDim); 
  for (i = 0; i < SpaceDim; ++i) numCells[i]=0;
  parser.queryarr("num_cells",numCells,0,SpaceDim);
      
  Real dx = domainLength/numCells[0];    
  for (i = 0; i < SpaceDim; ++i) m_sunXYZ[i] = ((int)floor(m_sunXYZ[i]/dx+0.5))*dx;
  
  m_R0dtIgnore = m_R0 - dx;

  Real safeDZ = 0.0;
  parser.query("RegSafeDZ",safeDZ);
  parser.query("RegionSafetyDZ",safeDZ);

  m_RegSafeZ   = dx*numCells[2] - m_sunXYZ[2] - safeDZ;

  m_netN = 0.0;
  if (nFluids()>1)  
  {
    if (!parser.contains("netN")) MayDay::Error("netN must be defined for multifluid calculations");    
  }
  parser.query( "netN",    m_netN  );
   
  
  m_verbosity = verbosity;

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "The SW interaction with the LISM problem input:" << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "lismN     = " << m_lismN     << endl;
    pout() << "lismV     = " << m_lismV     << endl;
    pout() << "sunN      = " << m_sunN      << endl;
    pout() << "sunV      = " << m_sunV      << endl;
    pout() << "sunT      = " << m_sunT      << endl;
    pout() << "initR     = " << m_initR     << endl;
    D_TERM(
    pout() << "XC        = " << m_sunXYZ[0]        << endl;,
    pout() << "YC        = " << m_sunXYZ[1]        << endl;,   
    pout() << "ZC        = " << m_sunXYZ[2]        << endl;); 
    pout() << "R0        = " << m_R0        << endl;
    pout() << "RegSafeZ  = " << m_RegSafeZ  << endl;
    if (nFluids()>1)
      pout() << "netN      = " << m_netN  << endl;
  }

  setFortranCommon();
  
      
}


// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void SWLISMProblem::setFortranCommon( )
{
  CH_assert(m_isFortranCommonSet == false);
  
  Real sunZC = 0.0;
  FORT_SETSWLISM( CHF_CONST_REAL( m_gamma ),
                  CHF_CONST_REAL( m_lismN ),
                  CHF_CONST_REAL( m_lismV ),
                  CHF_CONST_REAL( m_lismT ),
                  CHF_CONST_REAL( m_lismB ),
                  CHF_CONST_REAL( m_sunXYZ[0] ),
                  CHF_CONST_REAL( m_sunXYZ[1] ),
                  CHF_CONST_REAL( sunZC ),
                  CHF_CONST_REAL( m_R0    ),
                  CHF_CONST_REAL( m_sunN  ),
                  CHF_CONST_REAL( m_sunV  ),
                  CHF_CONST_REAL( m_sunT  ),
                  CHF_CONST_REAL( m_initR ),
                  CHF_CONST_REAL( m_netN  ),
                  CHF_CONST_REAL( m_TMLIM ),
                  CHF_CONST_REAL( m_RegSafeZ));

  if (nFluids()>=1)
  {                    
    Real scaleLen  = 1.5e+13;

    FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( m_lismV  ),
                         CHF_CONST_REAL( scaleLen ),
                         CHF_CONST_REAL( m_lismN  ) );  
  }
  
  if (m_subproblem == SW_TURB)
    FORT_SETSWLISMTURB( CHF_CONST_REAL( m_dLISMDeviation ),
                  CHF_CONST_REAL( m_dSWDeviation ),
                  CHF_CONST_REAL( m_R0Turb) );    
                         
  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void SWLISMProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this object.
 */
void SWLISMProblem::define(const ProblemDomain& a_domain,                           
                           const int            a_level)
{
  MultiFluidProblem::define(a_domain, a_level);
  
  int i;
  IntVect gridDim, LCorner, UCorner;
  
  Real dx = m_csh->dx(0,m_level);
  
  D_TERM(
    m_sunIJK[0] = (int)floor(m_sunXYZ[0]/dx+0.5);,
    m_sunIJK[1] = (int)floor(m_sunXYZ[1]/dx+0.5);,
    m_sunIJK[2] = (int)floor(m_sunXYZ[2]/dx+0.5););  
    
  CH_assert(fabs(m_sunIJK[0]*dx-m_sunXYZ[0])<1e-6);
    
  IntVect SunPos(m_sunIJK);    
  
  // Defining m_R0Box
  Real R0  = m_R0;
  int  iR0 = 1 + (int)floor(R0/dx+0.5);
  gridDim = IntVect( D_DECL( iR0,iR0,iR0 ));
    
  LCorner = SunPos; UCorner = SunPos;
  LCorner[0]-=gridDim[0];
  UCorner[0]+=gridDim[0]-1;        
  UCorner[1]+=gridDim[1]-1;        
  
  m_R0Box.define(LCorner,UCorner);    
  
  // Define data for path-through BC
  int fluids = nFluids();
  if (fluids > 1)
  {
    Box strip_a,strip_b;// '_a' means after the Sun, '_b' means before the Sun
    Vector<Box> vstrips_a, vstrips_b;
    for (i=0; i<iR0; i++)
    {
      strip_a.define(IntVect( D_DECL(SunPos[0]+i,    SunPos[1]    , SunPos[2])),
                     IntVect( D_DECL(SunPos[0]+i,    SunPos[1]+iR0, SunPos[2])));
                     
      strip_b.define(IntVect( D_DECL(SunPos[0]-i-1,  SunPos[1]    , SunPos[2])),
                     IntVect( D_DECL(SunPos[0]-i-1,  SunPos[1]+iR0, SunPos[2])));
                     
      vstrips_a.push_back(strip_a);
      vstrips_b.push_back(strip_b);
    }
    Vector<int> procs;
    LoadBalance(procs, vstrips_a);
    
    DisjointBoxLayout dbl_strips_a(vstrips_a, procs);
    DisjointBoxLayout dbl_strips_b(vstrips_b, procs);        
    
    int numVars = (fluids == 5 ? 2*UNUM_E : UNUM_E);
    
    m_lstrips.define(dbl_strips_b, numVars);
    m_ladjstr.define(dbl_strips_a, numVars);        
  }
  
  ParmParse parser("mhdam");            
  char buf[20];std::vector<Real> tmpVect;
  
  sprintf(buf, "kinetic_grid%ibb", a_level+1);
  if (parser.contains(buf))
  {
    parser.queryarr( buf, tmpVect, 0, CH_SPACEDIM);    
    parser.queryarr( buf, tmpVect, 0, CH_SPACEDIM);    
    RealVect rv(D_DECL(tmpVect[0], tmpVect[1], tmpVect[2]));
                    
    IntVect iLo(D_DECL((int)((m_sunXYZ[0]-rv[0])/dx), 0, 0));
    IntVect iHi(D_DECL((int)((m_sunXYZ[0]+rv[0])/dx), (int)(rv[1]/dx), 0));
    
    m_gridBox.define(iLo,iHi);   
    m_gridBox&=m_domain;     
        
  } else
  {
    sprintf(buf, "level_grid%ibb", a_level+1);
    if (parser.contains(buf))
    {
      parser.queryarr( buf, tmpVect, 0, 2*CH_SPACEDIM);    
      RealVect rLo(D_DECL(tmpVect[0], tmpVect[1], tmpVect[2]));
      RealVect rHi(D_DECL(tmpVect[CH_SPACEDIM], tmpVect[CH_SPACEDIM+1], tmpVect[CH_SPACEDIM+2]));
              
      rLo += m_sunXYZ;    
      rHi += m_sunXYZ;    
      
      IntVect iLo(D_DECL((int)(rLo[0]/dx), (int)(rLo[1]/dx), (int)(rLo[2]/dx)));
      IntVect iHi(D_DECL((int)(rHi[0]/dx), (int)(rHi[1]/dx), (int)(rHi[2]/dx)));
      
      m_gridBox.define(iLo,iHi);   
      m_gridBox&=m_domain;          
    }
  }
    
  
}


// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* SWLISMProblem::new_PhysProblem()
{
  SWLISMProblem* retval = new SWLISMProblem(m_physModel);
  
  retval->copy_PhysProblem(this);  
  
  return static_cast<PhysProblem*>(retval);
}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void SWLISMProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  const SWLISMProblem* PP = dynamic_cast<const SWLISMProblem*>(a_PP);
  if (PP == NULL) MayDay::Error("SWLISMProblem::copy_PhysProblem. Wrong argument");
  
  MultiFluidProblem::copy_PhysProblem(a_PP);
  
  if (PP->m_isFortranCommonSet == true)
  {
    this->m_sunXYZ     = PP->m_sunXYZ; 
    this->m_R0         = PP->m_R0;
    this->m_R0dtIgnore = PP->m_R0dtIgnore;

    this->m_lismN      = PP->m_lismN;
    this->m_lismM      = PP->m_lismM;
    this->m_lismV      = PP->m_lismV;
    this->m_lismT      = PP->m_lismT;
    this->m_lismB      = PP->m_lismB;
    this->m_sunN       = PP->m_sunN;
    this->m_sunV       = PP->m_sunV;
    this->m_sunT       = PP->m_sunT;
    this->m_initR      = PP->m_initR;
    this->m_stop_ref   = PP->m_stop_ref;
    this->m_max_levelTS = PP->m_max_levelTS;

    this->m_gamma      = PP->m_gamma;
    this->m_netN       = PP->m_netN;
    this->m_verbosity  = PP->m_verbosity;
    this->m_TMLIM      = PP->m_TMLIM;
    
    this->m_region_tracer = PP->m_region_tracer;
    
    this->m_subproblem = PP->m_subproblem;
    
    this->m_photoionize = PP->m_photoionize;
    this->m_const_H = PP->m_const_H;
    this->m_output_region = PP->m_output_region;
    
    if (m_subproblem == SW_TURB)
    {
      this->m_dTimeFactor       = PP->m_dTimeFactor;
      this->m_dTurbSWTimeStep   = PP->m_dTurbSWTimeStep;
      this->m_dTurbLISMTimeStep = PP->m_dTurbLISMTimeStep;
      this->m_dLISMDeviation    = PP->m_dLISMDeviation;
      this->m_dSWDeviation      = PP->m_dSWDeviation;  
      this->m_R0Turb            = PP->m_R0Turb;
    }
  
    this->setFortranCommonSet();
  }
            
}


                                                             // Fill ghost cells
void SWLISMProblem::fillGhostCells(          FArrayBox&      a_W,
                                       const FArrayBox&      a_U,
                                       const int&            a_dir,
                                       const Real&           a_time)
{
  fillGhostCellsDefault( a_W, a_U, a_dir, a_time);
  if (m_subproblem == SW_TURB)
  {
    fillGhostCellsZankTurb2007( a_W, a_U, a_dir, a_time);
  }
}                                       


                                                             // Fill ghost cells
void SWLISMProblem::fillGhostCellsDefault(          FArrayBox&      a_W,
                                       const FArrayBox&      a_U,
                                       const int&            a_dir,
                                       const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);
  
  int indW,indD,i;
  
  Box WBox = a_W.box();

                       // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  
  
  int jtop = m_domain.size(1)-1;
  
  int nGS;
  int nGSMax = 4; // How many two ghost cells must be filled
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  
  int iCP    = eqSys->correctionPotentialIndex();
  int iRegTr = -1;
  if  (m_region_tracer)
    iRegTr =  (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);  
  
  
  if (a_dir == 0)
  {
    // Right boundary
    indW = WBox.bigEnd( a_dir );
    indD =  tmp.bigEnd( a_dir );    
    
    nGS  = MIN(abs(indW-indD),nGSMax);
        
    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
    
                                                         // Fill the ghost cells      
      FORT_LISMINIT(
          CHF_FRA(a_W), 
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(iRegTr),
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_REAL(dx),
          CHF_BOX(boundaryBox));
      
      if (iCP > 0)
      FORT_SWLISMDEDNERGS(
          CHF_FRA(a_W), 
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_INT(iCP),            
          CHF_CONST_INT(jtop),            
          CHF_BOX(boundaryBox));
    }

    // Left boundary  
    indW = WBox.smallEnd( a_dir );
    indD =  tmp.smallEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
  
                                                         // Fill the ghost cells      
      FORT_LISMINIT(
          CHF_FRA(a_W), 
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(iRegTr),
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_REAL(dx),
          CHF_BOX(boundaryBox));  

      if (iCP > 0)
      FORT_SWLISMDEDNERGS(
          CHF_FRA(a_W), 
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_INT(iCP),     
          CHF_CONST_INT(jtop),                   
          CHF_BOX(boundaryBox));     
                
    }    
        
  }
  
  if (a_dir == 1)
  {
    indW = WBox.bigEnd( a_dir );
    indD =  tmp.bigEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    // Top boundary
    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
    
                                                         // Fill the ghost cells      
      FORT_LISMINIT(
          CHF_FRA(a_W), 
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(iRegTr),
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_REAL(dx),
          CHF_BOX(boundaryBox));
      
      if (iCP > 0)
      FORT_SWLISMDEDNERGS(
          CHF_FRA(a_W), 
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_INT(iCP),      
          CHF_CONST_INT(jtop),              
          CHF_BOX(boundaryBox));
      
    }
  
    indW = WBox.smallEnd( a_dir );
    indD =  tmp.smallEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
    
    // Axisymmetric case, bottom boundary
    if( indW != indD )
    {
      int sign         = -1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD-nGS, nGS );
  
      Box srcBox = WBox;
      Box desBox = WBox;
      int i_srcBox = indD;
      int i_desBox = indD - 1;
      
      for (i = 0; i < boundaryBox.size(a_dir); i++) 
      {
        srcBox.setRange( a_dir, i_srcBox, 1 );
        desBox.setRange( a_dir, i_desBox, 1 );
        a_W.copy(a_W, srcBox, 0, desBox, 0, a_W.nComp());        
        
        i_srcBox++;
        i_desBox--;
      }
                                                         
      a_W.mult(-1.0, boundaryBox, WVELY);
      a_W.mult(-1.0, boundaryBox, WBY);       
      
      for (int iFluid = 1; iFluid < fluids; ++iFluid)
      { 
        int iRho = eqSys->densityIndexPrim(iFluid);
        a_W.mult(-1.0, boundaryBox, iRho+WVELY);        
      }           
              
    }
  }
  #ifndef NDEBUG
    FORT_VIEWBOXDATA(
      CHF_FRA(a_W)
      );
  #endif
}

                                                            // Fill ghost cells
void SWLISMProblem::fillGhostCellsZankTurb2007(       FArrayBox&      a_W,
                                       const FArrayBox&      a_U,
                                       const int&            a_dir,
                                       const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  if ((a_dir == 1) || (m_dTurbLISMTimeStep < 0.0))
  {    
    return;
  }  
      
  int indW,indD,i;
  
  Real dx = m_csh->dx(0,m_level);
  
  Box WBox = a_W.box();

                       // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  
  
  int nGS;
  int nGSMax = 4; // How many two ghost cells must be filled
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  
  int iRegTr = -1;
  
  
  if (a_dir == 0)
  {              
    // Left boundary
    indW = WBox.smallEnd( a_dir );
    indD =  tmp.smallEnd( a_dir );
    
    nGS  = MIN(abs(indW-indD),nGSMax);
        
    if( indW != indD )        
    {   
      CH_assert(m_LISMTurb.nComp() == a_W.nComp());      
      
      int sign         =-1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
      
      // Calculate LISM params in one strip of ghost cells
      Box bStrip(
        IntVect(D_DECL(-1, tmp.loVect()[1],tmp.loVect()[2])),
        IntVect(D_DECL(-1, tmp.hiVect()[1],tmp.loVect()[2])));
              
      if( a_time >= m_dNextTurbLISMTime )
      {
        FORT_LISMBC_TURB(
            CHF_FRA(m_LISMTurb),        
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(sign),
            CHF_CONST_REAL(dx),  
            CHF_CONST_INT(iRhoN),     
            CHF_CONST_INT(fluids),
            CHF_BOX(bStrip));
            
        pout() << "Update LISM turbulence" << endl;    
      }
      
      // Duplicate the strip into all ghost cells        
      Box desBox = tmp;        
      int i_desBox = indD - 1;        
      for (i = 0; i < boundaryBox.size(a_dir); i++) 
      {          
        desBox.setRange( a_dir, i_desBox, 1 );
        a_W.copy(m_LISMTurb, bStrip, 0, desBox, 0, a_W.nComp());                            
        i_desBox--;
      }            
    }        
  }
  
  
  #ifndef NDEBUG
    FORT_VIEWBOXDATA(
      CHF_FRA(a_W)
      );
  #endif
}


// Set boundary fluxes
void SWLISMProblem::fluxBC(    FArrayBox&      a_F,
                                  FArrayBox&      a_Bn,
                            const FArrayBox&      a_WMinus,
                            const FArrayBox&      a_WPlus,
                            const int&            a_dir,
                            const Side::LoHiSide& a_side,
                            const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int fluids = nFluids();  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);  
  int iCP    = eqSys->correctionPotentialIndex();
  
  // In periodic case, this doesn't do anything
  if (!m_domain.isPeriodic(a_dir))
  {
    int sign;
    Box FBox = a_F.box();
    Box tmp = FBox;

    // Determine which side and thus shifting directions
    if (a_side == Side::Lo)
    {
      sign = -1;
    }
    else
    {
      sign = 1;
    }

    tmp.shiftHalf(a_dir,sign);

    // Is there a domain boundary next to this grid
    if (!m_domain.contains(tmp))
    {
      tmp &= m_domain;

      Box boundaryBox;

      // Find the strip of cells next to the domain boundary
      if (a_side == Side::Lo)
      {
        boundaryBox = bdryLo(tmp,a_dir);
      }
      else
      {
        boundaryBox = bdryHi(tmp,a_dir);
      }
      
      #ifndef NDEBUG
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WPlus)
        );
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WMinus)
        );          
      #endif
      
      // Pogorelov's boundary conditions      
      
      // boundaryBox is node centered. 
      // We need cell centered box that lies outside the problem domain.      
                      
      CH_assert(m_RS!=NULL);        
      if (fluids>1) CH_assert(m_RSGD!=NULL);
                  
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
      
      for (int iFluid = 1; iFluid < fluids; ++iFluid)
      { 
        int iRho = eqSys->densityIndexCons(iFluid);        
        m_RSGD->fluxes( a_F, a_WPlus, a_WMinus, a_dir, iRho, boundaryBox );
      }              
                   
                            
      FArrayBox& shiftWLeft  = (FArrayBox&)(a_WPlus);
      FArrayBox& shiftWRight = (FArrayBox&)(a_WMinus);
                          // Shift the left and right primitive variable boxes
                          // 1/2 cell so they are face centered
      shiftWLeft .shiftHalf(a_dir, 1);
      shiftWRight.shiftHalf(a_dir,-1);
      CH_assert( shiftWLeft.box().contains(boundaryBox) );
      CH_assert( shiftWRight.box().contains(boundaryBox) );        
      
      if (iCP == 0)
      FORT_SWLISMBC( CHF_FRA(a_F),
            CHF_FRA1(a_Bn,0),              
            CHF_CONST_FRA(shiftWLeft), 
            CHF_CONST_FRA(shiftWRight),
            CHF_CONST_INT(sign), 
            CHF_CONST_INT(a_dir), 
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),            
            CHF_BOX(boundaryBox) );
      else             
      {
        int iRho = URHO;
        //FORT_DEDNERBC(
        // CHF_FRA(a_F),
        // CHF_CONST_FRA(shiftWLeft), 
        // CHF_CONST_FRA(shiftWRight),
        // CHF_CONST_INT(sign), 
        // CHF_CONST_INT(a_dir),          
        // CHF_BOX(boundaryBox) );   

        //FORT_DEDNER2X2F(CHF_FRA(a_F),
        //      CHF_CONST_FRA(shiftWLeft),
        //      CHF_CONST_FRA(shiftWRight),
        //      CHF_CONST_INT(a_dir),
        //      CHF_CONST_INT(iRho),
        //      CHF_BOX(boundaryBox));

        FORT_SWLISMBCDEDNER( CHF_FRA(a_F),
            CHF_FRA1(a_Bn,0),              
            CHF_CONST_FRA(shiftWLeft), 
            CHF_CONST_FRA(shiftWRight),
            CHF_CONST_INT(sign), 
            CHF_CONST_INT(a_dir), 
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),            
            CHF_BOX(boundaryBox) );  
      }
                          // Shift the left and right primitive variable boxes
                          // back to their original position (no net change is made!)
      shiftWLeft .shiftHalf(a_dir,-1);
      shiftWRight.shiftHalf(a_dir, 1);
      
      // Axisymmetric case, bottom boundary
      if ((a_dir==1) && (a_side==Side::Lo)) a_Bn.setVal(0.0,boundaryBox,0); 
                                                             
    }
  }
}

// Set up initial conditions
void SWLISMProblem::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = m_csh->dx(0,m_level);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);  
  int iCP = eqSys->correctionPotentialIndex();  
  int iRegTr = -1;
  if  (m_region_tracer)
      iRegTr =  (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);  

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition on this grid    
    FORT_SWLISMINIT(CHF_CONST_FRA(U),
                    CHF_CONST_INT(iCP),
                    CHF_CONST_INT(iRegTr),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),  
                    CHF_BOX(uBox));
                            
  }
  
}


void SWLISMProblem::initialize(LevelData<FArrayBox>& a_U,
                                     Interval & a_comp)
{
}                                     


/// Calculate problem specific source terms
/**
 */
void SWLISMProblem::explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                                  BaseFab<int> & a_REG,
                               const Real      & a_dt,                               
                               const Box       & a_box)
{  
  int iBGN,iEND;
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  int iRhoN  = eqSys->densityIndexCons(1);
  int iFluids = nFluids();
  
  // Charge exchange source terms for multifluid models
  if (iFluids>= 2)
  {
    FArrayBox Wcartesian;
    FArrayBox & W = ( ((CoordinateSystem == CoordinateSystemHandler::CS_Cartesian   ) ||
                       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric) ||
                       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical )) ? (FArrayBox&)a_W : Wcartesian );

    if  ((CoordinateSystem == CoordinateSystemHandler::CS_Polar      ) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Spherical  ))
    {
      Wcartesian.define(a_box, a_W.nComp());
      Wcartesian.copy(a_W);
      m_csh->transCurvVectToCartesian(Wcartesian,a_box,m_level);
    }

    switch (iFluids)
    {
      case 2:
        {
          PickupIons* pPI = eqSys->getPickupIons();

          if( pPI != NULL )
          {
            int iRhoPIU  = pPI->consInterval().begin();
            int iRhoPIW  = pPI->primInterval().begin();

            FORT_CHARGE_EXCHANGE_2FPI_REG( CHF_FRA(a_S),
                                           CHF_CONST_FRA(W),
                                           CHF_CONST_FIA1(a_REG,0),
                                           CHF_CONST_REAL(a_dt),
                                           CHF_CONST_INT(iRhoN),
                                           CHF_CONST_INT(iRhoPIU),
                                           CHF_CONST_INT(iRhoPIW),
                                           CHF_BOX(a_box) );
          }
          else
          {
            FORT_CHARGE_EXCHANGE_MF( CHF_FRA(a_S),
                                     CHF_CONST_FRA(W),
                                     CHF_CONST_FIA1(a_REG,0),
                                     CHF_CONST_REAL(a_dt),
                                     CHF_CONST_INT(iRhoN),
                                     CHF_CONST_INT(iFluids),
                                     CHF_BOX(a_box) );
          }
        }
        break;
      case 3: FORT_CHARGE_EXCHANGE_MF( CHF_FRA(a_S),
                                       CHF_CONST_FRA(W),
                                       CHF_CONST_FIA1(a_REG,0),
                                       CHF_CONST_REAL(a_dt),
                                       CHF_CONST_INT(iRhoN),
                                       CHF_CONST_INT(iFluids),
                                       CHF_BOX(a_box) );
              break;
      case 4: 
        {
          PickupIons* pPI = eqSys->getPickupIons();

          if( pPI != NULL )
          {
            int iRhoPIU = pPI->consInterval().begin();
            int iRhoPIW = pPI->primInterval().begin();

            FORT_CHARGE_EXCHANGE_4FPI( CHF_FRA(a_S),
                                       CHF_CONST_FRA(W),
                                       CHF_CONST_FIA1(a_REG,0),
                                       CHF_CONST_REAL(a_dt),
                                       CHF_CONST_INT(iRhoN),
                                       CHF_CONST_INT(iRhoPIU),
                                       CHF_CONST_INT(iRhoPIW),
                                       CHF_BOX(a_box) );
          }
          else
          {
            FORT_CHARGE_EXCHANGE_4F( CHF_FRA(a_S),
                                     CHF_CONST_FRA(W),
                                     CHF_CONST_FIA1(a_REG,0),
                                     CHF_CONST_REAL(a_dt),
                                     CHF_CONST_INT(iRhoN), 
                                     CHF_BOX(a_box) );
          }
        }
              break;
      case 5: FORT_CHARGE_EXCHANGE_MF( CHF_FRA(a_S),
                                       CHF_CONST_FRA(W),
                                       CHF_CONST_FIA1(a_REG,0),
                                       CHF_CONST_REAL(a_dt),
                                       CHF_CONST_INT(iRhoN),
                                       CHF_CONST_INT(iFluids),
                                       CHF_BOX(a_box) );
              break;
    }

     // Add charge exchnage source terms for plasma
    iBGN   = UMOMX;
    iEND   = UENG;
    if (m_photoionize == true)
    {
      iBGN = URHO;
      if  ((CoordinateSystem == CoordinateSystemHandler::CS_Polar      ) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Spherical  ))
      {
        FORT_PHOTOIONIZE_SPH(
          CHF_FRA(a_S),
          CHF_CONST_FRA(W),
          CHF_CONST_FIA1(a_REG,0),
          CHF_CONST_REAL(a_dt),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(iFluids),
          CHF_BOX(a_box),
          CHF_CONST_INT(m_level));
      }
    }
//    FORT_ADDSOURCES( CHF_FRA(a_U),
//                     CHF_CONST_FRA(a_S),
//                     CHF_CONST_INT(iBGN),
//                     CHF_CONST_INT(iEND),
//                     CHF_BOX(a_box) );                                          
//    // Add charge exchnage source terms for neutrals
//    iBGN = eqSys->densityIndexCons(1);

    iEND = eqSys->densityIndexCons(iFluids-1)+UENG;

    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );
                     
  } else
  if (m_const_H == true)
  {
    Real Wnet[WNUM_E] = {m_netN/m_lismN, 1.0, 0.0, 0.0, eos_k*m_netN*m_lismT/(m_lismN*eos_mp*m_lismV*m_lismV)};
    int iWNUM_E = WNUM_E;
          
    
    FArrayBox Wcartesian;
    FArrayBox & W = ( ((CoordinateSystem == CoordinateSystemHandler::CS_Cartesian   ) ||
                       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric) ||
                       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical )) ? (FArrayBox&)a_W : Wcartesian );

    if  ((CoordinateSystem == CoordinateSystemHandler::CS_Polar      ) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) ||
         (CoordinateSystem == CoordinateSystemHandler::CS_Spherical  ))
    {
      Wcartesian.define(a_box, a_W.nComp());
      Wcartesian.copy(a_W);
      m_csh->transCurvVectToCartesian(Wcartesian,a_box,m_level);
      
      FORT_CHARGE_EXCHANGE_EXPSPH_CONSTH(
          CHF_FRA(a_S),
          CHF_CONST_FRA(W),        
          CHF_CONST_R1D(Wnet,iWNUM_E),  
          CHF_CONST_REAL(a_dt),          
          CHF_BOX(a_box),
          CHF_CONST_INT(m_level));
    } else
    {
      FORT_CHARGE_EXCHANGE_CONSTH(
          CHF_FRA(a_S),
          CHF_CONST_FRA(W),      
          CHF_CONST_R1D(Wnet,iWNUM_E),    
          CHF_CONST_REAL(a_dt),          
          CHF_BOX(a_box),
          CHF_CONST_INT(m_level));      
    }
    
    
    
    iBGN   = UMOMX;
    iEND   = UENG;
    
          
    FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );
  }
}
                                     
                                              // Problem specific postprocessing
void SWLISMProblem::postprocessing(       FArrayBox & a_U,
                                    const FArrayBox & a_W,
                                    const Real      & a_dt,
                                    const Real      & a_time,                                    
                                    const Box       & a_box       )
{
  Real dx = m_csh->dx(0,m_level);
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);  
  int iCP = eqSys->correctionPotentialIndex();
  int iRegTr = -1;
  if  (m_region_tracer)
      iRegTr =  (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);  
      
  if (iRegTr > 0)
  {
    FORT_REGTRACER_REINIT( CHF_FRA(a_U),                   
                      CHF_CONST_INT(iRegTr),
                      CHF_BOX(a_box) );
  }

  
  FORT_SWLISMREINIT(CHF_CONST_FRA(a_U),
                    CHF_CONST_INT(iCP),                    
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),  
                    CHF_BOX(a_box));   
                    
  if (m_subproblem == SW_TURB)
  {
    // Different "turbulent" regimes
  
    //FORT_SWLISMREINIT_TURB(CHF_CONST_FRA(a_U),
    //                  CHF_CONST_REAL(dx),
    //                  CHF_BOX(a_box));
       
    //FORT_SW_TURBR0(
    //       CHF_FRA(a_U),
    //       CHF_CONST_FRA(m_Phi),
    //       CHF_CONST_REAL(dx),
    //       CHF_BOX(a_box));
    
    BaseFab<int>  Region(a_box,1);  
    
    FORT_DEFINE_REGIONS_2F( CHF_CONST_FRA(a_W),
                            CHF_FIA1(Region,0),
                            CHF_CONST_REAL(dx),
                            CHF_BOX(a_box) );
                    
    //FORT_REG3_TURB1( CHF_FRA(a_U),
    //           CHF_CONST_FIA1(Region,0),     
    //           CHF_BOX(a_box));
               
    FORT_REG3_TURB2( CHF_FRA(a_U),
               CHF_CONST_FRA(m_Phi),
               CHF_CONST_FIA1(Region,0),     
               CHF_BOX(a_box));  
  }
                 
}


void SWLISMProblem::postTimeStep(LevelData<FArrayBox>&  a_U)
{
  if (m_subproblem == SW_TURB)
  {
    if (m_dTurbLISMTimeStep > 0) 
    {
      if( m_dTime >= m_dNextTurbLISMTime )
      {
        m_dNextTurbLISMTime += m_dTurbLISMTimeStep;      
      }
      if( m_dTime >= m_dNextTurbLISMTime )
      {
        m_dNextTurbLISMTime = m_dTime + m_dTurbLISMTimeStep;      
      }
    }
    
    if (m_dTurbSWTimeStep > 0) 
    {
      if( m_dTime >= m_dNextTurbSWTime )
      {
        m_dNextTurbSWTime += m_dTurbSWTimeStep;      
      }
      if( m_dTime >= m_dNextTurbSWTime )
      {
        m_dNextTurbSWTime = m_dTime + m_dTurbSWTimeStep;      
      }
    }
  }

  int fluids = nFluids();
  
  // Pogorelov's (passthrough) boundary conditions for neutrals  
  if (fluids > 1)
  { 
    Real dx = m_csh->dx(0,m_level);
    
    DataIterator dit1,dit2;    
    
    //for(dit1 = m_lstrips.dataIterator(); dit1.ok(); ++dit1) m_lstrips[dit1].setVal(-1.0,0);
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
    Interval intnH1(eqSys->densityIndexCons(1),eqSys->densityIndexCons(1) + UNUM_E - 1),intnH4;
    //int urho1,ueng1;
    //urho1 = eqSys->densityIndexCons(1);
    //ueng1 = urho1+UNUM_E-1;    
    a_U.copyTo(intnH1, m_lstrips, Interval(0,UNUM_E-1));
    a_U.copyTo(intnH1, m_ladjstr, Interval(0,UNUM_E-1));
    if (fluids>=5)
    {
      intnH4.define(eqSys->densityIndexCons(4), eqSys->densityIndexCons(4) + UNUM_E - 1);      
      a_U.copyTo(intnH4, m_lstrips, Interval(UNUM_E,2*UNUM_E-1));
      a_U.copyTo(intnH4, m_ladjstr, Interval(UNUM_E,2*UNUM_E-1));      
    }
          
    dit1 = m_lstrips.dataIterator();    
    dit2 = m_ladjstr.dataIterator();          
    for(; dit1.ok(); ++dit1,++dit2)
    {
      FArrayBox& from = m_lstrips[dit1];
      FArrayBox& to   = m_ladjstr[dit2];      
                                           
      FORT_NEUTRALS_SUNBC(
        CHF_CONST_FRA(from),
        CHF_FRA(to),       
        CHF_CONST_INTVECT(m_sunIJK),                          
        CHF_CONST_REAL(dx)); 
              
    }    
    m_ladjstr.copyTo(Interval(0,UNUM_E-1), a_U, intnH1);
    if (fluids>=5)
      m_ladjstr.copyTo(Interval(UNUM_E,2*UNUM_E-1), a_U, intnH4);
              
  }
  
}


void SWLISMProblem::defineRegions( const FArrayBox    & a_W,
                                            FArrayBox    & a_S,
                                            BaseFab<int> & a_R,
                                      const Box          & a_box)
{
  a_R.resize( a_box, 1 );
  
  Real dx = m_csh->dx(0,m_level);
    
  if  (m_region_tracer)
  {
    int iRegTr =  (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);  
    FORT_DEFINE_REGIONS_TRACER( CHF_CONST_FRA(a_W),
                          CHF_FIA1(a_R,0),
                          CHF_CONST_INT(iRegTr),
                          CHF_CONST_REAL(dx),
                          CHF_BOX(a_box) );
  }
  else
    FORT_DEFINE_REGIONS_2F( CHF_CONST_FRA(a_W),
                          CHF_FIA1(a_R,0),                          
                          CHF_CONST_REAL(dx),
                          CHF_BOX(a_box) );

}


// Number additional variables for writing to plot file
int SWLISMProblem::numPlotVars()
{
  int numPlotVars = 1;                         // Plasma temperature
  if  (m_output_region == true) numPlotVars++; // Plasma temperature + region
  return numPlotVars; 
}
  
// Names of the additional variables for writing to plot file  
Vector<string> SWLISMProblem::plotNames()
{
  Vector<string> retval;  
   
  if (numPlotVars() >= 1)
    retval.push_back("T");
  if (numPlotVars() >= 2)  
    retval.push_back("region");
  if (numPlotVars() > 2)  
  {
    retval.push_back("srho");
    retval.push_back("x-smom");
    retval.push_back("y-smom");
    retval.push_back("seng");    
  }    
        
  return retval;    
}
  
// Calculates variables for plotting using primitive variables
void SWLISMProblem::calcPlotVars(FArrayBox& a_vars,
                           int              a_comp,
                           const FArrayBox& a_W,
                           const Box&       a_box)
{
  if (numPlotVars() == 0) return;  
  
  BaseFab<int> Region;
  FArrayBox dummyFab;      
  defineRegions(a_W, dummyFab, Region, a_box);
  
  Real coeff = m_gamma*m_lismM*m_lismM*m_lismT;
      
  int iT      = a_comp+0;
  int iRegion = a_comp+1;
  
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    Real T = coeff*a_W(iv,WPRES)/a_W(iv,WRHO);
    if (numPlotVars() >= 1)  a_vars.set(iv, iT,      T);    
    if (numPlotVars() >= 2)  a_vars.set(iv, iRegion, (Real)(Region(iv, 0)) );    
  }
    
    
  if (numPlotVars() > 2)  
  {
    int iFluids = nFluids();
    Real dt = 1.0;
    FArrayBox S( a_box, a_W.nComp() );
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    int iRhoN  = eqSys->densityIndexPrim(1);  
    
    if (iFluids == 2)
      FORT_CHARGE_EXCHANGE_2F( CHF_FRA(S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_FIA1(Region,0),
                               CHF_CONST_REAL(dt),  
                               CHF_CONST_INT(iRhoN),
                               CHF_BOX(a_box) );
    if (iFluids == 4)
      FORT_CHARGE_EXCHANGE_4F( CHF_FRA(S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_FIA1(Region,0),
                               CHF_CONST_REAL(dt),
                               CHF_CONST_INT(iRhoN),
                               CHF_BOX(a_box) );
    
    a_vars.copy(S,0,a_comp+1,3);
    a_vars.copy(S,4,a_comp+4,1); // copy s-eng, s-momz is skipped
    
    #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(a_vars)
        );
    #endif
  }
}


void SWLISMProblem::primForPlot(       FArrayBox& a_W,
                                       const Box&       a_box)
{
  //return;
  
  Real   Bref = 1e+6*sqrt(m_lismN*eos_mp*m_lismV*m_lismV);
  
  double scalP = Bref*Bref;
  double scalV = 1.0e-5*m_lismV;
  
  //m_csh->transCartesianVectToCurv(a_W, a_box, m_level);

  a_W.mult( m_lismN, WRHO     );
  a_W.mult( scalV,   WVELX, 3 );  
  a_W.mult( scalP,   WPRES    );  
  a_W.mult( Bref,  WBX,   3 );
  
  int Fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  
  for (int iFluid = 1; iFluid < Fluids; ++iFluid)
  { 
    int iRho = eqSys->densityIndexPrim(iFluid);
    a_W.mult( m_lismN, iRho      );
    a_W.mult( scalV,   iRho+1, 3 );  
    a_W.mult( scalP,   iRho+4    );  
  }
}

//                            Return boundary condition flags for all boundaries
void SWLISMProblem::getBCFlags( eBoundaryConditions leftBC,
                                eBoundaryConditions rightBC,
                                eBoundaryConditions bottomBC,
                                eBoundaryConditions topBC,
                                eBoundaryConditions frontBC,
                                eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Continuous;
  bottomBC = BC_Axis;
  topBC    = BC_Continuous;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions SWLISMProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  case 1  : return (a_sd == Side::Lo) ? BC_Axis  : BC_Continuous;
  default : return BC_Undefined;
  }
}

/// Creates tagged cells for dynamic mesh refinement
/**
  Problem specific cells tagging
 */
void SWLISMProblem::tagCells(const FArrayBox&  a_U,
                                    const Box&        a_box,
                                          IntVectSet& a_tags)
{
  if (!m_gridBox.isEmpty())  
  {
    // Tag all m_gridBox
    BoxIterator bit(m_gridBox);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      a_tags |= iv;      
    }
  }
   
  if ((nFluids()>1) || (m_subproblem == SW_KINETIC)) // For MHD case R0 is too large to regrid it.
  if (!m_R0Box.isEmpty())  
  {
    Box R0Box(m_R0Box);
    R0Box.grow(m_R0Box.size(1)+m_R0Box.size(1)/4);
    
    // Tag all R0Box
    BoxIterator bit(R0Box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      a_tags |= iv;      
    }
  }

}                                         

/// Check geometrical/problem limitations for grid adaptation
/**
 */   
void SWLISMProblem::lockedCellsRegrid( BaseFab<int> & a_flag,
                                 const FArrayBox&  a_U,
                                 const Box&     a_box)
{
  a_flag.setVal(0);
  
  Real x;int dir;
  Real dx = m_csh->dx(0,m_level);
  
  BoxIterator bit(a_box);
  
  if( m_stop_ref > 0.0 ) 
  for (bit.begin(); bit.ok(); ++bit)
  {
    const IntVect& iv = bit();
    x = (iv[0] + 0.5)*dx;    
    if( x > m_stop_ref ) a_flag(iv,0) = 1;
  }
  
  if (m_max_levelTS <= m_level )
  {
    const Box& b = m_domain & a_U.box();
    FArrayBox W(b,8);
    
    FORT_CONSTOPRIM(
          CHF_FRA(W),
          CHF_CONST_FRA(a_U),
          CHF_BOX(b));
          
    BaseFab<int> Region;
    FArrayBox dummyFab;      
    defineRegions(W, dummyFab, Region, b);
    
    bool Reg2Present = false;
    bool Reg3Present = false;
    int  reg;
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      reg = Region(iv,0);
      if (reg == 2) Reg2Present = true;
      if (reg == 3) Reg3Present = true;      
    }
    if (Reg2Present && Reg3Present) a_flag.setVal(1);
    
    /*for (dir = 0; dir < SpaceDim; ++dir)  
    {    
      const Box bCenter = a_box & grow(m_domain,-BASISV(dir));
      const Box bLo     = a_box & adjCellLo(bCenter,dir);
      const int hasLo   = ! bLo.isEmpty();
      const Box bHi     = a_box & adjCellHi(bCenter,dir);
      const int hasHi   = ! bHi.isEmpty();
      
      FORT_LOCKEDCELLS_REGION(
         CHF_FIA1(a_flag,0),
         CHF_CONST_FIA1(Region,0),
         CHF_CONST_INT(dir),
         CHF_BOX(bLo),
         CHF_CONST_INT(hasLo),
         CHF_BOX(bHi),
         CHF_CONST_INT(hasHi),
         CHF_BOX(bCenter));     
   }*/
  } 
}

   

/// Compute problem specific dt and returns a cell with the minimum dt. 
/**
 */                               
Real SWLISMProblem::computeDt( const FArrayBox& a_U,
                               const FArrayBox& a_dt,
                               const Box&     a_box,
                               IntVect&       a_minDtCell)
{
    
  //return MultiFluidProblem::computeDt( a_U, a_dt, a_box, a_minDtCell);    
  
  if( m_R0dtIgnore <= 0.0 ) 
  {
    return MultiFluidProblem::computeDt( a_U, a_dt, a_box, a_minDtCell);    
  }

  Real D_DECL(x,y,z), dist2, R0dtIgnore2, dt;
  Real dx = m_csh->dx(0,m_level);
  
  R0dtIgnore2 = m_R0dtIgnore*m_R0dtIgnore;
  
  Real dtMin = numeric_limits<Real>::max();
  
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  {
    const IntVect& iv = bit();                      
    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]  );
      
    dist2 = D_TERM( x*x, + y*y, + z*z );

    dt = a_dt.get(iv,0);
    if ((dist2 > R0dtIgnore2 ) && (dt < dtMin))
    {
      dtMin = dt;
      a_minDtCell = iv;
    }      
  }
  
  return dtMin;
}

extern "C"
{

// Prototype for Fortran procedure fft2dmain
#define FORT_FFT2DMAIN FORTRAN_NAME( FFT2DMAIN ,fft2dmain )
void FORT_FFT2DMAIN
(
    int* const nx, 
    int* const ny, 
    Real* const dx,
    Real* const phi    
);

#define FORT_FFT2DMAIN2 FORTRAN_NAME( FFT2DMAIN2 ,fft2dmain2 )
void FORT_FFT2DMAIN2
(
    
);


}

// Calculate Phi for turbulence
void SWLISMProblem::CalculatePhi(const LevelData<FArrayBox>& a_U,
                                       PatchMHDAM& a_patch)
{
  //if (m_dTurbSWTimeStep < 0.0) return;
  
  // First, we need bounding box of the region 3.  
  Box bbReg3;   

  Real dx = m_csh->dx(0,m_level); 
  
  BaseFab<int> Region(m_domain.domainBox(),1);  

  int D_DECL(i,j,k);
  
  int reg;
  
  IntVect loVect = m_domain.domainBox().bigEnd();
  IntVect hiVect = m_domain.domainBox().smallEnd();
  
  IntVect aux_IntVect;
  
  FArrayBox dummyFAB;
  
  const DisjointBoxLayout& grids = a_U.disjointBoxLayout();  
  for (DataIterator dit = grids.dataIterator(); dit.ok(); ++dit)
  {
    // The current box
    Box curBox = grids.get(dit());
    // The current grid of conserved variables
    const FArrayBox& curU = a_U[dit()];    
    FArrayBox W(curBox,a_U.nComp());
    
    m_eqSys->stateToPrim(W, curU, curBox);
    
    defineRegions(W, dummyFAB, Region, curBox);    
    
    D_TERM(
      for (i=curBox.smallEnd(0);i<=curBox.bigEnd(0);i++),
      for (j=curBox.smallEnd(1);j<=curBox.bigEnd(1);j++),
      for (k=curBox.smallEnd(2);k<=curBox.bigEnd(2);k++))                                     
      { 
        D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
        Region.getVal(&reg,aux_IntVect,0,1);        
        if (reg == 3)
        {
          loVect.min(aux_IntVect);
          hiVect.max(aux_IntVect);          
        }
      }
  }  

  if (loVect != m_domain.domainBox().bigEnd())    
    CH_assert(loVect<hiVect);
  
  IntVect loVectGlobal,hiVectGlobal;  
  
  //MPI_Allreduce(loVect.dataPtr(),loVectGlobal.dataPtr(),CH_SPACEDIM, MPI_INT, MPI_MIN,Chombo_MPI::comm);
  //MPI_Allreduce(hiVect.dataPtr(),hiVectGlobal.dataPtr(),CH_SPACEDIM, MPI_INT, MPI_MAX,Chombo_MPI::comm);
  
  //loVect = loVectGlobal;
  //hiVect = hiVectGlobal;
  //assert(loVect<hiVect);  
        
  //int max_size = MAX(hiVect[0] - loVect[0],hiVect[1] - loVect[1]);
  //hiVect[0] = loVect[0] + max_size;
  //hiVect[1] = loVect[1] + max_size;
  
  // Using R0Turb
  //int iSun = (int)(m_sunXC/m_dx);
  //int jSun = (int)(m_sunYC/m_dx);
  int iSun = m_sunIJK[0];
  int jSun = m_sunIJK[1];
  int iR0Turb = (int)(m_R0Turb/dx) +1;
  
  loVect[0] = iSun - iR0Turb;
  loVect[1] = jSun - iR0Turb;
  
  hiVect[0] = iSun + iR0Turb;
  hiVect[1] = jSun + iR0Turb;
    
     
  bbReg3.define(loVect,hiVect);
  m_Phi.define(bbReg3,CH_SPACEDIM);// CH_SPACEDIM components for velocity disturbance.
  
  int nx = bbReg3.size(0);
  int ny = bbReg3.size(1);
  
  int curProc = procID();  
  
  Real* tmpPhi;
  
  if (curProc == 0)
  {    
    Real PhiMax;    
        
    tmpPhi = m_Phi.dataPtr(0); // x - component  
    FORT_FFT2DMAIN(&nx,&ny,&dx,tmpPhi);    
    PhiMax = MAX(fabs(m_Phi.max(0)),fabs(m_Phi.min(0)));  
    m_Phi.mult(1.0/PhiMax,0);
    
    tmpPhi = m_Phi.dataPtr(1); // y - component  
    FORT_FFT2DMAIN(&nx,&ny,&dx,tmpPhi);
    PhiMax = MAX(fabs(m_Phi.max(1)),fabs(m_Phi.min(1)));  
    m_Phi.mult(1.0/PhiMax,1);  
  }
  
  tmpPhi = m_Phi.dataPtr();

#ifdef CH_MPI  
  MPI_Bcast(tmpPhi,bbReg3.numPts()*m_Phi.nComp(),MPI_CH_REAL,0,Chombo_MPI::comm);
#endif

  static bool haswritten = false;
  if (!haswritten)
  {
    char file_name[1000];  
    sprintf(file_name,"phi%i.dat",curProc);
    FILE* tfile = OpenTecplotFile(file_name,"TITLE = \"PHI\"\n VARIABLES=\"X\" \"Y\" \"PHIX\" \"PHIY\"");
    WriteFArrayBoxToTecplotFile(tfile, m_Phi, m_Phi.box(), Interval(0,m_Phi.nComp()-1), m_csh->dx(0,m_level));
    CloseTecplotFile(tfile);
    haswritten = true;
  }

      
}

