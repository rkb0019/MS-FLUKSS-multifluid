#include <limits>

#include "BoxIterator.H"
#include "IntVectSet.H"
#include "LoadBalance.H"
#include "CH_HDF5.H"

#include "HeliosphericProblem.H"
#include "HeliosphericF_F.H"
#include "HeliosphericPlasmaBCF_F.H"
#include "SWLISMF_F.H"
#include "ChargeExchange2F_F.H"

#include "bedfordcxF_F.H"
#include "bedford4FPIREGF_F.H"

#include "LGintegrator.H"
#include "DebugF_F.H"
#include "TecplotIO.H"
#include "EosCommon.H"
#include "CSHSpherical.H"
#include "CSHSphericalF_F.H"
#include "MHDAMDefs.H"
#include "EqSysMHDMF.H"
#include "RiemannSolver.H"
#include "PatchIdealMHDF_F.H"
#include "TMBreechEtAl2008.H"
#include "HeliosphericTurbF_F.H"
#include "PITwoEquations.H"
#include "HeliosphericPIF_F.H"
#include "PatchMHDAMF_F.H"
#include "ShocksF_F.H"

//Extra includes
#include <time.h>
#include <sstream>
#include <fstream>
//Extra includes

//#define OUTPUT_DT
//#define OUTPUT_RELGRADB // output abs(grad |B|))


#define SPECIAL_MESH512

#ifdef OUTPUT_DT
#include "PatchEulerF_F.H"
#endif 

#ifdef OUTPUT_RELGRADB
#include "PatchMHDAMF_F.H"
#endif 


HeliosphericProblem::HeliosphericProblem()
 : MultiFluidProblem()
{
  m_isFortranCommonSet = false;
  m_physModel          = PP_Undefined;  

  m_adaptRmin          = 0.0;
  m_adaptRmax          = 0.0;

  m_helioAdapt         = 0;
  m_adaptXL1           = 1.0e10;
  m_adaptXL2           = 1.0e10;

  m_bSunGravity        = false;
  m_bSunHeating        = false;
  
  m_writeH5            = false;
  
  m_photoionize        = false;  
  m_const_H            = false;
  m_region_tracer      = false;
  
  m_output_vecSPH      = false;
}

// Input parameters
void HeliosphericProblem::input( ParmParse & parser, int verbosity )
{
  m_R0         = 0.1;
  m_concludeR0 = 10.0;

  m_lismN      = 0.06;
  m_lismV      = 2640000.0;
  m_lismT      = 6527.0;
  m_lismB      = 1.5e-6;
  m_lismBX     = 1.0;
  m_lismBY     = 0.0;
  m_lismBZ     = 0.0;
  m_lismUX     = -1.0;
  m_lismUY     = 0.0;
  m_lismUZ     = 0.0;
  m_netN       = 0.1;
  
  m_sunN       = 7.4;
  m_sunV       = 45000000.0;
  m_sunT       = 51100.0;
  m_sunB       = 0.0;
  m_sunTILT    = 0.0;  
  m_initR      = 80.0;
  m_subproblem = 0;
  m_dirBrNorth = 1.0;
  m_sunBmonopolar = -1;
  m_TMLIM      = 50000.0;
  m_writeDt    = 24.0;
  

  m_sunIntBCRadius = 1.0;

  int iFluids  = 1;

  m_sunZ2      = 3e13;
  m_sunSigmaC  = 0.60;
  m_sunLambda  = 0.008*1.5e+13;
  m_lismZ2     = 0.01*m_sunZ2;
  m_lismSigmaC =      m_sunSigmaC;
  m_lismLambda = 10.0*m_sunLambda;

  m_sunN_PI    = 0.001;
  m_sunT_PI    = 51000.0;
  m_lismN_PI   = 0.001;
  m_lismT_PI   = 6527.0;

  m_gamma      = 1.6666666667;

  m_verbosity  = verbosity;
  
  m_output_region = false;

  parser.query( "gamma",   m_gamma   );
  parser.query( "lismN",   m_lismN   );
  parser.query( "lismV",   m_lismV   );
  parser.query( "lismT",   m_lismT   );
  parser.query( "lismB",   m_lismB   );
  parser.query( "lismBX",  m_lismBX  );
  parser.query( "lismBY",  m_lismBY  );
  parser.query( "lismBZ",  m_lismBZ  );
  parser.query( "lismUX",  m_lismUX  );
  parser.query( "lismUY",  m_lismUY  );
  parser.query( "lismUZ",  m_lismUZ  );
  parser.query( "sunN",    m_sunN    );
  parser.query( "sunV",    m_sunV    );
  parser.query( "sunT",    m_sunT    );
  parser.query( "sunB",    m_sunB    );
  parser.query( "sunTILT", m_sunTILT );  
  parser.query( "initR",   m_initR   );  
  parser.query( "dirBrNorth",     m_dirBrNorth     );
  parser.query( "sunBmonopolar",  m_sunBmonopolar  );
  parser.query( "subproblem",     m_subproblem     );
  parser.query( "sunIntBCRadius", m_sunIntBCRadius );
  parser.query( "concludeR0",     m_concludeR0     );
  
  if (m_sunBmonopolar!=1) m_sunBmonopolar = 0;
  
  int H5write = 0;
  parser.query( "h5write",  H5write);
  m_writeH5 = (H5write == 1);
      
  if (m_writeH5 == true)  
  {   
    parser.query("h5write_file_name", m_writeH5File);
    parser.query("h5write_dt", m_writeDt);   // hours
    parser.query("h5write_start", m_writeT0);   
    
    int writeRotating = 1;
    parser.query("h5write_rotating_frame", writeRotating);        
    m_writeRotating = (writeRotating == 1);
    
    m_writeDt*=(m_lismV*3600.0)/eos_AU; // dimensionless    
  }

  Real rho        = m_lismN*eos_mp;
  Real pref       = rho*m_lismV*m_lismV;
  Real p          = 2.0*eos_k*m_lismN*m_lismT;
  Real lismP      = p/pref;
  m_lismM         = 1.0/sqrt( m_gamma*lismP );
  
  if (m_dirBrNorth>0.0) m_dirBrNorth = 1.0; else m_dirBrNorth = -1.0;

  parser.query( "fluids",  iFluids   );

  switch( iFluids ){
  case 2  : m_physModel  = PP_2FluidPM; break;
  case 3  : m_physModel  = PP_3FluidPM; break;
  case 4  : m_physModel  = PP_4FluidPM; break;
  case 5  : m_physModel  = PP_5FluidPM; break;
  default : m_physModel  = PP_MHDPM;
  }

  D_TERM(parser.query( "XC",      m_sunXYZ[0]    );,
         parser.query( "YC",      m_sunXYZ[1]    );,
         parser.query( "ZC",      m_sunXYZ[2]    ););
         
  parser.query( "R0",      m_R0    );

  parser.query( "adaptR1", m_adaptRmin );
  parser.query( "adaptR2", m_adaptRmax );

//                       Default value of m_helioAdapt depends on the subproblem
  switch( m_subproblem ) {
  case HPBC_LONGTAIL :
    m_helioAdapt = HPAD_LONGTAIL;
    break;
  case HPBC_HPINSTAB :
    m_helioAdapt = HPAD_HPINSTAB;
    break;
  case HPBC_SOLARCYCLE :
    m_helioAdapt = HPAD_TAILDETAIL;
    break;
  default :
    m_helioAdapt = HPAD_CARTESIAN;
  }

  parser.query( "helioAdaptation", m_helioAdapt );

//  DEfault value of m_adaptXL1 and m_adaptXL2 corresponds to tagCellsTailDetail
  m_adaptXL1 =-2000.0;
  m_adaptXL2 =- 800.0;
  parser.query( "adaptXL1", m_adaptXL1 );
  parser.query( "adaptXL2", m_adaptXL2 );

  int i;
  
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
  
  m_R0dtIgnore = m_R0 - 2.0*dx;
  
  m_R0dtIgnore = -1.0;

//        Set boundaries of solar wind spreading in Z direction for Cartesian CS
  Real safeDZ = 0.0;
  parser.query("RegSafeDZ",safeDZ);
  parser.query("RegionSafetyDZ",safeDZ);

  Real minZ =-m_sunXYZ[2];
  Real maxZ = dx*numCells[2] - m_sunXYZ[2];
  m_RegSafeZTop  = maxZ - safeDZ;
  m_RegSafeZBot  = minZ + safeDZ;
 
  int photoionize = 0;
  parser.query( "photoionize",   photoionize   );
  m_photoionize = (photoionize==1);
  
  int const_H = 0;
  parser.query( "const_H",   const_H);
  m_const_H = (const_H==1);
  
  if ((iFluids > 1) && (m_const_H))
  {
    MayDay::Error("constant neutrals should be used with the mhd model only");
  }  
  if ((nFluids()>1) || (m_const_H))
  {
    if (!parser.contains("netN")) MayDay::Error("mhdam.netN must be defined for multifluid calculations");    
  }
  parser.query( "netN",    m_netN  );  
   
  
  int iSunGravity = 0;
  parser.query( "sun_gravity",  iSunGravity);
  m_bSunGravity = (iSunGravity == 1);
  
  int iSunHeating = 0;
  parser.query( "sun_heating",  iSunHeating);    
  m_bSunHeating = (iSunHeating == 1);
  
  parser.query( "TMLIM",        m_TMLIM);
  
  int nRegionTracer = 1; // -1;
  parser.query("region_tracer",nRegionTracer);  
  m_region_tracer = (nRegionTracer > 0);
  
  int output_region = 1; // -1;
  parser.query( "output_region", output_region);
  m_output_region = output_region > 0;

  parser.query( "lismZ2",       m_lismZ2     );
  parser.query( "lismSigmaC",   m_lismSigmaC );
  parser.query( "lismLambda",   m_lismLambda );
  parser.query( "sunZ2",        m_sunZ2      );
  parser.query( "sunSigmaC",    m_sunSigmaC  );
  parser.query( "sunLambda",    m_sunLambda  );

  parser.query( "lismN_PI",     m_lismN_PI   );
  parser.query( "lismT_PI",     m_lismT_PI   );
  parser.query( "sunN_PI",      m_sunN_PI    );
  parser.query( "sunT_PI",      m_sunT_PI    );

  m_Bref = 1e+6*sqrt(m_lismN*eos_mp*m_lismV*m_lismV);
  
  std::string problemString;    
  parser.query("problem",problemString);
  if (problemString == "helioKinetic")    
  {
    m_subproblem = HPBC_KINETIC;
    if (iFluids != 1) MayDay::Error("Number of fluids should be 1 for the kinetic problem");    
  }
  
  int output_vecSPH = 0;
  parser.query ("output_vec_sph",output_vecSPH);
  m_output_vecSPH = (output_vecSPH == 1);
     
                                                         // Print the parameters
  if (( verbosity >= 2 ) || (procID() == 0))
  {
    pout() << "The 3D heliosperic problem input:" << endl;
    pout() << "subproblem= " << m_subproblem << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "lismN     = " << m_lismN     << endl;
    pout() << "lismV     = " << m_lismV     << endl;
    pout() << "lismB     = " << m_lismB     << endl;
    pout() << "lismBX    = " << m_lismBX    << endl;
    pout() << "lismBY    = " << m_lismBY    << endl;
    pout() << "lismBZ    = " << m_lismBZ    << endl;
    pout() << "sunN      = " << m_sunN      << endl;
    pout() << "sunV      = " << m_sunV      << endl;
    pout() << "sunT      = " << m_sunT      << endl;
    pout() << "sunB      = " << m_sunB      << endl;
    pout() << "sunIntBCRadius = " << m_sunIntBCRadius << endl;
    pout() << "initR     = " << m_initR     << endl;
    pout() << "netN      = " << m_netN      << endl;
    D_TERM(
      pout() << "XC        = " << m_sunXYZ[0]        << endl;,
      pout() << "YC        = " << m_sunXYZ[1]        << endl;,   
      pout() << "ZC        = " << m_sunXYZ[2]        << endl;); 
    pout() << "R0        = " << m_R0          << endl;
    pout() << "adaptR1   = " << m_adaptRmin   << endl;
    pout() << "adaptR2   = " << m_adaptRmax   << endl;

    pout() << "helioAdpt = " << m_helioAdapt  << endl;
    pout() << "adaptXL1  = " << m_adaptXL1    << endl;
    pout() << "adaptXL2  = " << m_adaptXL2    << endl;

	pout() << "RegSafeZT = " << m_RegSafeZTop << endl;
    pout() << "RegSafeZB = " << m_RegSafeZBot << endl;
  }

  setFortranCommon(  );
  
  if ((m_subproblem == HPBC_CIR) || (m_subproblem == HPBC_SOLARCYCLE))
  {
    Real Vfast=m_sunV,Vslow=m_sunV,Tfast=m_sunT,Tslow=m_sunT,Nfast=m_sunN,Nslow=m_sunN,fs_a=1.0,fs_b=45.0;
    Real tilt_min = 10.0, tilt_max = 70.0;
    Real slow_min = 15.0 ,slow_max = 80.0;
    
    parser.query( "Vfast", Vfast);
    parser.query( "Vslow", Vslow);
    parser.query( "Tfast", Tfast);
    parser.query( "Tslow", Tslow);
    parser.query( "Nfast", Nfast);
    parser.query( "Nslow", Nslow);
    parser.query( "fs_a",  fs_a);
    parser.query( "fs_b",  fs_b);
    
    parser.query( "tilt_min",  tilt_min);
    parser.query( "tilt_max",  tilt_max);
    parser.query( "slow_min",  slow_min);
    parser.query( "slow_max",  slow_max);
    
    int iCoordSys = 0;
    parser.query( "coord_sys", iCoordSys );
    CoordinateSystemHandler::eCoordinateSystem eCS = (CoordinateSystemHandler::eCoordinateSystem)(iCoordSys);
    if (eCS == CoordinateSystemHandler::CS_Cartesian)
    {
      domainLength = 1;
    }
    
  
    FORT_SETFASTSLOWSW(
      CHF_CONST_REAL(domainLength),
      CHF_CONST_REAL(Vfast),
      CHF_CONST_REAL(Vslow),
      CHF_CONST_REAL(Tfast),
      CHF_CONST_REAL(Tslow),
      CHF_CONST_REAL(Nfast),
      CHF_CONST_REAL(Nslow),
      CHF_CONST_REAL(fs_a),
      CHF_CONST_REAL(fs_b),
      CHF_CONST_REAL(tilt_min),
      CHF_CONST_REAL(tilt_max),
      CHF_CONST_REAL(slow_min),
      CHF_CONST_REAL(slow_max));
    
  }
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void HeliosphericProblem::setFortranCommon( )
{
  CH_assert(m_isFortranCommonSet == false);

  if( m_verbosity >= 2 ) pout() << "Enter into setFortranCommon" << endl;

  FORT_SETHELIOS( CHF_CONST_REAL( m_gamma   ),
                  CHF_CONST_REAL( m_lismN   ),
                  CHF_CONST_REAL( m_lismV   ),
                  CHF_CONST_REAL( m_lismUX  ),
                  CHF_CONST_REAL( m_lismUY  ),
                  CHF_CONST_REAL( m_lismUZ  ),
                  CHF_CONST_REAL( m_lismT   ),
                  CHF_CONST_REAL( m_lismB   ),
                  CHF_CONST_REAL( m_lismBX  ),
                  CHF_CONST_REAL( m_lismBY  ),
                  CHF_CONST_REAL( m_lismBZ  ),
                  CHF_CONST_REAL( m_sunXYZ[0]   ),
                  CHF_CONST_REAL( m_sunXYZ[1]   ),
                  CHF_CONST_REAL( m_sunXYZ[2]   ),
                  CHF_CONST_REAL( m_R0    ),
                  CHF_CONST_REAL( m_sunN    ),
                  CHF_CONST_REAL( m_sunV    ),
                  CHF_CONST_REAL( m_sunT    ),
                  CHF_CONST_REAL( m_sunB    ),
                  CHF_CONST_REAL( m_sunTILT ),
                  CHF_CONST_REAL( m_sunIntBCRadius ),
                  CHF_CONST_REAL( m_initR   ),
                  CHF_CONST_REAL( m_netN    ),
                  CHF_CONST_REAL( m_dirBrNorth),
                  CHF_CONST_REAL( m_TMLIM),
                  CHF_CONST_REAL( m_RegSafeZTop),
                  CHF_CONST_REAL( m_RegSafeZBot),
                  CHF_CONST_INT(  m_sunBmonopolar)  );

  if( nFluids() >= 1 )
  {
    Real scaleLen  = 1.5e+13;

    FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( m_lismV  ),
                           CHF_CONST_REAL( scaleLen ),
                           CHF_CONST_REAL( m_lismN  ) );
  }

  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void HeliosphericProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this object.
 */
void HeliosphericProblem::define(const ProblemDomain& a_domain,                                  
                                 const int            a_level)
{
  MultiFluidProblem::define(a_domain, a_level);
  
  int i;
  int fluids = nFluids();
      
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)
  {
    IntVect LCorner, UCorner;
    
    Real dx = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? m_csh->dx(0,m_level) : 1);
    
    D_TERM(
      m_sunIJK[0] = (int)floor(m_sunXYZ[0]/dx+0.5);,
      m_sunIJK[1] = (int)floor(m_sunXYZ[1]/dx+0.5);,
      m_sunIJK[2] = (int)floor(m_sunXYZ[2]/dx+0.5););  
    
    D_TERM(  
      CH_assert(fabs(m_sunIJK[0]*dx-m_sunXYZ[0])<1e-6);,
      CH_assert(fabs(m_sunIJK[1]*dx-m_sunXYZ[1])<1e-6);,
      CH_assert(fabs(m_sunIJK[2]*dx-m_sunXYZ[2])<1e-6););
      
    IntVect SunPos(m_sunIJK);    
    
    // Defining m_R0Box
    Real R0 = m_R0;
    int iR0 = 1 + (int)floor(R0/dx+0.5); 
    int gridDim = iR0;
      
    LCorner = SunPos; UCorner = SunPos;
    LCorner-=gridDim;
    UCorner+=gridDim-1;          
    
    m_R0dtIgnore = m_R0 - 2.0*dx;
    
    m_R0Box.define(LCorner,UCorner);        
    
    // Define data for path-through BC 
    if (fluids > 1)
    {
      Box strip_a,strip_b;// '_a' means after the Sun, '_b' means before the Sun
      Vector<Box> vstrips_a, vstrips_b;
      for (i=0; i<iR0; i++)
      {
        strip_a.define(IntVect( D_DECL(SunPos[0]+i,    SunPos[1]-iR0-1, SunPos[2]-iR0-1)),
                       IntVect( D_DECL(SunPos[0]+i,    SunPos[1]+iR0  , SunPos[2]+iR0  )));
                       
        strip_b.define(IntVect( D_DECL(SunPos[0]-i-1,  SunPos[1]-iR0-1, SunPos[2]-iR0-1)),
                       IntVect( D_DECL(SunPos[0]-i-1,  SunPos[1]+iR0  , SunPos[2]+iR0  )));
                       
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
    
    if (m_subproblem == HPBC_KINETIC)
    {
      //if (m_eqSys->numTrackingSurfaces()>0)
      //{
      //  MayDay::Error("kinetic problems must not be used with level set");
      //}
          
      sprintf(buf, "kinetic_grid%ibb", a_level+1);
      if (parser.contains(buf))
      {
        parser.queryarr( buf, tmpVect, 0, 2*CH_SPACEDIM);    
        RealVect rLo(D_DECL(tmpVect[0], tmpVect[1], tmpVect[2]));
        RealVect rHi(D_DECL(tmpVect[3], tmpVect[4], tmpVect[5]));
                
        rLo += m_sunXYZ;    
        rHi += m_sunXYZ;    
        
        IntVect iLo(D_DECL((int)(rLo[0]/dx), (int)(rLo[1]/dx), (int)(rLo[2]/dx)));
        IntVect iHi(D_DECL((int)(rHi[0]/dx), (int)(rHi[1]/dx), (int)(rHi[2]/dx)));
        
        m_kineticBox.define(iLo,iHi);        
      } 
    }
      
    sprintf(buf, "level_grid%ibb", a_level+1);
    if (parser.contains(buf))
    {
      parser.queryarr( buf, tmpVect, 0, 2*CH_SPACEDIM);    
      RealVect rLo(D_DECL(tmpVect[0], tmpVect[1], tmpVect[2]));
      RealVect rHi(D_DECL(tmpVect[3], tmpVect[4], tmpVect[5]));
              
      rLo += m_sunXYZ;    
      rHi += m_sunXYZ;    
      
      IntVect iLo(D_DECL((int)(rLo[0]/dx), (int)(rLo[1]/dx), (int)(rLo[2]/dx)));
      IntVect iHi(D_DECL((int)(rHi[0]/dx), (int)(rHi[1]/dx), (int)(rHi[2]/dx)));
      
      m_R0Box.define(iLo,iHi);        
    }
  }  
   
  // Define data for path-through BC (spherical coordinate system)
  if ((fluids > 1)&&(m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical))
  {
    const Box & dBox = m_domain.domainBox();   

    Vector<Box> vstrips, vadjstrips;
    Vector<int> procs;
    Box b;

    b.define( IntVect(D_DECL(dBox.smallEnd()[0], dBox.smallEnd()[1], dBox.smallEnd()[2])), 
              IntVect(D_DECL(dBox.smallEnd()[0], dBox.size(1)/4 - 1, dBox.bigEnd()[2])));               
    vstrips.push_back(b);
    
    b.define( IntVect(D_DECL(dBox.smallEnd()[0]-1,   dBox.smallEnd()[1],     dBox.smallEnd()[2])), 
              IntVect(D_DECL(dBox.smallEnd()[0]-1, 2*dBox.size(1)/4 - 1, dBox.bigEnd()[2])));               
    vadjstrips.push_back(b);
    
    b.define( IntVect(D_DECL(dBox.smallEnd()[0]-1, 2*dBox.size(1)/4,     dBox.smallEnd()[2])), 
              IntVect(D_DECL(dBox.smallEnd()[0]-1,   dBox.bigEnd()[1], dBox.bigEnd()[2])));               
    vadjstrips.push_back(b);
            
    b.define( IntVect(D_DECL(dBox.smallEnd()[0], 3*dBox.size(1)/4  , dBox.smallEnd()[2])), 
              IntVect(D_DECL(dBox.smallEnd()[0],   dBox.bigEnd()[1], dBox.bigEnd()[2]))); 
    vstrips.push_back(b);        
                                  
    procs.push_back(0);  
    procs.push_back(0);  
    
    DisjointBoxLayout strips(vstrips, procs/*, m_domain*/);
    DisjointBoxLayout adjacent(vadjstrips, procs/*, m_domain*/);   

    int numVars = (fluids == 5 ? 2*UNUM_E : UNUM_E);     
    
    m_lstrips.define(strips,   numVars);
    m_ladjstr.define(adjacent, numVars);                      
  }
  
  if (m_writeH5 == true)  
  {
    prepareForH5Writing();
  }
  
  m_ls_indices.m_iHCS   = -1;     
  m_ls_indices.m_iHCSb  = -1;    
  m_ls_indices.m_iRegTr = -1;   
  
  if (m_eqSys->numTrackingSurfaces()>0)
  {
    if (m_region_tracer == true)
    {
      if (m_eqSys->numTrackingSurfaces()>1) m_ls_indices.m_iHCS  = m_eqSys->lvlsStateInterval().begin();
      if (m_eqSys->numTrackingSurfaces()>2) m_ls_indices.m_iHCSb = m_ls_indices.m_iHCS + 1;
      m_ls_indices.m_iRegTr = m_eqSys->lvlsStateInterval().end();
    } else
    {
      if (m_eqSys->numTrackingSurfaces()>0) m_ls_indices.m_iHCS  = m_eqSys->lvlsStateInterval().begin();
      if (m_eqSys->numTrackingSurfaces()>1) m_ls_indices.m_iHCSb = m_ls_indices.m_iHCS + 1;
    }
  }
}
                                               // Set the Equation System flag
void HeliosphericProblem::setEquationSystem( EquationSystem* a_eqSys )
{
  if( a_eqSys != NULL )
  {
    m_eqSys = a_eqSys;

    TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
    if( pTurbMod != NULL )
    {
      if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
      {
        FORT_SETHELIOS_TM( CHF_CONST_REAL( m_lismZ2     ),
                           CHF_CONST_REAL( m_lismLambda ),
                           CHF_CONST_REAL( m_lismSigmaC ),
                           CHF_CONST_REAL( m_sunZ2      ),
                           CHF_CONST_REAL( m_sunLambda  ),
                           CHF_CONST_REAL( m_sunSigmaC  ) );

        const TMBreechEtAl2008* pBreech  = dynamic_cast<const TMBreechEtAl2008*>(pTurbMod);

        Real scaleLen  = 1.5e+13;
        pBreech->setScales( scaleLen, m_lismV );
      }
    }

    PickupIons * pPickUp = m_eqSys->getPickupIons();
    if( pPickUp != NULL )
    {
      if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
      {
        FORT_SETHELIOS_PI( CHF_CONST_REAL( m_lismN_PI ),
                           CHF_CONST_REAL( m_lismT_PI ),
                           CHF_CONST_REAL( m_sunN_PI  ),
                           CHF_CONST_REAL( m_sunT_PI  ) );

         //RKB          
         // Load Table for applying hybrid BC for PUI across TS
         FORT_LOAD_TABLE_BC();

      }
    }
  }
}

void HeliosphericProblem::defineMesh(const ProblemDomain & a_prob_domain,
                                     const Vector<Real>  & a_domainBox)
{
  if (m_csh->coordinateSystem() != CoordinateSystemHandler::CS_Spherical) return;    
  
  const Box & domain = a_prob_domain.domainBox();   
    
      
  int i; IntVect iv_off;Box boxD;
      
  ParmParse parser("mesh");      
      
  //return;
    
  if (parser.contains("hs_i"))
  {
    // Fill dr    
      
    Real * data_r = new Real[domain.size(0)+1];  
    Real R0 = a_domainBox[0];    
        
    int iHP  = 3*domain.size(0)/4;
    Real rHP = 200.0;
    Real cHP = 0.1;    
    Real c2 = 4.3;
    
    parser.query( "hs_i",   iHP  ); // How many points in 'r' direction before the HP
    parser.query( "hs_r",   rHP  ); // HP radius
    parser.query( "hs_c1",  cHP  ); // Expansion coefficient
    parser.query( "hs_c2",  c2  );  // Expansion coefficient after the HP
    
    int imax = iHP;
            
    // Fill first iHP cells
    for (i=0;i<iHP;i++)
    {
      data_r[i] = R0+(rHP-R0)*(exp(cHP*i/imax)-1)/(exp(cHP)-1);
    }
    
    // Fill rest cells      
    imax    = domain.size(0)-iHP;  
    for (i=iHP;i<=domain.size(0);i++)
    {
      data_r[i] = rHP+(a_domainBox[SpaceDim]-rHP)*(exp(c2*(i-iHP)/imax)-1)/(exp(c2)-1);
    }
    
    iv_off = IntVect::Zero;
    iv_off[0] = 1;
          
    boxD.define(a_prob_domain.domainBox().smallEnd()*iv_off, 
             a_prob_domain.domainBox().bigEnd()  *iv_off);
    FArrayBox dr(boxD,1);
             
    BoxIterator bit(boxD);
    for (bit.begin(); bit.ok(); ++bit)
    { 
      IntVect iv = bit();
                              
      dr.set(iv, 0, data_r[iv[0]+1]-data_r[iv[0]]);         
    }
    
    
    if ((procID() == 0) && (0))
    {
      for (i=0;i<=domain.size(0);i++)
        pout() << data_r[i] << endl;
    }
    
    delete[] data_r;
    
    m_csh->setGridSpacing(dr,0);
  }

  //return;
  
  // Fill theta
  Real dtheta_approx = d_PI/domain.size(2);
  Real nearZangle  = 14.0/180.0*d_PI;
  Real alpha       = 1.11;
       alpha       = 1.21; // Izmodenov test
  int  N           = (int) ceil(log((dtheta_approx-nearZangle*(1-alpha))/dtheta_approx)/log(alpha)-1);
  N = ((N+1)/2)*2; // N should be even
  
  Real * dtheta_arr = new Real[domain.size(2)];
  
  Real dtheta     = dtheta_approx;
  Real dtheta_sum = 0.0;
  for (i=0;i<N;i++)
  {
    int i1 = N-1-i;
    dtheta_arr[i1] = dtheta;
    
    int i2 = domain.size(2)-1-i1;
    dtheta_arr[i2] = dtheta;
            
    dtheta_sum += dtheta;
    dtheta     *= alpha;
  }
  
  dtheta = (d_PI_2-dtheta_sum)/(domain.size(2)/2-N);
  
  for (i=N;i<domain.size(2)/2;i++)
  {    
    dtheta_arr[i] = dtheta;        
    dtheta_arr[domain.size(2)-1-i] = dtheta_arr[i];                
  }
  
  iv_off    = IntVect::Zero;
  iv_off[2] = 1;
        
  boxD.define(a_prob_domain.domainBox().smallEnd()*iv_off, 
           a_prob_domain.domainBox().bigEnd()  *iv_off);
  FArrayBox dthetaFAB(boxD,1,dtheta_arr);
  
  m_csh->setGridSpacing(dthetaFAB,2);
  
  delete[] dtheta_arr;
      
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* HeliosphericProblem::new_PhysProblem()
{
  HeliosphericProblem* retval = new HeliosphericProblem();
  
  retval->copy_PhysProblem(this);  

  return static_cast<PhysProblem*>(retval);
}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void HeliosphericProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  const HeliosphericProblem* PP = dynamic_cast<const HeliosphericProblem*>(a_PP);

  if (PP == NULL) MayDay::Error("HeliosphericProblem::copy_PhysProblem. Wrong argument");

  if (PP->m_isFortranCommonSet == true)
  {
    this->m_sunXYZ     = PP->m_sunXYZ; 
    this->m_R0         = PP->m_R0;    
    this->m_concludeR0 = PP->m_concludeR0;

    this->m_lismN      = PP->m_lismN;
    this->m_lismM      = PP->m_lismM;
    this->m_lismV      = PP->m_lismV;
    this->m_lismT      = PP->m_lismT;
    this->m_lismB      = PP->m_lismB;
    this->m_lismBX     = PP->m_lismBX;
    this->m_lismBY     = PP->m_lismBY;
    this->m_lismBZ     = PP->m_lismBZ;
    this->m_lismUX     = PP->m_lismUX;
    this->m_lismUY     = PP->m_lismUY;
    this->m_lismUZ     = PP->m_lismUZ;
    this->m_Bref       = PP->m_Bref;
    
    this->m_sunN       = PP->m_sunN;
    this->m_sunV       = PP->m_sunV;
    this->m_sunT       = PP->m_sunT;
    this->m_sunB       = PP->m_sunB;
    this->m_sunTILT    = PP->m_sunTILT;    
    this->m_dirBrNorth = PP->m_dirBrNorth;    
    this->m_TMLIM      = PP->m_TMLIM;
    this->m_bSunGravity= PP->m_bSunGravity;
    this->m_bSunHeating= PP->m_bSunHeating;
    this->m_sunBmonopolar = PP->m_sunBmonopolar;

    this->m_initR          = PP->m_initR;
    this->m_sunIntBCRadius = PP->m_sunIntBCRadius;

    this->m_gamma      = PP->m_gamma;
    this->m_netN       = PP->m_netN;
    this->m_verbosity  = PP->m_verbosity;
    this->m_subproblem = PP->m_subproblem;

    this->m_sunZ2      = PP->m_sunZ2;
    this->m_sunSigmaC  = PP->m_sunSigmaC;
    this->m_sunLambda  = PP->m_sunLambda;
    this->m_lismZ2     = PP->m_lismZ2;
    this->m_lismSigmaC = PP->m_lismSigmaC;
    this->m_lismLambda = PP->m_lismLambda;

    this->m_sunN_PI    = PP->m_sunN_PI;
    this->m_sunT_PI    = PP->m_sunT_PI;
    this->m_lismN_PI   = PP->m_lismN_PI;
    this->m_lismT_PI   = PP->m_lismT_PI;
    
    this->m_writeH5     = PP->m_writeH5;
    this->m_writeH5File = PP->m_writeH5File;
    this->m_writeT0     = PP->m_writeT0;        
    this->m_writeDt     = PP->m_writeDt;     
    this->m_writeRotating = PP->m_writeRotating;   
    
    this->m_region_tracer = PP->m_region_tracer;
    
    this->m_output_region = PP->m_output_region;
    
    this->m_photoionize    = PP->m_photoionize;
    this->m_const_H        = PP->m_const_H;
    
    this->m_output_vecSPH  = PP->m_output_vecSPH;
    this->m_RegSafeZTop    = PP->m_RegSafeZTop;
    this->m_RegSafeZBot    = PP->m_RegSafeZBot;

    this->setFortranCommonSet();
  }

  this->m_adaptRmin  = PP->m_adaptRmin;
  this->m_adaptRmax  = PP->m_adaptRmax;

  this->m_helioAdapt = PP->m_helioAdapt;
  this->m_adaptXL1   = PP->m_adaptXL1;
  this->m_adaptXL2   = PP->m_adaptXL2;

  MultiFluidProblem::copy_PhysProblem(a_PP);
}

// Number additional variables for writing to plot file
int HeliosphericProblem::numPlotVars()
{
#ifdef OUTPUT_DT // output timesteps for each population
  return nFluids();
#endif

#ifdef OUTPUT_RELGRADB
  return 1;
#endif

  int numPlotVars = 1;                         // Plasma temperature
  if  (m_output_region == true) numPlotVars++; // Plasma temperature + region
  
  return numPlotVars; 
  
  if( nFluids() > 1 )
  {
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
    if (nFluids() == 2)
    {
      PickupIons* pPI = eqSys->getPickupIons();

      if( pPI != NULL )
      {
        return numPlotVars+6;  // Six source terms
      }
    }
    return numPlotVars+4;  // Four source terms
  }        
  return numPlotVars; 
}

// Names of the additional variables for writing to plot file  
Vector<string> HeliosphericProblem::plotNames()
{
  Vector<string> retval;  

#ifdef OUTPUT_DT // output timesteps for each population
  for (int i=0; i<numPlotVars(); i++)
  { 
    char buffer[40];
    sprintf(buffer,"dt%i", i);retval.push_back(buffer);
  }  
  return retval;    
#endif

#ifdef OUTPUT_RELGRADB
  retval.push_back(std::string("relGradB"));
  return retval;    
#endif



  if (numPlotVars() >= 1)
    retval.push_back("T");
  if (numPlotVars() >= 2)  
    retval.push_back("region");
  if (numPlotVars() >= 2+4)  
  {
    retval.push_back("srho");
    retval.push_back("x-smom");
    retval.push_back("y-smom");
    retval.push_back("seng");
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
    if (nFluids() == 2)
    {
      PickupIons* pPI = eqSys->getPickupIons();

      if( pPI != NULL )
      {
        retval.push_back("srho_pi");
        retval.push_back("seng_pi");
      }
    }        
  }    
        
  return retval;    
}
  
// Calculates variables for plotting using primitive variables
void HeliosphericProblem::calcPlotVars(FArrayBox& a_vars,
                           int              a_comp,
                           const FArrayBox& a_W,
                           const Box&       a_box)
{
  int nPlotVars = numPlotVars();
  if (nPlotVars == 0) return;  

#ifdef OUTPUT_DT // output timesteps for each population
  FArrayBox U(a_box, m_eqSys->numStates());
  m_eqSys->primToState(U,a_W,a_box);
  m_csh->transCartesianVectToCurv(U,a_box,m_level);
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  
  FArrayBox vars(a_vars.box(),nPlotVars);
  
  if ( (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian) ||
       (CoordinateSystem == CoordinateSystemHandler::CS_Cylindrical) ||
       (CoordinateSystem == CoordinateSystemHandler::CS_Axisymmetric))
  {
  
    vars.setVal(0.0,a_box,0,nPlotVars);

    Real dx = m_csh->dx(0,m_level);

    FORT_MAXWAVESPEED(CHF_FRA1(vars,0),
                      CHF_CONST_FRA(U),
                      CHF_BOX(a_box));

    for (int iFluid = 1; iFluid < nFluids(); ++iFluid)
    {
      int startRho = eqSys->densityIndexCons(iFluid);

      FORT_MAXWAVESPEED_E(CHF_FRA1(vars,iFluid),
                           CHF_CONST_FRA(U),
                           CHF_CONST_INT(startRho),
                           CHF_BOX(a_box));
    }


    vars.invert(dx);
  }
    
  if ( (CoordinateSystem == CoordinateSystemHandler::CS_Polar) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_PolarAxisym) || 
       (CoordinateSystem == CoordinateSystemHandler::CS_Spherical) )
  {        
    vars.setVal(numeric_limits<Real>::max(),a_box,0,nPlotVars);                
        
    FORT_MINDT_SPHERICAL(CHF_FRA1(vars,0),
                CHF_CONST_FRA(U),   
                CHF_CONST_INT(m_level),
                CHF_BOX(a_box));
                
    for (int iFluid = 1; iFluid < nFluids(); ++iFluid)
    { 
      int startRho = eqSys->densityIndexCons(iFluid);
       
      FORT_MINDT_SPHERICAL_E(CHF_FRA1(vars,iFluid),
                CHF_CONST_FRA(U),   
                CHF_CONST_INT(startRho),
                CHF_CONST_INT(m_level),
                CHF_BOX(a_box));
    } 
  }
  a_vars.copy(vars,0,a_comp,nPlotVars);
  return;
#endif   

#ifdef OUTPUT_RELGRADB
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  FArrayBox B(a_box,1);
  FArrayBox modGradB(a_box,1);
  B.setVal(0.0);
  modGradB.setVal(0.0);

  int iBGN = WBX;
  int iEND = WBZ;

  FORT_GETVECTMAGNITUDE(
            CHF_FRA1(B,0),
            CHF_CONST_FRA(a_vars),
            CHF_CONST_INT(iBGN),
            CHF_CONST_INT(iEND),
            CHF_BOX(a_box));
                
  if (CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
  {
    FORT_GETRELGRAD_SPHERICAL(
            CHF_FRA1(modGradB,0),
            CHF_CONST_FRA1(B,0),            
            CHF_BOX(a_box),
            CHF_CONST_INT(m_level));
    
  }
  a_vars.copy(modGradB,0,a_comp,1);
  
  #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(a_vars)
        );
    #endif
    
  return;
#endif 

  FArrayBox W(a_box, a_W.nComp());  
  W.copy(a_W);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
    
  BaseFab<int> Region;
  FArrayBox dummyFab;      
  defineRegions(W, dummyFab, Region, a_box);
      
  Real coeff = m_gamma*m_lismM*m_lismM*m_lismT;
      
  int iT      = a_comp+0;
  int iRegion = a_comp+1;
  
  BoxIterator bit(a_box);
  for (bit.begin(); bit.ok(); ++bit)
  { 
    IntVect iv = bit();
    Real T = coeff*W(iv,WPRES)/W(iv,WRHO);
    if (numPlotVars() >= 1)  a_vars.set(iv, iT,      T);    
    if (numPlotVars() >= 2)  a_vars.set(iv, iRegion, (Real)(Region(iv, 0)) );    
  }
  
  
  int isrho   = 2;
  if (numPlotVars() >= 2+4)  
  {
    int iFluids = nFluids();
    Real dt = 1.0;
    FArrayBox S( a_box, a_W.nComp() );
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    int iRhoN  = eqSys->densityIndexPrim(1);  
    
    if (iFluids == 2)
    {
      PickupIons* pPI = eqSys->getPickupIons();

      if( pPI != NULL )
      {
        int iRhoPIU  = pPI->consInterval().begin();
        int iRhoPIW  = pPI->primInterval().begin();

        FORT_CHARGE_EXCHANGE_2FPI_REG( CHF_FRA(S),
                                       CHF_CONST_FRA(a_W),
                                       CHF_CONST_FIA1(Region,0),
                                       CHF_CONST_REAL(dt),
                                       CHF_CONST_INT(iRhoN),
                                       CHF_CONST_INT(iRhoPIU),
                                       CHF_CONST_INT(iRhoPIW),
                                       CHF_BOX(a_box) );

        a_vars.copy( S, iRhoPIU, a_comp+isrho+4, 2 );
      }
      else
      {
        FORT_CHARGE_EXCHANGE_2F( CHF_FRA(S),
                                 CHF_CONST_FRA(a_W),
                                 CHF_CONST_FIA1(Region,0),
                                 CHF_CONST_REAL(dt),  
                                 CHF_CONST_INT(iRhoN),                           
                                 CHF_BOX(a_box) );
      }
    }
    if (iFluids == 4)
    {
      PickupIons* pPI = eqSys->getPickupIons();

      if( pPI != NULL)
      {
        int iRhoPIU = pPI->consInterval().begin();
        int iRhoPIW = pPI->primInterval().begin();

        FORT_CHARGE_EXCHANGE_4FPI( CHF_FRA(S),
                                   CHF_CONST_FRA(a_W),
                                   CHF_CONST_FIA1(Region,0),
                                   CHF_CONST_REAL(dt),
                                   CHF_CONST_INT(iRhoN),
                                   CHF_CONST_INT(iRhoPIU),
                                   CHF_CONST_INT(iRhoPIW),
                                   CHF_BOX(a_box) );

        a_vars.copy(S, iRhoPIU, a_comp+isrho+4, 2);
      }
      else
      {
        FORT_CHARGE_EXCHANGE_4F( CHF_FRA(S),
                                 CHF_CONST_FRA(a_W),
                                 CHF_CONST_FIA1(Region,0),
                                 CHF_CONST_REAL(dt),
                                 CHF_CONST_INT(iRhoN),
                                 CHF_BOX(a_box) );
      }
    }
            
    a_vars.copy(S,0,a_comp+isrho,  3);
    a_vars.copy(S,4,a_comp+isrho+3,1); // copy s-eng, s-momz is skipped
    
    #ifndef NDEBUG
      FORT_VIEWBOXDATA(
        CHF_FRA(a_vars)
        );
    #endif
  }
  
}

void HeliosphericProblem::primForPlot(       FArrayBox& a_W,
                                       const Box&       a_box)
{
  //return;
  double scalP = m_Bref*m_Bref;
  double scalV = 1.0e-5*m_lismV;
  
  //m_csh->transCartesianVectToCurv(a_W, a_box, m_level);

  a_W.mult( m_lismN, WRHO     );
  a_W.mult( scalV,   WVELX, 3 );  
  a_W.mult( scalP,   WPRES    );  
  a_W.mult( m_Bref,  WBX,   3 );
  
  if ((m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)&&(m_output_vecSPH == true))
  {
    Real dx = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? m_csh->dx(0,m_level) : 1);
    Real x,y,z = 0.0;
    
    Real vx,vy,vz,vr,vp,vt,bx,by,bz,br,bp,bt;
    Real r,cosT,sinT,cosF,sinF;
    
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    { 
      IntVect iv = bit();
      
      D_TERM(
        x = (iv[0]-m_sunIJK[0]+0.5)*dx;, 
        y = (iv[1]-m_sunIJK[1]+0.5)*dx;, 
        z = (iv[2]-m_sunIJK[2]+0.5)*dx;);
        
      r = sqrt(D_TERM(x*x,+y*y,+z*z));
      
      cosT = z/r;
      sinT = sqrt(1-cosT*cosT);
      cosF = x/sqrt(x*x+y*y);
      sinF = y/sqrt(x*x+y*y);
      
      vx   = a_W.get(iv,WVELX);
      vy   = a_W.get(iv,WVELY);
      vz   = a_W.get(iv,WVELZ);
      
      vr =  (vx*cosF + vy*sinF)*sinT + vz*cosT;
      vp =  -vx*sinF + vy*cosF;
      vt =  (vx*cosF + vy*sinF)*cosT - vz*sinT;
      
      a_W.set(iv,WVELX,vr);
      a_W.set(iv,WVELY,vp);
      a_W.set(iv,WVELZ,vt);
      
      bx   = a_W.get(iv,WBX);
      by   = a_W.get(iv,WBY);
      bz   = a_W.get(iv,WBZ);
      
      br =  (bx*cosF + by*sinF)*sinT + bz*cosT;
      bp =  -bx*sinF + by*cosF;
      bt =  (bx*cosF + by*sinF)*cosT - bz*sinT;
      
      a_W.set(iv,WBX,br);
      a_W.set(iv,WBY,bp);
      a_W.set(iv,WBZ,bt);

    }
  }
  
  int Fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  
  for (int iFluid = 1; iFluid < Fluids; ++iFluid)
  { 
    int iRho = eqSys->densityIndexPrim(iFluid);
    a_W.mult( m_lismN, iRho      );
    a_W.mult( scalV,   iRho+1, 3 );  
    a_W.mult( scalP,   iRho+4    );  
  }

  TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
  if( pTurbMod != NULL )
  {
    if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
    {
      int iZ2  = pTurbMod->primInterval().begin();
      double scalZ = scalV*scalV;
      a_W.mult( scalZ, iZ2 );
    }
  }

  PickupIons * pPickUp = m_eqSys->getPickupIons();
  if( pPickUp != NULL )
  {
    if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
    {
      int iRhoPI  = pPickUp->primInterval().begin();
      a_W.mult( m_lismN, iRhoPI     );
      a_W.mult( scalP,   iRhoPI + 1 );
    }
  }
}


/// Things to do before advancing one level by one time step.
// We perfrom path-through BC
void HeliosphericProblem::preTimeStep(LevelData<FArrayBox>&  a_U, Real a_time)
{ 
  double test_time1, test_time2;
  test_time1 = MPI_Wtime();
 
  int fluids = nFluids();
  
  if (m_verbosity >= 2)  pout() << "Enter HeliosphericProblem::preTimeStep"  << endl;
  
  if ((m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) && (fluids > 1))
  { 
    int sign      = 1;        
    int jsize     = m_domain.domainBox().size(1);
    
    DataIterator dit1,dit2;    
    for(dit1 = m_lstrips.dataIterator(); dit1.ok(); ++dit1) m_lstrips[dit1].setVal(-1.0,0);
    
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
    Interval intnH1(eqSys->densityIndexCons(1),eqSys->densityIndexCons(1) + UNUM_E - 1),intnH4;
    
    //int urho1 = eqSys->densityIndexCons(1);
    //int ueng1 = urho1+UNUM_E-1;
    
    a_U.copyTo(intnH1, m_lstrips, Interval(0,UNUM_E-1));
    if (fluids>=5)
    {
      intnH4.define(eqSys->densityIndexCons(4), eqSys->densityIndexCons(4) + UNUM_E - 1);      
      a_U.copyTo(intnH4, m_lstrips, Interval(UNUM_E,2*UNUM_E-1));
    }
          
    dit1 = m_lstrips.dataIterator();    
    dit2 = m_ladjstr.dataIterator();          
    for(; dit1.ok(); ++dit1,++dit2)
    {
      FArrayBox& from = m_lstrips[dit1];
      FArrayBox& to   = m_ladjstr[dit2];      
      
      //int vector[1] = {1};
      //m_csh->transCartesianVectToCurv(from,vector,1,from.box(),m_level);                    
                
      FORT_NEUTRALS_SUNBCSPHERICAL(
        CHF_CONST_FRA(from),
        CHF_FRA(to),       
        CHF_CONST_INT(jsize),       
        CHF_CONST_INT(sign) ); 
        
      sign++;
    }    
    m_ladjstr.copyTo(Interval(0,UNUM_E-1), a_U, intnH1);
    if (fluids>=5)
      m_ladjstr.copyTo(Interval(UNUM_E,2*UNUM_E-1), a_U, intnH4);                    
  }
  
  if (m_writeH5 == true)
  {
    
    if (a_time >=m_writeT0)
    {
      if (m_lastInd < 0) 
        writeSphericalSlice(a_U, a_time);
      else
      if (m_lastWrittenTime + m_writeDt < a_time)
        writeSphericalSlice(a_U, a_time);      
    } 
  }
  
  if (m_verbosity >= 2)  pout() << "Leave HeliosphericProblem::preTimeStep"  << endl;

/*  double time_diff, time_diff_max;
  test_time2 = MPI_Wtime();
  time_diff = test_time2 - test_time1;
  MPI_Reduce(&time_diff,&time_diff_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  int myrank; //int rank = 0, total_processes, myrank;
  //MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  stringstream stream_tal;
  //int curStep_tal = getCurStep();
  //stream_tal << "HeliosphericProblem::explicitSource time taken for rank " << myrank << " at iteration " << curStep_tal << " : " << test_time2-test_time1 <<endl;
  if(myrank==0){
    stream_tal << "HeliosphericProblem::preTimeStep max time taken at iteration ?: " << time_diff_max <<endl;
    //cout << stream_tal.str();
  }; // else {
     //stream_tal << "HeliosphericProblem::preTimeStep time taken for rank " << myrank << " at iteration ?: " << time_diff <<endl;
    //cout << stream_tal.str();
  //};*/
  
}

/// Calculate problem specific source terms
/**
 */
void HeliosphericProblem::explicitSource(       FArrayBox & a_U,
                                     FArrayBox & a_S,
                               const FArrayBox & a_W,
                                  BaseFab<int> & a_REG,
                               const Real      & a_dt,                               
                               const Box       & a_box)
{
  double test_time1, test_time2;
  test_time1 = MPI_Wtime();
  
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
            FORT_CHARGE_EXCHANGE_2F( CHF_FRA(a_S),
                                     CHF_CONST_FRA(W),
                                     CHF_CONST_FIA1(a_REG,0),
                                     CHF_CONST_REAL(a_dt),
                                     CHF_CONST_INT(iRhoN),
                                     CHF_BOX(a_box) );
          }
        }
        break;
      case 3: FORT_CHARGE_EXCHANGE_3F( CHF_FRA(a_S),
                                       CHF_CONST_FRA(W),
                                       CHF_CONST_FIA1(a_REG,0),
                                       CHF_CONST_REAL(a_dt),
                                       CHF_CONST_INT(iRhoN),
                                       CHF_BOX(a_box) );
              break;
      case 4: 
        {
          PickupIons* pPI = eqSys->getPickupIons();

          if( pPI != NULL)
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

    // Add charge exchange source terms for plasma
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
//    // Add charge exchange source terms for neutrals
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
    Real Wnet[WNUM_E] = {m_netN/m_lismN, m_lismUX, m_lismUY, m_lismUZ, eos_k*m_netN*m_lismT/(m_lismN*eos_mp*m_lismV*m_lismV)};
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

  
  
  
  if (m_bSunGravity == true)
  {
    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
    {            
      FORT_SUNGRAVITYSPH( CHF_FRA(a_S),
                       CHF_CONST_FRA(a_W),
                       CHF_CONST_REAL(a_dt),                     
                       CHF_CONST_INT(m_level),
                       CHF_BOX(a_box) );
                       
      iBGN   = UMOMX;
      iEND   = UENG;    
          
      FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );
    }
  }
  
  if (m_bSunHeating == true)
  {
    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
    {            
      FORT_SUNHEATINGSPH( CHF_FRA(a_S),
                       CHF_CONST_FRA(a_W),
                       CHF_CONST_REAL(a_dt),                     
                       CHF_CONST_INT(m_level),
                       CHF_BOX(a_box) );
                       
      iBGN   = UENG;
      iEND   = UENG;    
          
      FORT_ADDSOURCES( CHF_FRA(a_U),
                     CHF_CONST_FRA(a_S),
                     CHF_CONST_INT(iBGN),
                     CHF_CONST_INT(iEND),
                     CHF_BOX(a_box) );
    }
  }

/*  double time_diff, time_diff_max;
  test_time2 = MPI_Wtime();
  time_diff = test_time2 - test_time1;
  MPI_Reduce(&time_diff,&time_diff_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  int myrank;   //int rank = 0, total_processes, myrank;
  //MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  stringstream stream_tal; */ //uncomment this and above for data
//  int curStep_tal = getCurStep();
/*  if(myrank==0){
    stream_tal << "HeliosphericProblem::explicitSource max time taken at iteration " << curStep_tal << " : " << time_diff_max <<endl;
    cout << stream_tal.str();
  } else {
    stream_tal << "HeliosphericProblem::explicitSource time taken for rank " << myrank << " at iteration " << curStep_tal << " : " << time_diff <<endl;
    cout << stream_tal.str();
  };*/

/*  if(myrank==0){
    stream_tal << "HeliosphericProblem::explicitSource max time taken at iteration ?: " << time_diff_max <<endl;
    cout << stream_tal.str();
  } else {
    stream_tal << "HeliosphericProblem::explicitSource time taken for rank " << myrank << " at iteration ?: " << time_diff <<endl;
    cout << stream_tal.str();
  }; */

}                               
                                                             // Fill ghost cells
void HeliosphericProblem::fillGhostCells(       FArrayBox&      a_W,
                                          const FArrayBox&      a_U,
                                          const int&            a_dir,
                                          const Real&           a_time)
{
  double test_time1, test_time2;
  test_time1 = MPI_Wtime();

  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();

                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;

  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  
  int fluids = nFluids();  
  

  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    int indW = WBox.bigEnd( a_dir );
    int indD =  tmp.bigEnd( a_dir );  
  
    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                           // Fill the ghost cells
      FORT_HELIOGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_CONST_INT(m_ls_indices.m_iRegTr),
                    CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2  = pTurbMod->primInterval().begin();
            FORT_HELIOGS_TM( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOGS_PI( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }
    }

    indW     = WBox.smallEnd( a_dir );
    indD     =  tmp.smallEnd( a_dir );

    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indW, indD - indW );

                                                           // Fill the ghost cells
      FORT_HELIOGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_CONST_INT(m_ls_indices.m_iRegTr),
                    CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2  = pTurbMod->primInterval().begin();
            FORT_HELIOGS_TM( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOGS_PI( CHF_FRA(a_W),
                             CHF_CONST_INT(sign),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }
    }
  } 
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {
    CSHSpherical* csh_sph = static_cast<CSHSpherical*>(m_csh);
    
    int indW,indD,nGS,jsize;    
    
    jsize = m_domain.domainBox().size(1);
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain;  
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = abs(indW-indD);
            
        
    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
      
      if( a_dir == 0 )
      {
        FORT_HELIOGSSPHERICAL(
            CHF_FRA(a_W),
            CHF_CONST_FRA(a_U),
            CHF_CONST_INT(sign),
            CHF_CONST_INT(a_dir),
            CHF_CONST_INT(jsize),
            CHF_CONST_INT(iRhoN),
            CHF_CONST_INT(fluids),
            CHF_CONST_INT(m_ls_indices.m_iHCS),
            CHF_CONST_INT(m_ls_indices.m_iHCSb),
            CHF_CONST_INT(m_ls_indices.m_iRegTr),
            CHF_CONST_INT(m_level),
            CHF_CONST_REAL(a_time),
            CHF_BOX(boundaryBox));

        if( m_eqSys != NULL )
        {
          TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
          if( pTurbMod != NULL )
          {
            if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
            {
              int iZ2  = pTurbMod->primInterval().begin();
              FORT_HELIOGSSPHERICAL_TM( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_HELIOGSSPHERICAL_PI( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }
        }
      }
            
      if (a_dir == 2)
      {
        m_eqSys->stateToPrim(a_W, a_U, boundaryBox);
        csh_sph->zaxisBC(a_W, boundaryBox, sign);

        //FORT_ZAXISFIRSTORDER(
        // CHF_FRA(a_W),            
        // CHF_CONST_INT(sign),
        // CHF_CONST_INT(indD),
        // CHF_CONST_INT(nGS),
        // CHF_BOX(boundaryBox));
      }
    }

    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
  
      if (a_dir == 0)
      {
        FORT_HELIOGSSPHERICAL(
          CHF_FRA(a_W),           
          CHF_CONST_FRA(a_U),           
          CHF_CONST_INT(sign),
          CHF_CONST_INT(a_dir),
          CHF_CONST_INT(jsize),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(m_ls_indices.m_iHCS),
          CHF_CONST_INT(m_ls_indices.m_iHCSb),
          CHF_CONST_INT(m_ls_indices.m_iRegTr),
          CHF_CONST_INT(m_level),
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox)); 
          
        // Plasma parameters at the inner boundary
//why is this here?
/*        
        FORT_HELIOGSPLASMASPHERICAL(
          CHF_FRA(a_W),           
          CHF_CONST_FRA(a_U), 
          CHF_CONST_INT(m_subproblem),           
          CHF_CONST_INT(m_ls_indices.m_iHCS),
          CHF_CONST_INT(m_ls_indices.m_iHCSb),
          CHF_CONST_INT(m_ls_indices.m_iRegTr),
          CHF_CONST_INT(m_level),
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox)); 
*/
        if( m_eqSys != NULL )
        {
          TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
          if( pTurbMod != NULL )
          {
            if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
            {
              int iZ2  = pTurbMod->primInterval().begin();
              FORT_HELIOGSSPHERICAL_TM( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_HELIOGSSPHERICAL_PI( CHF_FRA(a_W),
                                        CHF_CONST_FRA(a_U),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_CONST_REAL(a_time),
                                        CHF_BOX(boundaryBox) );
            }
          }
        }
      }
      if (a_dir == 2)
      {
        m_eqSys->stateToPrim(a_W, a_U, boundaryBox);
        csh_sph->zaxisBC(a_W, boundaryBox, sign);
                                                                   
        //FORT_ZAXISFIRSTORDER(
        // CHF_FRA(a_W),            
        // CHF_CONST_INT(sign),
        // CHF_CONST_INT(indD),
        // CHF_CONST_INT(nGS),
        // CHF_BOX(boundaryBox));
      }            
    }    
        
  }

/*  double time_diff, time_diff_max;
  test_time2 = MPI_Wtime();
  time_diff = test_time2 - test_time1;
  MPI_Reduce(&time_diff,&time_diff_max,1,MPI_DOUBLE,MPI_MAX,0,MPI_COMM_WORLD);
  int myrank; //int rank = 0, total_processes, myrank;
  //MPI_Comm_size(MPI_COMM_WORLD, &total_processes);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  stringstream stream_tal;
  //int curStep_tal = getCurStep();
  //stream_tal << "HeliosphericProblem::explicitSource time taken for rank " << myrank << " at iteration " << curStep_tal << " : " << test_time2-test_time1 <<endl;
  if(myrank==0){
    stream_tal << "HeliosphericProblem::fillGhostCells max time taken at iteration ?: " << time_diff_max <<endl;
    //cout << stream_tal.str();
  };*/

}

                                                          // Set boundary fluxes
void HeliosphericProblem::fluxBC(       FArrayBox&      a_F,
                                        FArrayBox&      a_Bn,
                                  const FArrayBox&      a_WMinus,
                                  const FArrayBox&      a_WPlus,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side,
                                  const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  int sign;
  Box FBox = a_F.box();
  Box tmp = FBox;
  
  Real dx = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? m_csh->dx(0,m_level) : 1);
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  
    

                            // Determine which side and thus shifting directions
  if( a_side == Side::Lo )
  {
    sign = -1;
  }
  else
  {
    sign = 1;
  }

  tmp.shiftHalf(a_dir,sign);

                                 // Is there a domain boundary next to this grid
  if( !m_domain.contains( tmp ) )
  {
    tmp &= m_domain;

    Box boundaryBox;

                          // Find the strip of cells next to the domain boundary
    if( a_side == Side::Lo )
    {
      boundaryBox = bdryLo(tmp,a_dir);
    }
    else
    {
      boundaryBox = bdryHi(tmp,a_dir);
    }
    
    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
    {
      if (a_dir == 0)
      {
        m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
        
        for (int iFluid = 1; iFluid < fluids;++iFluid)        
          m_RSGD->fluxes( a_F, a_WPlus, a_WMinus, a_dir, iRhoN+(iFluid-1)*UNUM_E, boundaryBox );
        
      }
      if (a_dir == 2)
        a_F.setVal(0.0, boundaryBox, 0, a_F.nComp()); 
    } else
    {
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
        
      for (int iFluid = 1; iFluid < fluids;++iFluid)        
        m_RSGD->fluxes( a_F, a_WPlus, a_WMinus, a_dir, iRhoN+(iFluid-1)*UNUM_E, boundaryBox );
              
    }

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
    FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
    FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
    shiftWLeft .shiftHalf(a_dir, 1);
    shiftWRight.shiftHalf(a_dir,-1);        

    
                                                      // Set the boundary fluxes
    if( m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian )
    {
      FORT_HELIOBC( CHF_FRA(a_F),
                    CHF_FRA1(a_Bn,0),
                    CHF_CONST_FRA(shiftWLeft),
                    CHF_CONST_FRA(shiftWRight),
                    CHF_CONST_INT(sign),
                    CHF_CONST_REAL(dx),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_BOX(boundaryBox) );

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2    = pTurbMod->primInterval().begin();
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOBC_TM( CHF_FRA(a_F),
                             CHF_CONST_FRA(shiftWLeft),
                             CHF_CONST_FRA(shiftWRight),
                             CHF_CONST_INT(sign),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iZ2),
                             CHF_CONST_INT(iRhoZ2),
                             CHF_BOX(boundaryBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            FORT_HELIOBC_PI( CHF_FRA(a_F),
                             CHF_CONST_FRA(shiftWLeft),
                             CHF_CONST_FRA(shiftWRight),
                             CHF_CONST_INT(sign),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_INT(a_dir),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(boundaryBox) );
          }
        }
      }
    }

    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)     
    {
      if( a_dir == 0 )
      {
        FORT_HELIOBCSPHERICAL( CHF_FRA(a_F),
                    CHF_FRA1(a_Bn,0),
                    CHF_CONST_FRA(shiftWLeft),
                    CHF_CONST_FRA(shiftWRight),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(iRhoN),
                    CHF_CONST_INT(fluids),
                    CHF_BOX(boundaryBox) );

        if( m_eqSys != NULL )
        {
          TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
          if( pTurbMod != NULL )
          {
            if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
            {
              int iZ2    = pTurbMod->primInterval().begin();
              int iRhoZ2 = pTurbMod->consInterval().begin();
              FORT_HELIOBCSPHERICAL_TM( CHF_FRA(a_F),
                                        CHF_CONST_FRA(shiftWLeft),
                                        CHF_CONST_FRA(shiftWRight),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iZ2),
                                        CHF_CONST_INT(iRhoZ2),
                                        CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_HELIOBCSPHERICAL_PI( CHF_FRA(a_F),
                                        CHF_CONST_FRA(shiftWLeft),
                                        CHF_CONST_FRA(shiftWRight),
                                        CHF_CONST_INT(sign),
                                        CHF_CONST_INT(a_dir),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_BOX(boundaryBox) );
            }
          }
        }
      }
    }

    shiftWLeft .shiftHalf(a_dir,-1);
    shiftWRight.shiftHalf(a_dir, 1);
  }
}

                                                    // Set up initial conditions
void HeliosphericProblem::initialize( LevelData<FArrayBox>& a_U )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian ? m_csh->dx(0,m_level) : 1.0);
  
  int ibox = 0;
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN = eqSys->densityIndexCons(1);  

                                          // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
                                                     // Storage for current grid
    FArrayBox& U = a_U[dit()];

                                                          // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    if (m_verbosity >= 4)    
    {
      pout() << "   Box: " << ibox << " "; uBox.p(); pout().flush();
    }
        
                      // Set up initial condition in this grid
    if( m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian )
    {
      FORT_HELIOINIT( CHF_FRA(U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_BOX(uBox) );
                      
      if (m_subproblem == HPBC_SOLARCYCLE)
      {
        Real t = 0.0;
        FORT_HELIOREINIT_CYCLE( CHF_FRA(U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(t),
                      CHF_CONST_REAL(m_initR),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(uBox) );
      }

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOINIT_TM( CHF_FRA(U),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(iRhoZ2),
                               CHF_BOX(uBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->consInterval().begin();
            FORT_HELIOINIT_PI( CHF_FRA(U),
                               CHF_CONST_REAL(dx),
                               CHF_CONST_INT(iRhoPI),
                               CHF_BOX(uBox) );
          }
        }
      }
    }
    if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)  
    {
      FORT_HELIOINITSPHERICAL( CHF_FRA(U),                      
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iHCSb),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox));

//why is this here?
/*      FORT_HELIOINITPLASMASPHERICAL( CHF_FRA(U),
                      CHF_CONST_INT(m_subproblem),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(m_ls_indices.m_iHCSb),
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );
*/

      if (m_subproblem == HPBC_SOLARCYCLE)
      {
        Real t = 0.0;
        FORT_HELIOREINIT_CYCLE_SPHERICAL( CHF_FRA(U),
                      CHF_CONST_REAL(t),
                      CHF_CONST_REAL(m_initR),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );
      }


      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iRhoZ2 = pTurbMod->consInterval().begin();
            FORT_HELIOINITSPHERICAL_TM( CHF_FRA(U),
                                        CHF_CONST_INT(iRhoZ2),
                                        CHF_CONST_INT(m_level),
                                        CHF_BOX(uBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->consInterval().begin();
            FORT_HELIOINITSPHERICAL_PI( CHF_FRA(U),
                                        CHF_CONST_INT(iRhoPI),
                                        CHF_CONST_INT(m_level),
                                        CHF_BOX(uBox) );
          }
        }
      }
    }

    ibox++;
  }
}

void HeliosphericProblem::initialize(LevelData<FArrayBox>& a_U,
                                                Interval & a_comp)
{
  DataIterator dit = a_U.boxLayout().dataIterator();
    
  Real pref       = m_lismN*eos_mp*m_lismV*m_lismV;
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    
  // Population 1 is added
  if (a_comp.contains(eqSys->densityIndexPrim(1)))
  {
    int iRho = eqSys->densityIndexPrim(1);
    Real p   = (eos_k*m_netN*m_lismT)/pref;
    Real eng = p/(m_gamma-1.0) + 0.5*m_netN*m_lismV*m_lismV/(m_lismN*m_lismV*m_lismV);
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& U = a_U[dit()];
      U.setVal(m_netN/m_lismN,iRho);
      U.setVal(m_netN*m_lismUX/(m_lismN*m_lismV),iRho+UMOMX);
      U.setVal(m_netN*m_lismUY/(m_lismN*m_lismV),iRho+UMOMY);
      U.setVal(m_netN*m_lismUZ/(m_lismN*m_lismV),iRho+UMOMZ);
      U.setVal(eng,iRho+UENG);      
    }
  }
  
  int iFluid, Fluids = nFluids();  
  
  for (iFluid = 2; iFluid < Fluids; ++iFluid)
  if (a_comp.contains(eqSys->densityIndexPrim(iFluid)))
  { 
    int iRho = eqSys->densityIndexPrim(iFluid);
    for (dit.begin(); dit.ok(); ++dit)
    {
      FArrayBox& U = a_U[dit()];
    
      U.setVal(1.0e-7,iRho);
      U.setVal(0.0,iRho+UMOMX);
      U.setVal(0.0,iRho+UMOMY);
      U.setVal(0.0,iRho+UMOMZ);
      U.setVal(1.0e-7/(m_gamma-1.0),iRho+UENG);              
    }
  }
        
      

                                          // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    FArrayBox& U = a_U[dit()];
    if (a_comp.contains(m_ls_indices.m_iHCS))
      U.setVal(1.0,m_ls_indices.m_iHCS);
    if (a_comp.contains(m_ls_indices.m_iHCSb))
      U.setVal(1.0,m_ls_indices.m_iHCSb);
  }
}


                                              // Problem specific postprocessing
void HeliosphericProblem::postprocessing(       FArrayBox & a_U,
                                          const FArrayBox & a_W,
                                          const Real      & a_dt,
                                          const Real      & a_time,
                                          const Box       & a_box       )
{
  Real newTime = a_dt + a_time;  
  
/*  if (m_ls_indices.m_iRegTr > 0)
  {
    FORT_REGTRACER_REINIT( CHF_FRA(a_U),                   
                      CHF_CONST_INT(m_ls_indices.m_iRegTr),
                      CHF_BOX(a_box) );
  }*/
    
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {  
    Real dx = m_csh->dx(0,m_level);
    EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
    int iRhoN = eqSys->densityIndexCons(1);  
    int fluids = nFluids();
    
    if (m_subproblem == HPBC_SOLARCYCLE)
      FORT_HELIOREINIT_CYCLE( CHF_FRA(a_U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(newTime),
                      CHF_CONST_REAL(m_R0),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(a_box) );
    else        
      FORT_HELIOREINIT_DEFAULT( CHF_FRA(a_U),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_REAL(newTime),
                      CHF_CONST_REAL(m_R0),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(a_box) );

 
    PickupIons* pPI = eqSys->getPickupIons();//FF with the cycle, assumes spherically symm PI 
    if( pPI != NULL )
    {
     if( pPI->modelID() == PickupIons::PI_TWO_EQNS )
       {
       int iRhoPI  = pPI->consInterval().begin();
       FORT_HELIOREINIT_PI( CHF_FRA(a_U),
                             CHF_CONST_REAL(dx),
                             CHF_CONST_REAL(newTime),
                             CHF_CONST_INT(iRhoPI),
                             CHF_BOX(a_box) );
       }
    }

  }
  else if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)     
  {
    if (m_subproblem == HPBC_SOLARCYCLE)
    {
      EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
      int iRhoN = eqSys->densityIndexCons(1);  
      int fluids = nFluids();

      FORT_HELIOREINIT_CYCLE_SPHERICAL( CHF_FRA(a_U),
                      CHF_CONST_REAL(newTime),
                      CHF_CONST_REAL(m_R0),
                      CHF_CONST_INT(m_ls_indices.m_iHCS),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(a_box) );
    }
  }
}

// Pogorelov's (passthrough) boundary conditions for neutrals. Works for Cartesian coordinate system only 
void HeliosphericProblem::postTimeStep(LevelData<FArrayBox>&  a_U)
{
  int fluids = nFluids();
  
  // Pogorelov's (passthrough) boundary conditions for neutrals  
  if ((fluids > 1) && (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian))
  { 
    if( m_verbosity >= 3 )
    {
      pout() << "HeliosphericProblem::postTimeStep, pass through bc" << endl;
    }
  
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
//    for(dit1.end(),--dit1,dit2.begin(); dit1.ok(); --dit1,++dit2)
    for(; dit1.ok(); ++dit1,++dit2)
    {
      FArrayBox& from = m_lstrips[dit1];
      FArrayBox& to   = m_ladjstr[dit2];      
      
      //FArrayBox& to     = m_lstrips[dit1];
      //FArrayBox& from   = m_ladjstr[dit2];                  
                                           
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

                                                               // Define Regions
void HeliosphericProblem::defineRegions( const FArrayBox    & a_W,
                                               FArrayBox    & a_S,
                                               BaseFab<int> & a_R,
                                         const Box          & a_box)
{
  a_R.resize( a_box, 1 );
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();  
  Real dx = (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian ? m_csh->dx(0,m_level) : -1.0);
  
  if  (m_region_tracer)
  {    
/*    FORT_DEFINE_REGIONS_TRACER( CHF_CONST_FRA(a_W),
                          CHF_FIA1(a_R,0),
                          CHF_CONST_INT(m_ls_indices.m_iRegTr),
                          CHF_CONST_REAL(dx),
                          CHF_BOX(a_box) );*/

//FF 2022 this one does not require REGTRACER_REINIT-surface0 is not overridden and remains in [-1,1].
// Should work with lsreinit=1 in the inputs file. Works with Cartesian or spherical grids.
    FORT_DEFINE_REGIONS_TRACER_4F( CHF_CONST_FRA(a_W),
                          CHF_FIA1(a_R,0),
                          CHF_CONST_INT(m_ls_indices.m_iRegTr),
                          CHF_CONST_REAL(dx),
                          CHF_CONST_INT(m_level),
                          CHF_BOX(a_box) );


  } else
  {    
    FORT_DEFINE_REGIONS_2F( CHF_CONST_FRA(a_W),
                            CHF_FIA1(a_R,0),
                            CHF_CONST_REAL(dx),
                            CHF_BOX(a_box) );
  }                          
    
}

//                              Creates tagged cells for dynamic mesh refinement
void HeliosphericProblem::tagCells( const FArrayBox&  a_U,
                                    const Box&        a_box,
                                          IntVectSet& a_tags )
{

  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    switch( m_helioAdapt ) {
    case HPAD_LONGTAIL :
      tagCellsLongTail( a_U, a_box, a_tags );
      break;
    case HPAD_HPINSTAB :
      tagCellsHPInstab( a_U, a_box, a_tags );
      break;
    case HPAD_TAILDETAIL :
      tagCellsTailDetail( a_U, a_box, a_tags );
      break;
    default :
//      tagCellsCartesian( a_U, a_box, a_tags );
      tagCellsTailHP_Federico( a_U, a_box, a_tags );
    }
  }

  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {
    // SW-LISM with tilt
    //tagCellsTilt(a_U,a_box,a_tags);   
    
    //if (m_eqSys->numTrackingSurfaces()>0) tagCellsTiltLS(a_U,a_box,a_tags); // SW-LISM with level set
    
   //  tagCellsTiltLSV1V2(a_U,a_box,a_tags);
     tagCellsSpherical(a_U,a_box,a_tags);
  }
}

void HeliosphericProblem::tagCellsCartesian( const FArrayBox&  a_U,
                                             const Box&        a_box,
                                                   IntVectSet& a_tags )
{
  Real dx = m_csh->dx(0,m_level);  

  const Box& b = a_box;

  BoxIterator bit( b );

  /*if( m_adaptRmin > 0.0 )
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();

    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );

    Real dist2 = D_TERM( x*x, + y*y, + z*z );

    if( dist2 < m_adaptRmin*m_adaptRmin )
    {
      a_tags |= iv;
    }
  }*/
  
  //if (nFluids()>1) // For MHD case R0 is too large to regrid it.
  
  // Tag all m_R0Box
  if (!m_R0Box.isEmpty())  
  {
    Box R0Box(m_R0Box);    
    R0Box.grow(1);
        
    BoxIterator bit(R0Box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      a_tags |= iv;      
    }
  } 
    
  Real Rmin  = m_R0 - 2.0*dx;
  Real Rmin2 = Rmin*Rmin;
  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    
    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );

    Real dist2 = D_TERM( x*x, + y*y, + z*z );

    
    int reg = Region(iv,0);
    
    // Add extra patches to the heliosheath
    if ((reg == 2) && (dx > 3.5) && (x > -50.0)) a_tags |= iv;      
    
    if ((reg == 3) && (dx > 1.6) && (dist2 > Rmin2)) a_tags |= iv;      
    
    if (m_subproblem == HPBC_KINETIC) 
    {    
      //if (reg == 3) a_tags |= iv;    
      if ((reg == 2) && (dx > 7.0) && (iv[0]*dx > 200.0)) a_tags |= iv;
    }
  }
  
  if (m_subproblem == HPBC_KINETIC)
  {    
    a_tags |= m_kineticBox;
  }
}

// Level 1 covers all heliosheath.
// Level 2 heliosheath within 500 AU
// Level 3 near the Sun
void HeliosphericProblem::tagCellsTailDetail(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags)
{
  Real dx = m_csh->dx(0,m_level);  

  const Box& b = a_box;

  BoxIterator bit( b );
  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);
  
  // Tag box near the Sun
  if (!m_R0Box.isEmpty())  
  {
    Box R0Box(m_R0Box);    
    R0Box.grow(1);
        
    BoxIterator bit(R0Box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                  
      a_tags |= iv;      
    }
  } 
  
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();

    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );
    
    Real dist2 = D_TERM( x*x, + y*y, + z*z );
	
	int reg = Region(iv,0);
	
	if( (reg == 2) || (reg == 3) ) {
      if( (m_level == 0) && (x > m_adaptXL1) ) {
        a_tags |= iv;
	  } else {
        if( (m_level == 1) && (x > m_adaptXL2) ) {
          a_tags |= iv;
		}
	  }
	}
  }       
}


void HeliosphericProblem::tagCellsLongTail(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags)
{
  Real dx = m_csh->dx(0,m_level);  

  Real R = 200;
  if( m_level == 0 ) R   = 500; else
  if( m_level == 1 ) R   = 150; else
  if( m_level == 2 ) R   = 50;

  const Box& b = a_box;

  BoxIterator bit( b );
  
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();

    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );
    
    Real dist2 = D_TERM( x*x, + y*y, + z*z );

    if( dist2 < R*R )
    {
      a_tags |= iv;
    }
  }       
}

void HeliosphericProblem::tagCellsTailHP_Federico(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags)
{
#if CH_SPACEDIM == 3  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);
  
  Real dx = m_csh->dx(0,m_level);  


  // FF Tag box near the Sun
  if (!m_R0Box.isEmpty())
  {
    Box R0Box(m_R0Box);
    R0Box.grow(4);
        
    BoxIterator bit(R0Box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();
      a_tags |= iv;
    }
  }



  const Box& b = a_box;

  BoxIterator bit( b );
  IntVectSet hp_tags;
  
  Real zmax = dx*m_domain.domainBox().size(2);
  
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();

    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );
    
    if (z > (zmax - 100.0)) continue;
    
    Real dist2 = D_TERM( x*x, + y*y, + z*z );
    
    int reg = Region(iv,0);
    
    Real R = 4000;  //  160;    
    if (m_level == 0) 
    {      
      if (reg == 3) a_tags |= iv;
      else 
      {
        if ((reg == 2) && (dist2 < R*R)) a_tags |= iv;
      }
    }

    // nose refinement 
    if (m_level == 0)
    {
        if (sqrt((x-200)*(x-200)+y*y+z*z)  < 200) a_tags |= iv;
    }
 
    
    if (m_level == 1) 
    {
      if ( ((reg == 2) || (reg == 3)) && (dist2 < 200*200)) a_tags |= iv;
    }
        
    
    bool frontHP = true;
    Real fromtHPdist = 500.0;        
    if (m_level >=2) fromtHPdist = 200.0;
        
            
    Real stopHP_refinement = 2000; 
    if (m_ls_indices.m_iRegTr > 0)
    {
      Real reg_tr = a_U.get(iv,m_ls_indices.m_iRegTr);
      if (fabs(reg_tr)<0.7)
      {
        if ((m_level == 0) && (dist2 < stopHP_refinement*stopHP_refinement)) hp_tags |= iv;
        if ((m_level == 0) && (frontHP == true))
        {
          if ((dist2 < fromtHPdist*fromtHPdist ) && (x > 0.0))
          {
            hp_tags |= iv;
          }          
        } else
        {
          if (m_level <= 1 ) hp_tags |= iv;
          if ((m_level == 2) && (dist2 < stopHP_refinement*stopHP_refinement)) hp_tags |= iv;
        }
      }
    }
    
  }
  
  hp_tags.grow(2);
  a_tags |= hp_tags;
#endif       
}

void HeliosphericProblem::tagCellsHPInstab(const FArrayBox&  a_U,
                         const Box&        a_box,
                               IntVectSet& a_tags)
{
#if CH_SPACEDIM == 3  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);
  
  Real dx = m_csh->dx(0,m_level);  

  const Box& b = a_box;

  BoxIterator bit( b );
  IntVectSet hp_tags;
  
  Real zmax = dx*m_domain.domainBox().size(2);
  
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();

    Real D_DECL(x,y,z);

    D_EXPR(
      x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
      y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
      z = (iv[2] + 0.5)*dx - m_sunXYZ[2]
    );
    
    if (z > (zmax - 100.0)) continue;
    
    Real dist2 = D_TERM( x*x, + y*y, + z*z );
    
    int reg = Region(iv,0);
    
    Real R = 160;    
    if (m_level == 0) 
    {      
      if (reg == 3) a_tags |= iv;
      else 
      {
        if ((reg == 2) && (dist2 < R*R)) a_tags |= iv;
      }
    }
    
    if (m_level == 1) 
    {
      if ( ((reg == 2) || (reg == 3)) && (x > 0.0)) a_tags |= iv;
    }
        
    
    bool frontHP = true;
    Real fromtHPdist = 500.0;        
    if (m_level >=2) fromtHPdist = 300.0;
        
            
    Real stopHP_refinement = 800.0;
    if (m_ls_indices.m_iRegTr > 0)
    {
      Real reg_tr = a_U.get(iv,m_ls_indices.m_iRegTr);
      if (fabs(reg_tr)<0.5)
      {
        if ((m_level == 0) && (dist2 < stopHP_refinement*stopHP_refinement)) hp_tags |= iv;
        if (frontHP == true)
        {
          if ((dist2 < fromtHPdist*fromtHPdist ) && (x > 0.0))
          {
            hp_tags |= iv;
          }          
        } else
        {
          if (dist2 < stopHP_refinement*stopHP_refinement) hp_tags |= iv;
        }
      }
    }
    
  }
  
  hp_tags.grow(2);
  a_tags |= hp_tags;
#endif       
}


// Tags HCS region (no level set)                              
void HeliosphericProblem::tagCellsTilt(const FArrayBox&  a_U,
                    const Box&        a_box,
                          IntVectSet& a_tags)
{
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);

  BoxIterator bit( a_box ); RealVect cv_coord,cs_coord;

  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);
    
    int reg = Region(iv,0);
    Real theta1 = (90.0-(m_sunTILT + 5.0))*d_PI_180;
    Real theta2 = (90.0+(m_sunTILT + 5.0))*d_PI_180;
    
    Real th_axis = 8.0;
    Real theta3 = (th_axis)*d_PI_180;
    Real theta4 = (180.0-th_axis)*d_PI_180;
    
    //Real phi    = 6.0;
    Real phi    = 9.0;
    Real phi1   = (90.0  - phi)*d_PI_180;
    Real phi2   = (270.0 + phi)*d_PI_180;
    
    if (m_level == 0)
    {        
              
      if ((reg == 3) && 
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      if ((reg == 2) && 
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
              
    }    
    if (m_level == 1)
    {  
      theta1 = (90.0-(m_sunTILT + 3.0))*d_PI_180;
      theta2 = (90.0+(m_sunTILT + 3.0))*d_PI_180;
    
      phi    = 15.0;
      phi1   = (90.0  - phi)*d_PI_180;
      phi2   = (270.0 + phi)*d_PI_180;
    
      if ((reg == 3) && (cv_coord[0]>30.0) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
    
      theta3 = (30.0)*d_PI_180;
      theta4 = (90.0+(m_sunTILT + 10.0))*d_PI_180;
      if ((reg == 2) && 
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }
    if (m_level == 2)
    {
      Real phi    = 30.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;
      
      if ((reg == 3) && (cv_coord[0]>70.0) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
    
      theta3 = (30.0)*d_PI_180;
      theta4 = (90.0+(m_sunTILT + 10.0))*d_PI_180;
      if ((reg == 2) && 
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }

    if (m_level == 3)
    {
      Real phi    = 45.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-(m_sunTILT + 0.5))*d_PI_180; 
      if (cv_coord[0]>100.0) theta3 = (90.0-(m_sunTILT + 2.5))*d_PI_180; 
      if (cv_coord[0]>110.0) theta3 = (90.0-(m_sunTILT + 5.0))*d_PI_180; 
      theta4 = (90.0+(m_sunTILT + 0.5))*d_PI_180;
      if ((reg == 2) && //(cv_coord[0]>90.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }

    if (m_level == 4)
    {
      Real phi    = 60.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-20.0)*d_PI_180;         
      theta4 = (90.0+20.0)*d_PI_180;
      if ((reg == 2) && (cv_coord[0]>108.0) && (cv_coord[0]<117.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    } 
    
  }
}
 
// Tags HCS region using level set data
void HeliosphericProblem::tagCellsTiltLS(const FArrayBox&  a_U,
                      const Box&        a_box,
                      IntVectSet& a_tags)
{  
  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);

  BoxIterator bit( a_box ); RealVect cv_coord,cs_coord;

  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);
    
    int reg = Region(iv,0);
    Real theta1 = (90.0-(m_sunTILT + 5.0))*d_PI_180;
    Real theta2 = (90.0+(m_sunTILT + 5.0))*d_PI_180;
    
    Real th_axis = 20.0;
    Real theta3 = (th_axis)*d_PI_180;
    Real theta4 = (180.0-th_axis)*d_PI_180;
    
    //Real phi    = 6.0;
    Real phi    = 15.0;
    Real phi1   = (90.0  - phi)*d_PI_180;
    Real phi2   = (270.0 + phi)*d_PI_180;
    
    Real lsHCS = fabs(a_U(iv,m_ls_indices.m_iHCS));
    
    if (m_level == 0)
    {        
              
      if  ( ((reg == 3) || (reg == 2)) && (lsHCS<0.9) &&            
             (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
            ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }                
              
    }    
    if (m_level == 1)
    {  
      theta1 = theta3;
      theta2 = (90.0+(m_sunTILT + 20.0))*d_PI_180;
    
      phi    = 35.0;
      phi1   = (90.0  - phi)*d_PI_180;
      phi2   = (270.0 + phi)*d_PI_180;
    
      if  (  (reg == 2) && (lsHCS<0.9) &&            
             (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
            ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
    
      //theta3 = (30.0)*d_PI_180;
      //theta4 = (90.0+(m_sunTILT + 10.0))*d_PI_180;
      //if ((reg == 2) && 
      //    (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
      //    ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }
    if (m_level == 2)
    {
      theta1 = (35.0)*d_PI_180;
      theta2 = (90.0+(m_sunTILT + 20.0))*d_PI_180;
    
      phi    = 35.0;
      phi1   = (90.0  - phi)*d_PI_180;
      phi2   = (270.0 + phi)*d_PI_180;
    
      if  (  (reg == 2) && (lsHCS<0.9) &&            
             (cv_coord[0]>107.0) &&
             (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
            ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
    }

    if (m_level == 3)
    {
      Real phi    = 45.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-(m_sunTILT + 0.5))*d_PI_180; 
      if (cv_coord[0]>100.0) theta3 = (90.0-(m_sunTILT + 2.5))*d_PI_180; 
      if (cv_coord[0]>110.0) theta3 = (90.0-(m_sunTILT + 5.0))*d_PI_180; 
      theta4 = (90.0+(m_sunTILT + 0.5))*d_PI_180;
      if ((reg == 2) && //(cv_coord[0]>90.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }

    if (m_level == 4)
    {
      Real phi    = 60.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-20.0)*d_PI_180;         
      theta4 = (90.0+20.0)*d_PI_180;
      if ((reg == 2) && (cv_coord[0]>108.0) && (cv_coord[0]<117.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    } 
    
  }
}

  
//// Tags for spherical grid
//void HeliosphericProblem::tagCellsSpherical(const FArrayBox&  a_U,
//                          const Box&        a_box,
//                          IntVectSet& a_tags)
//{
//  
//  
//  BaseFab<int> Region;
//  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
//  m_eqSys->stateToPrim(W, a_U, a_box);
//  m_csh->transCartesianVectToCurv(W,a_box,m_level);
//  defineRegions(W, dummyFab, Region, a_box);
//
//  BoxIterator bit( a_box ); 
//  RealVect cv_coord,cs_coord;
//  IntVectSet hp_tags; 
//
//
//  bool reg1_box = false;
//  bool reg2_box = false;
//  bool reg3_box = false;
//  bool reg4_box = false;
//  for( bit.begin(); bit.ok(); ++bit )
//  {
//    const IntVect& iv = bit();
//    int reg = Region(iv,0);
//    reg1_box |= (reg==1);
//    if (reg1_box) return;
//
//    reg2_box |= (reg==2);
//    reg3_box |= (reg==3);
//    reg4_box |= (reg==4);
//
//  }
//
//
//
//
//
////  FArrayBox modGradB(a_box,1);
///*  D.setVal(0.0);
//  D.copy(W,WRHO,0, 1);
//
//
//
//
//  modGradD.setVal(0.0);
//
//
//  FORT_GETRELGRAD_SPHERICAL(
//          CHF_FRA1(modGradD,0),
//          CHF_CONST_FRA1(D,0),            
//          CHF_BOX(a_box),
//          CHF_CONST_INT(m_level));    
//*/  
//
//  for( bit.begin(); bit.ok(); ++bit )
//  {
//    const IntVect& iv = bit();
//    m_csh->getCellCenter(cv_coord,iv,m_level);
//    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);
//    
//    int reg = Region(iv,0);
//    Real rDist = cv_coord[0];   
//    
//    if ((m_level == 0) )
//    {   
//      
//      if ( (reg == 2)  ) 
//      {
//        a_tags |= iv;
////        continue;
//      }
//    } 
//     
///*      Real mgb = modGradB(iv,0);
//      
//      if ((reg == 4) && (mgb > 3.0) )
//      {
//        a_tags |= iv;
//        continue;
//      }*/
//      
//  
//
//
//   bool frontHP = true;
////    Real fromtHPdist = 10000.0;
////    if (m_level >=2) fromtHPdist = 300.0;
//
//
//    Real stopHP_refinement = 10000.0;
//    if (m_ls_indices.m_iRegTr > 0)
//    {
//      Real reg_tr = a_U.get(iv,m_ls_indices.m_iRegTr);
//      if (fabs(reg_tr)<0.8)
//      {
//
//        if ((m_level == 1) && (rDist < stopHP_refinement)) hp_tags |= iv;
// 
////        pout() << "here" <<  endl;
//
////        if (frontHP == true)
////        {
////          if ((dist2 < fromtHPdist*fromtHPdist ) && (x > 0.0))
////          {
////            hp_tags |= iv;
////          }
////        } else
////        {
////          if (dist2 < stopHP_refinement*stopHP_refinement) hp_tags |= iv;
////        }
//      }
//    }
//
//
//
//  }
//
//
//  hp_tags.grow(2);
//  a_tags |= hp_tags;
//
//}



// Tags for spherical grid (FF)
void HeliosphericProblem::tagCellsSpherical(const FArrayBox&  a_U,
                          const Box&        a_box,
                          IntVectSet& a_tags)
{


  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);

  BoxIterator bit( a_box );
  RealVect cv_coord,cs_coord;
  IntVectSet hp_tags;


  bool reg1_box = false;
  bool reg2_box = false;
  bool reg3_box = false;
  bool reg4_box = false;
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    int reg = Region(iv,0);
    reg1_box |= (reg==1);
    if (reg1_box) return;

    reg2_box |= (reg==2);
    reg3_box |= (reg==3);
    reg4_box |= (reg==4);

  }


  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);

    int reg = Region(iv,0);
    Real rDist = cv_coord[0];
    bool frontHP = false;
    Real frontHPdist = 200.0;
    Real stopHP_refinement = 1000000.0;


    if (m_level == 0 )
    {

      if ( reg == 2  )
      {
        a_tags |= iv;
      }

      if ((frontHP == true) && (reg == 4) && (rDist < frontHPdist )
         && (cv_coord[1] <= 3.14 || cv_coord[1] >= 2.36 ))
      {
        hp_tags |= iv;
      }

    }



   if (m_ls_indices.m_iRegTr > 0)
   {
     Real reg_tr = a_U.get(iv,m_ls_indices.m_iRegTr);
     if (fabs(reg_tr)<0.7)
     {

       if ((m_level >= 1) && (rDist < stopHP_refinement)) hp_tags |= iv;

     }
   }



  hp_tags.grow(2);
  a_tags |= hp_tags;

  }

}



  
// Tags HCS region between V1 and V2
void HeliosphericProblem::tagCellsTiltLSV1V2(const FArrayBox&  a_U,
                          const Box&        a_box,
                          IntVectSet& a_tags)
{
  
  
  BaseFab<int> Region;
  FArrayBox dummyFab,W(a_box, m_eqSys->numPrimitives());      
  m_eqSys->stateToPrim(W, a_U, a_box);
  m_csh->transCartesianVectToCurv(W,a_box,m_level);
  defineRegions(W, dummyFab, Region, a_box);

  BoxIterator bit( a_box ); RealVect cv_coord,cs_coord;
  
  Real V1phi = 5.76; // 6.260588792;
  Real V2phi = 0.656506916; // 0.6302118174;
  
  bool reg1_box = false;
  bool reg2_box = false;
  bool reg3_box = false;
  bool bppos = false;
  bool bpneg = false;
  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();        
    int reg = Region(iv,0);
    Real bp = W(iv,WBP);
    reg1_box |= (reg==1);    
    if (reg1_box) return;
    
    reg2_box |= (reg==2);    
    reg3_box |= (reg==3);    
    
    if (bp>0) bppos = true;
    else if (bp<0) bpneg = true;
    
  }
  
  FArrayBox B(a_box,1);
  FArrayBox modGradB(a_box,1);
  B.setVal(0.0);
  modGradB.setVal(0.0);

  int iBGN = WBX;
  int iEND = WBZ;

  FORT_GETVECTMAGNITUDE(
            CHF_FRA1(B,0),
            CHF_CONST_FRA(a_U),
            CHF_CONST_INT(iBGN),
            CHF_CONST_INT(iEND),
            CHF_BOX(a_box));
                
  
  FORT_GETRELGRAD_SPHERICAL(
          CHF_FRA1(modGradB,0),
          CHF_CONST_FRA1(B,0),            
          CHF_BOX(a_box),
          CHF_CONST_INT(m_level));    
  

  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);
    
    int reg = Region(iv,0);
    Real theta1 = (90.0-(m_sunTILT + 2.0))*d_PI_180;
    Real theta2 = (90.0+(m_sunTILT + 2.0))*d_PI_180;
    
    Real th_axis = 20.0;
    Real theta3 = (th_axis)*d_PI_180;
    Real theta4 = (180.0-th_axis)*d_PI_180;
    
    //Real phi    = 6.0;
    Real phi    = 15.0;
    Real phi1   = (90.0  - phi)*d_PI_180;
    Real phi2   = (270.0 + phi)*d_PI_180;
    
    //Real lsHCS = fabs(a_U(iv,m_ls_indices.m_iHCS));
    
    Real lsHCS = 1.0;
    
    Real V1b   = V1phi - 10.0*d_PI_180;
    Real V2b   = V2phi + 10.0*d_PI_180;
    
    Real Biv = m_Bref*B(iv,0);
            
    
    if ((m_level == 0) )
    {   
      
      if ( (reg == 3) && 
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      if ( (reg == 2) && (reg3_box) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      theta3 = (90.0-45)*d_PI_180;
      theta4 = (90.0+45)*d_PI_180;
      
      Real mgb = modGradB(iv,0);
      
      if ((reg == 2) && (mgb > 3.0) && (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && (cv_coord[0] < 160.0) &&
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      {
        a_tags |= iv;
        continue;
      }
      
      //if ((reg == 2) && bppos && bpneg && (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
      //    ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      //{
      //  a_tags |= iv;
      //  continue;
      //}
      
                
      if  ( (reg == 2) && (lsHCS<0.8) &&            
             (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
            ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }                              
    }    
    
    if (m_level == 1)
    {   
      
      if ( (reg == 3) && (cv_coord[0]>70.0) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;      
      }
      
      if ( (reg == 2) && (reg3_box) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      theta3 = (90.0-45)*d_PI_180;
      theta4 = (90.0+45)*d_PI_180;
                  
      Real mgb = modGradB(iv,0);
      if ( (reg == 2) && (mgb > 3.0) && (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && (cv_coord[0] < 160.0) &&
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      {
        a_tags |= iv;
        continue;
      }
      
      if ( (reg == 2) && (Biv < 2.0) && (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && (cv_coord[0] < 170.0) &&
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      {
        a_tags |= iv;
        continue;
      }
                                               
    }    
    
    if (m_level == 2)
    {   
            
      
      theta3 = (90.0-35)*d_PI_180;
      theta4 = (90.0+33)*d_PI_180;
      
      V1b   = V1phi - 10.0*d_PI_180;
      V2b   = V2phi +  0.0*d_PI_180;

      
      bool level2_criteria = ((cv_coord[0] > 95.0) && (cv_coord[0] < 120.0));
      
      Real mgb = modGradB(iv,0);
      if ((level2_criteria) && (reg == 2) && (mgb > 3.0) && (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && (cv_coord[0] < 160.0) &&
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      {
        a_tags |= iv;
        continue;
      }
      
      if ((level2_criteria) && (reg == 2) && (Biv < 2.0) && (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && (cv_coord[0] < 170.0) &&
          ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)))
      {
        a_tags |= iv;
        continue;
      }
                                               
    }    
    
    
    
    /*if (m_level == 1)
    {  
      th_axis = 30.0;
      theta3 = (th_axis)*d_PI_180;
      theta4 = (180.0-th_axis)*d_PI_180;
      
      Real c[3]={-30,-5,0};
      
      Real r = sqrt(
        (cs_coord[0]-c[0])*(cs_coord[0]-c[0])+
        (cs_coord[1]-c[1])*(cs_coord[1]-c[1])+
        (cs_coord[2]-c[2])*(cs_coord[2]-c[2]));
      
      if  ( (reg == 2) && (lsHCS<0.8) && (r<140.0) &&   
             (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
            ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
    }
    if (m_level == 2)
    {
      th_axis = 35.0;
      theta3 = (th_axis)*d_PI_180;
      theta4 = (180.0-th_axis)*d_PI_180;
      
      Real c[3]={-75,-20,21};
      
      Real r = sqrt(
        (cs_coord[0]-c[0])*(cs_coord[0]-c[0])+
        (cs_coord[1]-c[1])*(cs_coord[1]-c[1])+
        (cs_coord[2]-c[2])*(cs_coord[2]-c[2]));
      
      if  ( (reg == 2) && (lsHCS<0.8) && (r>200.0) &&   
             (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
            ((cv_coord[1]<V2b) || (cv_coord[1]>V1b)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
    }

    if (m_level == 3)
    {
      Real phi    = 45.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-(m_sunTILT + 0.5))*d_PI_180; 
      if (cv_coord[0]>100.0) theta3 = (90.0-(m_sunTILT + 2.5))*d_PI_180; 
      if (cv_coord[0]>110.0) theta3 = (90.0-(m_sunTILT + 5.0))*d_PI_180; 
      theta4 = (90.0+(m_sunTILT + 0.5))*d_PI_180;
      if ((reg == 2) && //(cv_coord[0]>90.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    }

    if (m_level == 4)
    {
      Real phi    = 60.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;        
          
      theta3 = (90.0-20.0)*d_PI_180;         
      theta4 = (90.0+20.0)*d_PI_180;
      if ((reg == 2) && (cv_coord[0]>108.0) && (cv_coord[0]<117.0) &&
          (cv_coord[2]>theta3) && (cv_coord[2]<theta4) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) a_tags |= iv;
      
    } */
    
  }

}


//                             Check geometrical limitations for grid adaptation
void HeliosphericProblem::lockedCellsRegrid( BaseFab<int> & a_flag,
                                       const FArrayBox&  a_U,
                                       const Box&     a_box)
{
  a_flag.setVal(0);
      
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)
  {
    if( m_adaptRmax <= 0.0 ) return;

    Real D_DECL(x,y,z), dist2, adaptR2max;    
    Real dx =  m_csh->dx(0,m_level);
    
    adaptR2max = m_adaptRmax*m_adaptRmax;
    
    BoxIterator bit(a_box);
    for (bit.begin(); bit.ok(); ++bit)
    {
      const IntVect& iv = bit();                      
      D_EXPR(
        x = (iv[0] + 0.5)*dx - m_sunXYZ[0],
        y = (iv[1] + 0.5)*dx - m_sunXYZ[1],
        z = (iv[2] + 0.5)*dx - m_sunXYZ[2]  );
        
      dist2 = D_TERM( x*x, + y*y, + z*z );

      if( dist2 > adaptR2max ) a_flag(iv,0) = 1;  
    }
  }
  
/*  BoxIterator bit( a_box ); RealVect cv_coord,cs_coord;

  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    m_csh->getCellCenter(cv_coord,iv,m_level);
    m_csh->transCurvCoordsToCartesian(cs_coord,cv_coord);
    
          
    if ((cs_coord[0]>-450.0)) a_flag(iv,0) = 1;  
    
  }*/
  
  
    
} 

Real HeliosphericProblem::computeDt( const FArrayBox& a_U,
                               const FArrayBox& a_dt,
                               const Box&     a_box,
                               IntVect&       a_minDtCell)
{
    
  //return MultiFluidProblem::computeDt( a_U, a_dt, a_box, a_minDtCell);    
  
  if (m_csh->coordinateSystem() != CoordinateSystemHandler::CS_Cartesian)
  {
    return MultiFluidProblem::computeDt( a_U, a_dt, a_box, a_minDtCell);    
  }
  
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
    
    //if (dist2 <= R0dtIgnore2 )
    //{ Needed to create breal point
    //  dt = a_dt.get(iv,0);
    //}
    
    if ((dist2 > R0dtIgnore2 ) && (dt < dtMin))
    {
      dtMin = dt;
      a_minDtCell = iv;
    }      
  }
  
  return dtMin;
}
                           

//                            Return boundary condition flags for all boundaries
void HeliosphericProblem::getBCFlags( eBoundaryConditions leftBC,
                                      eBoundaryConditions rightBC,
                                      eBoundaryConditions bottomBC,
                                      eBoundaryConditions topBC,
                                      eBoundaryConditions frontBC,
                                      eBoundaryConditions behindBC )
{
  leftBC   = BC_Continuous;
  rightBC  = BC_Fixed;
  bottomBC = BC_Continuous;
  topBC    = BC_Continuous;
  frontBC  = BC_Continuous;
  behindBC = BC_Continuous;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions HeliosphericProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    return ( a_sd == Side::Lo ) ? BC_Continuous : BC_Fixed;
  } else {
    return BC_Continuous;
  }
}

int HeliosphericProblem::lsIndexField(int a_s)
{
  //CH_assert((a_s==0)||(a_s==1));
  return WVELX;  
}


void HeliosphericProblem::prepareForH5Writing()
{
  int curProc = procID();      

  if (m_level>0) return;
  
  Real r0 = m_concludeR0;
  m_csh->getClosestIndex(m_concludeI,r0,0,m_level);  
  if (m_concludeI<0) MayDay::Error("mhdam.concludeR0 is outside of the problem domain");
          
  //ir0 = dBox.bigEnd()[0]-2;

  RealVect rvr0;
  m_csh->getCellCenter(rvr0,IntVect(D_DECL(m_concludeI,0,0)),m_level);
  r0 = rvr0[0];        
  m_concludeR0 = r0;
  
  hid_t h5f,dataspace,attr,dataset;
  hsize_t dimsf;herr_t status;  
  
  //H5E_auto_t efunc; void* edata; // turn auto error messaging off

#ifdef H516
  H5E_auto_t efunc; void* edata;
  H5Eget_auto(&efunc, &edata);
  H5Eset_auto(NULL, NULL);
  H5Eset_auto(efunc, edata);
#else
  H5E_auto2_t efunc; void* edata;
  H5Eget_auto2(H5E_DEFAULT, &efunc, &edata);
  H5Eset_auto2(H5E_DEFAULT, NULL, NULL);
  H5Eset_auto2(H5E_DEFAULT, efunc, edata);
#endif

  h5f = H5Fopen(m_writeH5File.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  
  if (h5f<0)
  {
    if (curProc == 0)
    {
      h5f = H5Fcreate(m_writeH5File.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      m_csh->writeGeomInfo(h5f);        
      
      dimsf = 1;    
      dataspace = H5Screate_simple(1, &dimsf, NULL);    

#ifdef H516      
  attr = H5Acreate(h5f, "r0", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(h5f, "r0", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#endif

      H5Awrite(attr, H5T_NATIVE_DOUBLE, &m_concludeR0);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);
      
      dimsf = 1;
      dataspace = H5Screate_simple(1, &dimsf, NULL);  

#ifdef H516        
  attr = H5Acreate(h5f, "ind_r0", H5T_NATIVE_INT, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(h5f, "ind_r0", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#endif


      H5Awrite(attr, H5T_NATIVE_INT, &m_concludeI);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);                                                     
      
      const Box & dBox = m_domain.domainBox();   
      int domain[2]={dBox.size(1),dBox.size(2)};
      dimsf = 2;
      dataspace = H5Screate_simple(1, &dimsf, NULL);   

#ifdef H516       
  attr = H5Acreate(h5f, "domain", H5T_NATIVE_INT, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(h5f, "domain", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); 
#endif

      H5Awrite(attr, H5T_NATIVE_INT, domain);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);                                                     
      
      int nComp = UNUM+m_eqSys->numTrackingSurfaces();
      dimsf = 1;
      dataspace = H5Screate_simple(1, &dimsf, NULL);   

#ifdef H516       
  attr = H5Acreate(h5f, "num_components", H5T_NATIVE_INT, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(h5f, "num_components", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); 
#endif

      H5Awrite(attr, H5T_NATIVE_INT, &nComp);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);
      
      int nDataSets = 0;
      dimsf = 1;
      dataspace = H5Screate_simple(1, &dimsf, NULL);   

#ifdef H516       
  attr = H5Acreate(h5f, "num_datasets", H5T_NATIVE_INT, dataspace, H5P_DEFAULT); 
#else
  attr = H5Acreate2(h5f, "num_datasets", H5T_NATIVE_INT, dataspace, H5P_DEFAULT, H5P_DEFAULT); 
#endif

      H5Awrite(attr, H5T_NATIVE_INT, &nDataSets);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);
            
      hsize_t chunksize = 16384;
      hsize_t maxdim = H5S_UNLIMITED;
      dimsf = 1;
      dataspace    = H5Screate_simple(1, &dimsf, &maxdim);   
      hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk ( cparms, 1, &chunksize);      // chuncked dataset to be able to add elements

#ifdef H516
  dataset = H5Dcreate(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
  dataset = H5Dcreate2(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

//      dataset   = H5Dcreate2(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);
      H5Sclose(dataspace);                                                               
      H5Dclose(dataset);                                                               
      
      H5Fclose(h5f);
    }
    
    m_lastWrittenTime = -1.0;
    m_lastInd = -1;    
  } else
  {    
    int nDataSets = 0;
    attr = H5Aopen(h5f,"num_datasets",H5P_DEFAULT);
    H5Aread(attr, H5T_NATIVE_INT, &nDataSets);
    H5Aclose(attr);        
    
    Real time = 0.0;
    ParmParse parser("mhdam");      
    if (parser.contains("restart_file"))
    {
      std::string chkFile;      
      parser.query("restart_file",chkFile);
      hid_t h5fchk = H5Fopen(chkFile.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
      attr = H5Aopen(h5fchk,"time",H5P_DEFAULT);
      H5Aread(attr, H5T_NATIVE_DOUBLE, &time);
      H5Aclose(attr);        
      H5Fclose(h5fchk);
    }
        
    if (nDataSets > 0)
    {

#ifdef H516
  dataset = H5Dopen(h5f, "datasets_time"); 
#else
  dataset = H5Dopen2(h5f, "datasets_time", H5P_DEFAULT);
#endif  

      dataspace = H5Dget_space(dataset);    /* dataspace handle */      
      H5Sget_simple_extent_dims(dataspace, &dimsf, NULL);
      if (dimsf!=nDataSets) MayDay::Warning("datasets_time size and number of data sets are different");

      Real * datasets_time = new Real[dimsf];      
      status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasets_time);      
      H5Dclose(dataset);     
      H5Sclose(dataspace);
      
      Real Tref      = eos_AU/(m_lismV*3600.0*24.0);
      Real timeDays  = (time-m_writeT0)*Tref; // days      
      
      Real * res = std::lower_bound(datasets_time, datasets_time + nDataSets,timeDays);
      int   iRes = res - datasets_time;
      CH_assert((iRes>=0)&&(iRes<=nDataSets));
      
      m_lastInd = iRes - 1;      
      if (m_lastInd>=0)
      {
        char dataspace_name[20];
        sprintf(dataspace_name,"data%d",m_lastInd);

#ifdef H516
  dataset = H5Dopen(h5f, dataspace_name);     
  attr = H5Aopen_name(dataset,"time");
#else
  dataset = H5Dopen2(h5f, dataspace_name, H5P_DEFAULT);     
  attr = H5Aopen_by_name(dataset, ".", "time", H5P_DEFAULT, H5P_DEFAULT);
#endif

        H5Aread(attr, H5T_NATIVE_DOUBLE, &m_lastWrittenTime);      
        H5Aclose(attr);  
        H5Dclose(dataset);     
        
        m_lastWrittenTime = m_writeT0 + m_lastWrittenTime*(m_lismV*3600.0*24.0)/eos_AU;
        CH_assert(m_lastWrittenTime<=time);        
      } else
      {
        m_lastWrittenTime = -1.0;
      }
      if ((m_lastInd<nDataSets-1) && (procID() == 0))
        pout() << "slices starting from " << m_lastInd+1 << " will be ovewritten" << endl;
    
    } else
    {
      m_lastWrittenTime = -1.0;
      m_lastInd = -1;    
    }  
    H5Fclose(h5f);
  }      
}


void HeliosphericProblem::writeSphericalSlice(const LevelData<FArrayBox>& a_U, Real a_time)
{
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();  
  
  if (CoordinateSystem != CoordinateSystemHandler::CS_Spherical) return;
  if (m_level>0) return;
  
  if (abs(m_verbosity) >= 2 ) pout() << "Enter HeliosphericProblem::writeSphericalSlice" << endl;
        
  
  const Box & dBox = m_domain.domainBox();       
                           
  
  Vector<Box> vstrips;
  Vector<int> procs;
  Box b;

  b.define( IntVect(D_DECL(m_concludeI, dBox.smallEnd()[1], dBox.smallEnd()[2])), 
            IntVect(D_DECL(m_concludeI, dBox.bigEnd()[1],   dBox.bigEnd()[2])));                         
  vstrips.push_back(b);                                         
  procs.push_back(0);      
  
  DisjointBoxLayout strips(vstrips, procs, m_domain);  
  
  LevelData<FArrayBox> lstrips(strips, UNUM+m_eqSys->numTrackingSurfaces());    
  
  a_U.copyTo(Interval(URHO,UBZ), lstrips, Interval(URHO,UBZ));
  
  if (m_eqSys->numTrackingSurfaces()>0)
    a_U.copyTo(m_eqSys->lvlsStateInterval(), lstrips, Interval(UBZ+1,UBZ+m_eqSys->numTrackingSurfaces()));
    
  if (m_lastInd < 0) m_lastInd = 0;
  else m_lastInd++;
        
  int curProc = procID();      
  if (curProc==0)
  {
    const FArrayBox & U = lstrips[lstrips.dataIterator()];
    FArrayBox W(U.box(),U.nComp());
    
    FORT_CONSTOPRIM(CHF_FRA(W),
                    CHF_CONST_FRA(U),
                    CHF_BOX(b));    
                    
    if (m_eqSys->numTrackingSurfaces()>0)
      W.copy(U,b,UBZ+1,b,UBZ+1,m_eqSys->numTrackingSurfaces());
      
            
              
    int vects[2]={WVELX,WBX};
    m_csh->transCartesianVectToCurv(W,vects,2,b,m_level);
            
    W.mult(m_lismN, WRHO);
    W.mult(m_lismV, WVELX, 3);  
    W.mult(m_Bref*m_Bref, WPRES);  
    W.mult(m_Bref , WBX  , 3);
    
    FArrayBox WRot(W.box(),W.nComp());
    if (m_writeRotating ==  true)
      FORT_TRANSFORMTOROTATINGFRAME(CHF_FRA(WRot),
                    CHF_CONST_FRA(W),
                    CHF_CONST_REAL(m_writeT0),
                    CHF_CONST_REAL(a_time),
                    CHF_CONST_INT(m_level),
                    CHF_BOX(b));
    else
      WRot.copy(W);
    
        
    
    hsize_t dimsf;
    hid_t dataspace,attr,dataset;herr_t status;

    
    RealVect cv_coord;
    m_csh->getCellCenter(cv_coord,b.smallEnd(),m_level);

#ifdef CH_USE_HDF5
    hid_t h5f = H5Fopen(m_writeH5File.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
    
    
    int nDataSets = m_lastInd + 1;
    attr = H5Aopen(h5f,"num_datasets",H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, &nDataSets);   
    H5Aclose(attr);    
           
    char dataspace_name[30];
    sprintf(dataspace_name,"data%i",m_lastInd);
    dimsf = b.numPts()*WRot.nComp();  
    
    Real timeDays = (a_time - m_writeT0)*eos_AU/(m_lismV*3600.0*24.0);
    
    if (H5Lexists(h5f, dataspace_name,H5P_DEFAULT) > 0)
    {
      pout() << "overwriting slice " << m_lastInd << endl;
      // rewrite current dataset

#ifdef H516
  dataset   = H5Dopen(h5f, dataspace_name);  
#else
  dataset   = H5Dopen2(h5f, dataspace_name, H5P_DEFAULT);  
#endif

      status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, WRot.dataPtr());            
            
      dimsf = 1;          
      attr = H5Aopen(dataset, "time", H5P_DEFAULT); 
      H5Awrite(attr, H5T_NATIVE_DOUBLE, &timeDays);   
      H5Aclose(attr);                                                           
      H5Dclose(dataset);      
    } else
    {  
      dataspace = H5Screate_simple(1, &dimsf, NULL);   

#ifdef H516
  dataset   = H5Dcreate(h5f, dataspace_name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
  dataset   = H5Dcreate2(h5f, dataspace_name, H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
  
      status    = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, WRot.dataPtr());      
      H5Sclose(dataspace);    
            
      dimsf = 1;    
      dataspace = H5Screate_simple(1, &dimsf, NULL);    

#ifdef H516      
  hid_t attr = H5Acreate(dataset, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
  hid_t attr = H5Acreate2(dataset, "time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT);
#endif

      H5Awrite(attr, H5T_NATIVE_DOUBLE, &timeDays);   
      H5Aclose(attr);                                                     
      H5Sclose(dataspace);                                                         
      H5Dclose(dataset);
    }
    
    dimsf = 1; // data space for one real number
    dataspace = H5Screate_simple(1, &dimsf, NULL);          
    dimsf = nDataSets;

#ifdef H516
  dataset = H5Dopen(h5f, "datasets_time");  
#else
  dataset = H5Dopen2(h5f, "datasets_time", H5P_DEFAULT);  
#endif
 
    H5Dset_extent(dataset, &dimsf);
    hid_t filespace = H5Dget_space(dataset);
    hsize_t offset = m_lastInd;        
    hsize_t count  = 1;
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, &offset, NULL, &count, NULL);      
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, dataspace, filespace, H5P_DEFAULT, &timeDays);        
    H5Dclose (dataset);
    H5Sclose (dataspace);
    H5Sclose (filespace);
    
    
    H5Fclose(h5f);
#endif

    if (0)
    {
      NodeFArrayBox coords(WRot.box(),3);          
      m_csh->getNodeCoordsCartesian(coords, m_level);         
      NodeFArrayBox coords_cv(WRot.box(),3);          
      m_csh->getNodeCoords(coords_cv, m_level);         
          
      FILE* tfile = OpenTecplotFile(strcat(dataspace_name,".dat"),W.nComp()+SpaceDim);    
      WriteFArrayBoxToTecplotFile(tfile, WRot, WRot.box(), Interval(0, WRot.nComp()-1), coords.getFab(), coords_cv.getFab()); 
      CloseTecplotFile(tfile);      
    }
    
  }  
  m_lastWrittenTime = a_time;
}

/// Things to do before after calculations are complete
void HeliosphericProblem::conclude(const LevelData<FArrayBox>& a_U,
                        Real a_time,
                        int  a_cur_step,                        
                        const std::string & a_inputfile)
{
  //return;
  
  if (m_writeH5 == false) return;
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();  
  
  if (CoordinateSystem != CoordinateSystemHandler::CS_Spherical) return;
  if (m_level>0) return;
  
  // We did not calculate anything, just want to write spherical distribution
  if ((a_time >=m_writeT0) && (m_lastInd < 0))    
    writeSphericalSlice(a_U, a_time);
  
  int curProc = procID();      
  if (curProc!=0) return;
  
  int i;char dataspace_name[20];
  
  
#ifdef CH_USE_HDF5
  hid_t h5f = H5Fopen(m_writeH5File.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);

  hsize_t dimsf;      
  hsize_t chunksize = 16384;
  hsize_t maxdim = H5S_UNLIMITED;
  
  hid_t   attr,dataset,dataspace;
  
  int nDataSets = 0;
  attr = H5Aopen(h5f,"num_datasets",H5P_DEFAULT);
  H5Aread(attr, H5T_NATIVE_INT, &nDataSets);
  H5Aclose(attr);    
  
  if ((nDataSets > 0) && (0))
  {
    Real * dsTimes = new Real[nDataSets];
    for (i=0;i<nDataSets;i++)
    {        
      sprintf(dataspace_name,"data%d",i);
  
#ifdef H516
  dataset = H5Dopen(h5f, dataspace_name);   
#else
  dataset = H5Dopen2(h5f, dataspace_name, H5P_DEFAULT);   
#endif
      attr = H5Aopen(dataset,"time",H5P_DEFAULT);
      H5Aread(attr, H5T_NATIVE_DOUBLE, &dsTimes[i]);
      H5Aclose(attr);
      H5Dclose(dataset);     
    }
    
//    H5E_auto_t efunc; void* edata; // turn auto error messaging off
//    H5Eget_auto(&efunc, &edata);
//    H5Eset_auto(NULL, NULL);
    //H5Ldelete(h5f, "datasets_time",H5P_DEFAULT);
    //dataset   = H5Dopen2(h5f, "datasets_time", H5P_DEFAULT);  
//    H5Eset_auto(efunc, edata);

    dimsf = nDataSets;
    if (H5Lexists(h5f, "datasets_time",H5P_DEFAULT) > 0)
    {
#ifdef H516
  dataset = H5Dopen(h5f, "datasets_time");  
#else
  dataset = H5Dopen2(h5f, "datasets_time", H5P_DEFAULT);  
#endif
 
      H5Dset_extent(dataset, &dimsf);
      herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dsTimes);        
    } else
    {              
      dataspace    = H5Screate_simple(1, &dimsf, &maxdim);   
      hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
      H5Pset_chunk ( cparms, 1, &chunksize);
      
#ifdef H516
    dataset   = H5Dcreate(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
#else
    dataset   = H5Dcreate2(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

//      dataset   = H5Dcreate2(h5f, "datasets_time", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT, cparms, H5P_DEFAULT);  
      herr_t status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dsTimes);        
      H5Sclose(dataspace);                                                               
    }
    H5Dclose(dataset);    
    delete[] dsTimes;
  }  
  
  if (H5Lexists(h5f, "input",H5P_DEFAULT) == 0)
  {
    hid_t s_type = H5Tcopy(H5T_C_S1);
    H5Tset_size(s_type, a_inputfile.length()); //extra requirement for strings
    hid_t aid  = H5Screate(H5S_SCALAR);

#ifdef H516
    hid_t problemdataset   = H5Dcreate(h5f, "input",  s_type, aid, H5P_DEFAULT);
#else
    hid_t problemdataset   = H5Dcreate2(h5f, "input",  s_type, aid, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    if(problemdataset < 0) MayDay::Error("HeliosphericProblem::conclude");  

    char* tmp = (char*)a_inputfile.c_str();
    H5Dwrite(problemdataset, s_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, tmp);

    H5Sclose(aid);
    H5Tclose(s_type);
    H5Dclose(problemdataset);
  }
       
  H5Fclose(h5f);
#endif
  
}

void HeliosphericProblem::shockBC( const FArrayBox & a_W,
                                         FArrayBox & a_U,
                                      BaseFab<int> & a_REG,
                                         const int   a_level)
{

  Box WBox = a_REG.box();
  WBox.grow(-1);
  BaseFab<int> boundary(WBox,1);
  boundary.setVal(0);

  FORT_SHOCKBOUNDARY( CHF_CONST_FIA1(a_REG,0),
                      CHF_FIA1(boundary,0),
                      CHF_BOX(WBox));

  BaseFab<Real> normal(WBox,4);
  normal.setVal(0);

//  FORT_SHOCKNORMAL( CHF_CONST_FRA(a_W),
//                    CHF_BOX(WBox),
//                    CHF_FRA(normal),
//                    CHF_CONST_FIA1(boundary,0));

  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoPIW = eqSys->getPickupIons()->primInterval().begin();
  int iRhoPIU = eqSys->getPickupIons()->consInterval().begin();

//changed from ZANKSCATTER2 to ZANKSCATTER
//  FORT_ZANKSCATTER( CHF_CONST_FRA(a_W),
//                     CHF_FRA(a_U),
//                      CHF_CONST_INT(iRhoPIW),
//                      CHF_CONST_INT(iRhoPIU),
//                      CHF_CONST_FIA1(boundary,0),
//                      CHF_CONST_FRA(normal),
//                      CHF_BOX(WBox),
//                      CHF_CONST_INT(a_level));

/*
//RKB
// Using VADIMBC with approximated functions
  FORT_VADIMBC( CHF_CONST_FRA(a_W),
                     CHF_FRA(a_U),
                      CHF_CONST_INT(iRhoPIW),
                      CHF_CONST_INT(iRhoPIU),
                      CHF_CONST_FIA1(boundary,0),
                      CHF_CONST_FRA(normal),
                      CHF_BOX(WBox),
                      CHF_CONST_INT(a_level));


*/
 
//RKB--> Using Table from hybrid Sim directly
  FORT_VADIMBC_TABLE( CHF_CONST_FRA(a_W),
                     CHF_FRA(a_U),
                      CHF_CONST_INT(iRhoPIW),
                      CHF_CONST_INT(iRhoPIU),
                      CHF_CONST_FIA1(boundary,0),
                      CHF_CONST_FRA(normal),
                      CHF_BOX(WBox),
                      CHF_CONST_INT(a_level));



   //exit(1);

}

