#include <iostream>
#include <iomanip>
using std::ifstream;
using std::ios;
#include "DebugOut.H"

#include "LoHiSide.H"

#include "RiemannProblemMF.H"
#include "RiemannProblemMFF_F.H"
#include "RiemannF_F.H"
#include "ChargeExchange2F_F.H"
#include "EqSysMHDMF.H"


#include "LoHiCenter.H"

// Null constructor
RiemannProblemMF::RiemannProblemMF()
{
  m_isFortranCommonSet = false;
}

// Input parameters
void RiemannProblemMF::input( ParmParse & parser, int verbosity )
{
  m_densityL     = 1.0;
  m_temperatureL = 5000.0;  
  m_velxL        = 25000.0;
  m_velyL        = 0.0;
  m_velzL        = 0.0;
  m_BxL          = 0.0;
  m_ByL          = 0.0;
  m_BzL          = 0.0;
  
  m_velyR        = 0.0;
  m_velzR        = 0.0;
  m_BxR          = 0.0;
  m_ByR          = 0.0;
  m_BzR          = 0.0;
  
  
  m_netDenL      = 1.0;
  m_startX       = 0.5;  
  m_gamma        = 1.6666666667;
  m_NInitDistr   = 0;

  parser.query( "gamma",            m_gamma        );
  
  parser.query( "netN",             m_netDenL);
  parser.query( "startX",           m_startX       );
  parser.query( "neutrals_initial_distribution", m_NInitDistr);
    

  
  parser.query( "densityL",         m_densityL      );
  parser.query( "temperatureL",     m_temperatureL );
  parser.query( "velxL",            m_velxL        );
  parser.query( "velyL",            m_velyL        );
  parser.query( "velzL",            m_velzL        );
  parser.query( "BxL",              m_BxL          );
  parser.query( "ByL",              m_ByL          );
  parser.query( "BzL",              m_BzL          );
            
  if (!parser.contains( "densityR"))
  {
    setFortranCommon( m_gamma,
                    m_densityL,
                    m_temperatureL,
                    m_velxL,
                    m_BxL,
                    m_netDenL,
                    m_startX,
                    m_NInitDistr    );
    
  } else
  {
    m_densityR     = m_densityL;
    m_temperatureR = m_temperatureL;    
    m_velxR        = m_velxL;
    m_velyR        = m_velyL;
    m_velzR        = m_velzL;
    m_BxR          = m_BxL;
    m_ByR          = m_ByL;
    m_BzR          = m_BzL;
    m_TMLIM        = 0.5*(m_temperatureR + m_temperatureL);
    
    m_netTemperatureL = m_temperatureL;
    m_netVelx = m_velxL;
    m_netVely = m_velyL;
    m_netVelz = m_velzL;
  
    parser.query( "densityR",         m_densityR      );
    parser.query( "temperatureR",     m_temperatureR );
    parser.query( "velxR",            m_velxR        );
    parser.query( "velyR",            m_velyR        );
    parser.query( "velzR",            m_velzR        );
    parser.query( "BxR",              m_BxR          );
    parser.query( "ByR",              m_ByR          );
    parser.query( "BzR",              m_BzR          );
    
    parser.query( "TMLIM",        m_TMLIM);
        
    parser.query( "netT",  m_netTemperatureL);
    parser.query( "netUX", m_netVelx);
    parser.query( "netUY", m_netVely);
    parser.query( "netUZ", m_netVelz);
    
    
    setFortranCommonLR();
  }
  
  
  
  int nFluids = 3;
  parser.query( "fluids",  nFluids   );

  switch( nFluids ){
  case 2  : m_physModel  = PP_2FluidPM; break;
  case 3  : m_physModel  = PP_3FluidPM; break;
  case 4  : m_physModel  = PP_4FluidPM; break;
  default : m_physModel  = PP_MHDPM;
  }


  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Two fluid 1D Riemann problem with neutrals input:"      << endl;
    pout() << "gamma        = " << m_gamma        << endl;
    pout() << "densityL     = " << m_densityL      << endl;
    pout() << "temperatureL = " << m_temperatureL << endl;
    pout() << "velxL        = " << m_velxL        << endl;
    pout() << "BxL          = " << m_BxL          << endl;
    pout() << "m_netDenL    = " << m_netDenL      << endl;
    pout() << "startX       = " << m_startX       << endl;
    pout() << "neutrals_initial_distribution       = " << m_NInitDistr       << endl;
  }

  
}

// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RiemannProblemMF::setFortranCommon( const Real&     a_gamma,
                                          const Real&     a_densityL,
                                          const Real&     a_temperatureL,
                                          const Real&     a_velxL,
                                          const Real&     a_BxL,
                                          const Real&     a_netnumL,
                                          const Real&     a_startX,
                                          const int&      a_NInitDistr  )
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETRIEMANNMF( CHF_CONST_REAL( a_gamma         ),
                        CHF_CONST_REAL( a_densityL       ),
                        CHF_CONST_REAL( a_temperatureL  ),
                        CHF_CONST_REAL( a_velxL         ),
                        CHF_CONST_REAL( a_BxL           ),
                        CHF_CONST_REAL( a_netnumL       ),
                        CHF_CONST_REAL( a_startX        ),
                        CHF_CONST_INT ( a_NInitDistr) );

  Real scaleLen  = 1.5e+13;
  Real scaleNDen = m_densityL;

  FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( m_velxL   ),
                         CHF_CONST_REAL( scaleLen  ),
                         CHF_CONST_REAL( scaleNDen ) );

    m_isFortranCommonSet = true;
}


// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void RiemannProblemMF::setFortranCommonLR()
{
  CH_assert(m_isFortranCommonSet == false);
  
  
  FORT_SETRIEMANNLR(   CHF_CONST_REAL( m_gamma     ),
                     CHF_CONST_REAL( m_densityL  ),
                     CHF_CONST_REAL( m_densityR  ),
                     CHF_CONST_REAL( m_temperatureL ),
                     CHF_CONST_REAL( m_temperatureR ),
                     CHF_CONST_REAL( m_velxL     ),
                     CHF_CONST_REAL( m_velxR     ),
                     CHF_CONST_REAL( m_velyL     ),
                     CHF_CONST_REAL( m_velyR     ),
                     CHF_CONST_REAL( m_velzL     ),
                     CHF_CONST_REAL( m_velzR     ),
                     CHF_CONST_REAL( m_BxL       ),
                     CHF_CONST_REAL( m_BxR       ),
                     CHF_CONST_REAL( m_ByL       ),
                     CHF_CONST_REAL( m_ByR       ),
                     CHF_CONST_REAL( m_BzL       ),
                     CHF_CONST_REAL( m_BzR       ),
                     CHF_CONST_REAL( m_netDenL   ),
                     CHF_CONST_REAL( m_netTemperatureL ),
                     CHF_CONST_REAL( m_netVelx ),
                     CHF_CONST_REAL( m_netVely ),
                     CHF_CONST_REAL( m_netVelz ),
                     CHF_CONST_REAL( m_startX  ),
                     CHF_CONST_INT ( m_NInitDistr),
                     CHF_CONST_REAL( m_TMLIM  )  );

    

  Real scaleLen  = 1.5e+13;
  Real scaleNDen = m_densityL;

  FORT_SETCHARGEEX_PARS( CHF_CONST_REAL( m_netVelx   ),
                         CHF_CONST_REAL( scaleLen  ),
                         CHF_CONST_REAL( scaleNDen ) );

    m_isFortranCommonSet = true;
}


// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void RiemannProblemMF::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void RiemannProblemMF::copy_PhysProblem(const PhysProblem* a_PP)
{
  const RiemannProblemMF* PP = dynamic_cast<const RiemannProblemMF*>(a_PP);
  
  MultiFluidProblem::copy_PhysProblem(a_PP);
  
  if (PP == NULL) MayDay::Error("RiemannProblemMF::copy_PhysProblem. Wrong argument");
  
  if( PP->m_isFortranCommonSet == true )
  {
    m_gamma        = PP->m_gamma;
    m_densityL   =  PP->m_densityL;
    m_densityR   = PP->m_densityR;
    m_temperatureL  = PP->m_temperatureL;
    m_temperatureR  = PP->m_temperatureR;
    m_velxL      = PP->m_velxL;
    m_velxR      = PP->m_velxR;
    m_velyL      = PP->m_velyL;
    m_velyR      = PP->m_velyR;
    m_velzL      = PP->m_velzL;
    m_velzR      = PP->m_velzR;
    m_BxL        = PP->m_BxL;
    m_BxR        = PP->m_BxR;
    m_ByL        = PP->m_ByL;
    m_ByR        = PP->m_ByR;
    m_BzL        = PP->m_BzL;
    m_BzR        = PP->m_BzR;
    
    m_netDenL         = PP->m_netDenL;
    m_netTemperatureL = PP->m_netTemperatureL;
    m_netVelx         = PP->m_netVelx;
    m_netVely         = PP->m_netVely;
    m_netVelz         = PP->m_netVelz;
    m_NInitDistr      = PP->m_NInitDistr;
    
    m_TMLIM      = PP->m_TMLIM;


    setFortranCommonSet();
  }  
  
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* RiemannProblemMF::new_PhysProblem()
{
  RiemannProblemMF* retval = new RiemannProblemMF();
  
  retval->copy_PhysProblem(this);

  return static_cast<PhysProblem*>(retval);
}

// Set boundary fluxes
void RiemannProblemMF::fluxBC(       FArrayBox&      a_F,
                                      FArrayBox&      a_Bn,
                                const FArrayBox&      a_W,
                                const FArrayBox&      a_Wextrap,
                                const int&            a_dir,
                                const Side::LoHiSide& a_side,
                                const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  Real dx = m_csh->dx(0,m_level);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  

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

      // Shift things to all line up correctly
      boundaryBox.shiftHalf(a_dir,-sign);
      a_F.shiftHalf(a_dir,-sign);

      // Set the boundary fluxes
      FORT_FLUXBCMF( CHF_FRA(a_F),
                      CHF_CONST_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_REAL(dx),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(boundaryBox) );

      // Shift returned fluxes to be face centered
      a_F.shiftHalf(a_dir,sign);
    }
  }
}

// Set up initial conditions
void RiemannProblemMF::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  Real dx = m_csh->dx(0,m_level);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    // Set up initial condition in this grid
    FORT_RIEMANNINITMF( CHF_CONST_FRA(U),
                         CHF_CONST_REAL(dx),
                         CHF_CONST_INT(iRhoN),
                         CHF_CONST_INT(fluids),
                         CHF_BOX(uBox));
  }
}
                                           // Define regions for charge exchange
void RiemannProblemMF::defineRegions( const FArrayBox    & a_W,
                                             FArrayBox    & a_S,
                                             BaseFab<int> & a_R,
                                       const Box          & a_box)
{
  a_R.resize( a_box, 1 );

  FORT_DEFINE_REGIONS_3F( CHF_CONST_FRA(a_W),
                               CHF_FIA1(a_R,0),
                               CHF_BOX(a_box) );
}
                                                             // Fill ghost cells
void RiemannProblemMF::fillGhostCells(       FArrayBox&      a_W,
                                        const FArrayBox&      a_U,
                                        const int&            a_dir,
                                        const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int fluids = nFluids();
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexPrim(1);  

                                   // In periodic case, this doesn't do anything
  if( !m_domain.isPeriodic(a_dir) )
  {
    Box WBox = a_W.box();

                         // See if this chops off the high side of the input box
    Box tmp  = WBox;
    tmp     &= m_domain;
    

    int indW = WBox.bigEnd( a_dir );
    int indD =  tmp.bigEnd( a_dir );

    if( indW != indD )
    {
      int sign         = 1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indD + 1, indW - indD );

                                                         // Fill the ghost cells
      FORT_RIEMANNGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(boundaryBox) );
    }

    indW     = WBox.smallEnd( a_dir );
    indD     =  tmp.smallEnd( a_dir );

    if( indW != indD )
    {
      int sign         =-1;
      Box boundaryBox  = WBox;
      boundaryBox.setRange( a_dir, indW, indD - indW );

                                                         // Fill the ghost cells
      FORT_RIEMANNGS( CHF_FRA(a_W),
                      CHF_CONST_INT(sign),
                      CHF_CONST_INT(a_dir),
                      CHF_CONST_INT(iRhoN),
                      CHF_CONST_INT(fluids),
                      CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void RiemannProblemMF::getBCFlags( eBoundaryConditions leftBC,
                                    eBoundaryConditions rightBC,
                                    eBoundaryConditions bottomBC,
                                    eBoundaryConditions topBC,
                                    eBoundaryConditions frontBC,
                                    eBoundaryConditions behindBC )
{
  leftBC   = BC_Fixed;
  rightBC  = BC_Continuous;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Undefined;
  behindBC = BC_Undefined;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions RiemannProblemMF::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  switch( a_dir ){
  case 0  : return (a_sd == Side::Lo) ? BC_Fixed : BC_Continuous;
  case 1  : return BC_Periodic;
  default : return BC_Undefined;
  }
}
