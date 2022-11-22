
#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "ShearFlowProblem.H"
#include "ShearFlowF_F.H"
#include "RiemannSolver.H"
#include "LGintegrator.H"

// Null constructor
ShearFlowProblem::ShearFlowProblem()
{
    m_isFortranCommonSet = false;
}


// Input parameters
void ShearFlowProblem::input( ParmParse & parser, int verbosity )
{
  m_N1  = 7;
  m_V1  = 450e+5;
  m_T1  = 51100;
  
  m_N2  = 7;
  m_V2  = 450e+5;
  m_T2  = 51100;
  Real fs_a = -1;
  
  m_physModel = PP_EulerPM;
    
  m_gamma      = 1.6666666667;

  parser.query( "gamma",     m_gamma     );

  parser.query( "N1",  m_N1  );
  parser.query( "V1",  m_V1  );
  parser.query( "T1",  m_T1  );
  
  parser.query( "N2",  m_N2  );
  parser.query( "V2",  m_V2  );
  parser.query( "T2",  m_T2  );
  
  parser.query( "fs_a",  fs_a);
  
  int physModel;
  
  if (parser.contains("phys_model"))
  {
    parser.query( "phys_model",  physModel);    
    if ((physModel == PP_EulerPM) || (physModel == PP_MHDPM))
      m_physModel = (ePhysicalModel)(physModel);
    else
      MayDay::Warning("Unsuported phys_model");    
  }
  
  

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "Shear flow problem input:"      << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "N1  = " << m_N1 << endl;
    pout() << "V1  = " << m_V1 << endl;
    pout() << "T1  = " << m_T1 << endl;
    pout() << "N2  = " << m_N2 << endl;
    pout() << "V2  = " << m_V2 << endl;
    pout() << "T2  = " << m_T2 << endl;    
  }

  setFortranCommon( m_gamma,
                    m_N1,m_V1,m_T1,m_N2,m_V2,m_T2, fs_a);
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void ShearFlowProblem::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}


// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void ShearFlowProblem::setFortranCommon( const Real&     a_gamma,
                                     const Real&     a_N1,
                                     const Real&     a_V1,
                                     const Real&     a_T1,
                                     const Real&     a_N2,
                                     const Real&     a_V2,
                                     const Real&     a_T2,
                                     const Real&     fs_a)
{
    CH_assert(m_isFortranCommonSet == false);

    FORT_SETSHEARFLOW( CHF_CONST_REAL( a_gamma    ),
                   CHF_CONST_REAL( a_N1  ),
                   CHF_CONST_REAL( a_V1  ),
                   CHF_CONST_REAL( a_T1  ),
                   CHF_CONST_REAL( a_N2  ),
                   CHF_CONST_REAL( a_V2  ),
                   CHF_CONST_REAL( a_T2  ),                   
                   CHF_CONST_REAL( fs_a  ) );

    m_isFortranCommonSet = true;
}


// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* ShearFlowProblem::new_PhysProblem()
{
  ShearFlowProblem* retval = new ShearFlowProblem();
  
  retval->copy_PhysProblem(this);

  if( m_isFortranCommonSet == true )
  {
    retval->m_gamma      = this->m_gamma;
    retval->m_N1   = this->m_N1;
    retval->m_V1   = this->m_V1;
    retval->m_T1   = this->m_T1;
    retval->m_N2   = this->m_N2;
    retval->m_V2   = this->m_V2;
    retval->m_T2   = this->m_T2;
    
    retval->setFortranCommonSet();
  }  

  return static_cast<PhysProblem*>(retval);
}

// Set up initial conditions
void ShearFlowProblem::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
    
  DataIterator dit = a_U.boxLayout().dataIterator();
  
  int jmiddle = m_domain.size(1)/2;
  Real dx = m_csh->dx(0,m_level);

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    
    // Set up initial condition in this grid
    FORT_SHEARFLOWINIT( CHF_CONST_FRA(U),
                    CHF_CONST_INT(jmiddle),
                    CHF_CONST_REAL(dx),
                    CHF_BOX(uBox));
  }
}


// Set boundary fluxes
void ShearFlowProblem::fluxBC(       FArrayBox&      a_F,
                                 FArrayBox&      a_Bn,
                           const FArrayBox&      a_WMinus,
                           const FArrayBox&      a_WPlus,
                           const int&            a_dir,
                           const Side::LoHiSide& a_side,
                           const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

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
       
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
     
    }
  }
}

                                                             // Fill ghost cells
void ShearFlowProblem::fillGhostCells(       FArrayBox&      a_W,
                                   const FArrayBox&      a_U,
                                   const int&            a_dir,
                                   const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  int jmiddle = m_domain.size(1)/2;

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
      FORT_SHEARFLOWGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(jmiddle),
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
      FORT_SHEARFLOWGS( CHF_FRA(a_W),
                    CHF_CONST_INT(sign),
                    CHF_CONST_INT(a_dir),
                    CHF_CONST_INT(jmiddle),
                    CHF_BOX(boundaryBox) );
    }
  }
}

//                            Return boundary condition flags for all boundaries
void ShearFlowProblem::getBCFlags( eBoundaryConditions leftBC,
                               eBoundaryConditions rightBC,
                               eBoundaryConditions bottomBC,
                               eBoundaryConditions topBC,
                               eBoundaryConditions frontBC,
                               eBoundaryConditions behindBC )
{
  leftBC   = BC_Continuous;
  rightBC  = BC_Fixed;
  bottomBC = BC_Periodic;
  topBC    = BC_Periodic;
  frontBC  = BC_Periodic;
  behindBC = BC_Periodic;
}

//                         Return the boundary condition flag for given boundary
PhysProblem::eBoundaryConditions ShearFlowProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    
    return (a_sd == Side::Lo ) ? BC_Continuous : BC_Fixed;
  } else {
    return BC_Periodic;
  }
}
