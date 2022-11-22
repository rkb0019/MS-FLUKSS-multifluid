#include <iostream>
#include <iomanip>
#include <stdlib.h>

using std::ifstream;
using std::ios;


#include "LoHiSide.H"
#include "LoHiCenter.H"

#include "FreeStreamSpherical.H"
#include "FreeStreamF_F.H"
#include "LGintegrator.H"
#include "TecplotIO.H"
#include "MHDAMDefs.H"
#include "DebugOut.H"
#include "DebugF_F.H"
#include "RiemannSolver.H"
//#include "HeliosphericF_F.H"
#include "LoadBalance.H"

// Null constructor
FreeStreamSpherical::FreeStreamSpherical()
  : PhysProblem()
{
  m_isFortranCommonSet = false;  
  
}

// Input parameters
void FreeStreamSpherical::input( ParmParse & parser, int verbosity )
{  
  // Read problem specific parameters here
  
  parser.query( "gamma",   m_gamma );
  parser.query( "M",   m_M );
  parser.query( "aA",  m_aA );    
 
  m_verbosity = verbosity;

  // Print the parameters
  if( verbosity >= 2 )
  {
    pout() << "The FreeStreamSpherical input:" << endl;
    pout() << "gamma     = " << m_gamma     << endl;
    pout() << "M     = " << m_M     << endl;
    pout() << "aA    = " << m_aA    << endl;        
  }

  setFortranCommon(); 
      
}


// Sets parameters in a common block used by Fortran routines:
//   a_gamma          - Gamma for polytropic, gamma-law gas
void FreeStreamSpherical::setFortranCommon( )
{
  CH_assert(m_isFortranCommonSet == false);  
  
  FORT_SETFREESTREAM( CHF_CONST_REAL( m_gamma ),
                      CHF_CONST_REAL( m_M ),
                      CHF_CONST_REAL( m_aA )
                    );
                    
                         
  m_isFortranCommonSet = true;
}

// Set the flag m_isFortranCommonSet to true so that new IBCs made with
// new_PhysProblem() will have this flag set without calling setFortranCommon()
// (this is a clumsy design and should be improved).
void FreeStreamSpherical::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

/// Define the object
/**
   Set the problem domain index space and the grid spacing for this object.
 */
void FreeStreamSpherical::define(const ProblemDomain& a_domain,                                  
                                  const int            a_level)
{
  PhysProblem::define(a_domain, a_level);  
  
}


// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* FreeStreamSpherical::new_PhysProblem()
{
  FreeStreamSpherical* retval = new FreeStreamSpherical();
  
  retval->copy_PhysProblem(this);  
  
  return static_cast<PhysProblem*>(retval);
}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void FreeStreamSpherical::copy_PhysProblem(const PhysProblem* a_PP)
{
  const FreeStreamSpherical* PP = dynamic_cast<const FreeStreamSpherical*>(a_PP);
  if (PP == NULL) MayDay::Error("FreeStreamSpherical::copy_PhysProblem. Wrong argument");
  
  PhysProblem::copy_PhysProblem(a_PP);  
  
  if (PP->m_isFortranCommonSet == true)
  {
    this->m_gamma = PP->m_gamma;
    this->m_M     = PP->m_M;
    this->m_aA    = PP->m_aA;    
  
    this->setFortranCommonSet();
  }  
}
                                                             // Fill ghost cells
void FreeStreamSpherical::fillGhostCells(       FArrayBox&      a_W,
                                       const FArrayBox&      a_U,
                                       const int&            a_dir,
                                       const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  int indW,indD,i;
  
  Box WBox = a_W.box();

                       // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  
  
  int nGS;
  int nGSMax = 4; // How many two ghost cells must be filled
  
  
  // Outer boundary
  indW = WBox.bigEnd( a_dir );
  indD =  tmp.bigEnd( a_dir );    
  
  nGS  = MIN(abs(indW-indD),nGSMax);
      
  if( indW != indD )
  {
    int sign         = 1;
    Box boundaryBox  = WBox;
    boundaryBox.setRange( a_dir, indD + 1, nGS );
  
                                                       // Fill the ghost cells   
    if (CoordinateSystem == CoordinateSystemHandler::CS_Polar)
    FORT_FREESTREAMPOLARW(
        CHF_FRA(a_W),           
        CHF_BOX(boundaryBox),
        CHF_CONST_INT(m_level));
        
    if (CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {      
      FORT_FREESTREAMSPHERICALW(
          CHF_FRA(a_W),           
          CHF_BOX(boundaryBox),
          CHF_CONST_INT(m_level));
            
    }
    if (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian)
    {
      FORT_FREESTREAMCARTESIANW(
          CHF_FRA(a_W),           
          CHF_BOX(boundaryBox));
    }
  }

  // Inner boundary  
  indW = WBox.smallEnd( a_dir );
  indD =  tmp.smallEnd( a_dir );
  
  nGS  = MIN(abs(indW-indD),nGSMax);
  
  if( indW != indD )
  {
    int sign         =-1;
    Box boundaryBox  = WBox;
    boundaryBox.setRange( a_dir, indD-nGS , nGS );

                                                       // Fill the ghost cells      
    if (CoordinateSystem == CoordinateSystemHandler::CS_Polar)
    FORT_FREESTREAMPOLARW(
        CHF_FRA(a_W),           
        CHF_BOX(boundaryBox),
        CHF_CONST_INT(m_level));
        
    if (CoordinateSystem == CoordinateSystemHandler::CS_Spherical)
    {      
      FORT_FREESTREAMSPHERICALW(
          CHF_FRA(a_W),           
          CHF_BOX(boundaryBox),
          CHF_CONST_INT(m_level));
            
    }
    if (CoordinateSystem == CoordinateSystemHandler::CS_Cartesian)
    {
      FORT_FREESTREAMCARTESIANW(
          CHF_FRA(a_W),           
          CHF_BOX(boundaryBox));
    }
  }    
        
 
  #ifndef NDEBUG
    FORT_VIEWBOXDATA(
      CHF_FRA(a_W)
      );
  #endif
}


// Set boundary fluxes
void FreeStreamSpherical::fluxBC(    FArrayBox&      a_F,
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
      
      #ifndef NDEBUG
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WPlus)
        );
      FORT_VIEWBOXDATACONST(
        CHF_CONST_FRA(a_WMinus)
        );          
      #endif
           
      
      // boundaryBox is node centered. 
      // We need cell centered box that lies outside the problem domain.
      Box mainBox = boundaryBox;
      mainBox.shiftHalf(a_dir,sign);                
      CH_assert(mainBox.type(a_dir) == IndexType::CELL);              
                      
      CH_assert(m_RS!=NULL);              
                  
      m_RS->fluxes( a_F, a_WPlus, a_WMinus,  a_dir, WRHO, boundaryBox );
        
     
                                                             
    }
  }
}

// Set up initial conditions
void FreeStreamSpherical::initialize(LevelData<FArrayBox>& a_U)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();

  // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
    // Storage for current grid
    FArrayBox& U = a_U[dit()];

    // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain;
    
    FORT_FREESTREAMSPHERICALU(
          CHF_FRA(U),           
          CHF_BOX(uBox));  
  }
  
}
                      

//                            Return boundary condition flags for all boundaries
void FreeStreamSpherical::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions FreeStreamSpherical::getBCFlags( int a_dir, Side::LoHiSide a_sd )
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
void FreeStreamSpherical::tagCells(const FArrayBox&  a_U,
                                    const Box&        a_box,
                                          IntVectSet& a_tags)
{
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();
  
  if (CoordinateSystem == CoordinateSystemHandler::CS_Polar)     
  {  
    Real PI = 3.141592654;
    Real dphi = m_csh->dx(1,m_level);
    
    IntVect lo,hi; Box b1,bc;
  
    lo[0] = (int)(m_domain.domainBox().size(0)/3);
    lo[1] = (int)(  ((70.0/180.0)*PI)/dphi );
    
    hi[0] = (int)(m_domain.domainBox().size(0)/2);
    hi[1] = (int)( ((110.0/180.0)*PI)/dphi );
    
    lo[0] = hi[0];
    lo[1] = hi[1];
    
    b1.define(lo,hi);
    bc = b1&a_box;

    if (!bc.isEmpty())  
    {
      // Tag all bc
      BoxIterator bit(bc);
      for (bit.begin(); bit.ok(); ++bit)
      {
        const IntVect& iv = bit();                  
        a_tags |= iv;      
      }
    }
  }
  else
  {
    IntVect iv(D_DECL(
      m_domain.domainBox().size(0)/2,
      m_domain.domainBox().size(1)/2,
      m_domain.domainBox().size(2)/2));
    
    iv += m_domain.domainBox().smallEnd();
    a_tags |= iv;      
    
    //a_tags.grow(4);
  }
}                                         

                                              // Problem specific postprocessing
void FreeStreamSpherical::postprocessing(       FArrayBox & a_U,
                                    const FArrayBox & a_W,
                                    const Real      & a_dt,
                                    const Real      & a_time,                                    
                                    const Box       & a_box       )
{
  return;
  Box updBox = m_domain.domainBox();
  updBox.setSmall(2,0);
  updBox.setBig(2,0);  
  updBox &= a_box;
  
  FORT_FREESTREAMSPHERICALU(
          CHF_FRA(a_U),           
          CHF_BOX(updBox));  
          
  updBox = m_domain.domainBox();
  updBox.setSmall(2,m_domain.domainBox().bigEnd(2));
  updBox.setBig(2,m_domain.domainBox().bigEnd(2));  
  updBox &= a_box;
  
  FORT_FREESTREAMSPHERICALU(
          CHF_FRA(a_U),           
          CHF_BOX(updBox));  
                                     
}

