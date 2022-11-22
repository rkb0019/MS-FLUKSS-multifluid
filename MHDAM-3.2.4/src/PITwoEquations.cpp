#include "PITwoEquations.H"
#include "PatchMHDAMF_F.H"
#include "PickupIonsF_F.H"
#include "LGintegrator.H"
#include "parstream.H"

//                                                                   Constructor
PITwoEquations::PITwoEquations( void )
{
  m_consStateInterval.define( 0, 1 );
  m_primStateInterval.define( 0, 1 );

  m_dGamma     = 1.6666666666666666667;
  m_dSmallVal  = 1.0e-12;

  m_isFortranCommonSet = false;
}

//                                                                    Destructor
PITwoEquations::~PITwoEquations( void )
{
}

//                                                              Input parameters
void PITwoEquations::input( ParmParse & parser, int a_verbosity )
{
  m_verbosity  = a_verbosity;

  m_dGamma     = 1.6666666666666666667;
  m_dSmallVal  = 1.0e-12;

  parser.query( "gammaPI",   m_dGamma    );
  parser.query( "smallPI",   m_dSmallVal );

                                                         // Print the parameters
  if( m_verbosity >= 2 )
  {
    pout() << "Pickup ions parametrs input:"     << endl;
    pout() << "gammaPI   = " << m_dGamma      << endl;
    pout() << "smallPI   = " << m_dSmallVal  << endl;
  }

  setFortranCommon( m_dSmallVal,  m_dGamma );
}

//                    Sets parameters in a common block used by Fortran routines
void PITwoEquations::setFortranCommon( const Real& a_dSmallVal,
                                       const Real& a_dGamma     )
{
  CH_assert(m_isFortranCommonSet == false);

  FORT_SETCONST_PI( CHF_CONST_REAL( a_dGamma  ),
                    CHF_CONST_REAL( a_dSmallVal ) );

  m_isFortranCommonSet = true;
}

//                                     Set the flag m_isFortranCommonSet to true
void PITwoEquations::setFortranCommonSet()
{
  m_isFortranCommonSet = true;
}

//                               Factory method - this object is its own factory
PickupIons * PITwoEquations::new_PickupIons( void )
{
  PITwoEquations * retval = new PITwoEquations();

  return static_cast<PickupIons*>(retval);
}

//                                  Copy internal data to another PITwoEquations
void PITwoEquations::copyTo( PickupIons * pPIons ) const
{
  PickupIons::copyTo( pPIons );

  PITwoEquations* pPI  = dynamic_cast<PITwoEquations*>(pPIons);

  pPI->m_dGamma    = m_dGamma;
  pPI->m_dSmallVal = m_dSmallVal;
}

//                                 Generate names for the conservative variables
void PITwoEquations::AddConsNames( Vector<string> * pNames )
{
  pNames->push_back( "density_PI" );
  pNames->push_back( "pressure_PI" );
}

//                                    Generate names for the primitive variables
void PITwoEquations::AddPrimNames( Vector<string> * pNames )
{
  pNames->push_back( "density_PI" );
  pNames->push_back( "pressure_PI" );
}

//     Compute the primitive variables from the conserved variables within a_box
void PITwoEquations::stateToPrim(       FArrayBox & a_W,
                                  const FArrayBox & a_U,
                                  const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  a_W.copy( a_U, a_box, iRhoPIU, a_box, iRhoPIW, 2 );
}

//     Compute the conserved variables from the primitive variables within a_box
void PITwoEquations::primToState(       FArrayBox & a_U,
                                  const FArrayBox & a_W,
                                  const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  a_U.copy( a_W, a_box, iRhoPIW, a_box, iRhoPIU, 2 );
}

//                                     Calculate fluxes of pickup ions variables
//                                     this is old for ugradp
void PITwoEquations::primToFlux(       FArrayBox & a_F,
                                 const FArrayBox & a_W,
                                 const int &       a_dir,
                                 const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  FORT_FLUXE_PI( CHF_FRA(a_F),
                 CHF_CONST_FRA(a_W),
                 CHF_CONST_INT(iRhoPIU),
                 CHF_CONST_INT(iRhoPIW),
                 CHF_CONST_INT(a_dir),
                 CHF_BOX(a_box) );
}

//this is new for pdivu
void PITwoEquations::primToFlux(       FArrayBox & a_F,
                                 const FArrayBox & a_W,
                            const BaseFab<int>   & a_REG,
                                 const int &       a_dir,
                                 const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  FORT_FLUXE_PI_REG( CHF_FRA(a_F),
                     CHF_CONST_FRA(a_W),
                     CHF_CONST_FIA1(a_REG,0),
                     CHF_CONST_INT(iRhoPIU),
                     CHF_CONST_INT(iRhoPIW),
                     CHF_CONST_INT(a_dir),
                     CHF_BOX(a_box) );
}

//               Compute an upwind fluxes at the faces for pickup ions variables
//               this is old for ugradp
void PITwoEquations::upwindFluxes(       FArrayBox & a_F,
                                   const FArrayBox & a_WPlus,
                                   const FArrayBox & a_WMinus,
                                   const int &       a_dir,
                                   const int &       a_iRho,
                                   const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
  FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

                                   // Riemann solver computes a_F for all edges
                                      // that are not on the physical boundary.
  FORT_UPWINDSCALARFLUXES_PI( CHF_FRA(a_F),
                              CHF_CONST_FRA(shiftWLeft),
                              CHF_CONST_FRA(shiftWRight),
                              CHF_CONST_INT(a_dir),
                              CHF_CONST_INT(a_iRho),
                              CHF_CONST_INT(iRhoPIU),
                              CHF_CONST_INT(iRhoPIW),
                              CHF_BOX(a_box) );

                            // Shift the left and right primitive variable boxes
                            // back to their original position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}

//new fluxes with pdivu
void PITwoEquations::upwindFluxes(       FArrayBox & a_F,
                                   const FArrayBox & a_WPlus,
                                   const FArrayBox & a_WMinus,
                                const BaseFab<int> & a_REG,
                                   const int &       a_dir,
                                   const int &       a_iRho,
                                   const Box &       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

         // Cast away "const" inputs so their boxes can be shifted left or right
         // 1/2 cell and then back again (no net change is made!)
  FArrayBox& shiftWLeft  = (FArrayBox&)a_WPlus;
  FArrayBox& shiftWRight = (FArrayBox&)a_WMinus;

                            // Shift the left and right primitive variable boxes
                            // 1/2 cell so they are face centered
  shiftWLeft .shiftHalf(a_dir, 1);
  shiftWRight.shiftHalf(a_dir,-1);

                                   // Riemann solver computes a_F for all edges
                                      // that are not on the physical boundary.
  FORT_UPWINDSCALARFLUXES_PI_REG( CHF_FRA(a_F),
                                  CHF_CONST_FRA(shiftWLeft),
                                  CHF_CONST_FRA(shiftWRight),
                                  CHF_CONST_FIA1(a_REG,0),
                                  CHF_CONST_INT(a_dir),
                                  CHF_CONST_INT(a_iRho),
                                  CHF_CONST_INT(iRhoPIU),
                                  CHF_CONST_INT(iRhoPIW),
                                  CHF_BOX(a_box) );
/*  FORT_UPWINDSCALARFLUXES_PI( CHF_FRA(a_F),
                              CHF_CONST_FRA(shiftWLeft),
                              CHF_CONST_FRA(shiftWRight),
                              CHF_CONST_INT(a_dir),
                              CHF_CONST_INT(a_iRho),
                              CHF_CONST_INT(iRhoPIU),
                              CHF_CONST_INT(iRhoPIW),
                              CHF_BOX(a_box) );*/

                            // Shift the left and right primitive variable boxes
                            // back to their original position
  shiftWLeft .shiftHalf(a_dir,-1);
  shiftWRight.shiftHalf(a_dir, 1);
}
//                                   Calculate pickup ions specific source terms
//                                   this is old for ugradp
void PITwoEquations::explicitSource(       FArrayBox & a_U,
                                           FArrayBox & a_S,
                                     const FArrayBox & a_W,
                                     const Real      & a_dt,
                                     const int       & a_level,
                                     const Box       & a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  pout() << "OLD_ExplicitSource" <<  endl;//FF
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();//m_csh included in PickupIons.H
  if ( CoordinateSystem == CoordinateSystemHandler::CS_Cartesian)
  {
       MayDay::Error("ERROR: For Cartesian grid, the old form (U gradP) of pressure equation for PUIs is not implemented. SOURCE_PICKUP_IONS_CARTESIAN not present");
  }



  FORT_SOURCE_PICKUP_IONS( CHF_FRA(a_S),
                           CHF_CONST_FRA(a_W),
                           CHF_CONST_REAL(a_dt),
                           CHF_CONST_INT(iRhoPIU),
                           CHF_CONST_INT(iRhoPIW),
                           CHF_CONST_INT(a_level),
                           CHF_BOX(a_box));

  int iPressPIU  = iRhoPIU + 1;

  FORT_ADDSOURCES( CHF_FRA(a_U),
                   CHF_CONST_FRA(a_S),
                   CHF_CONST_INT(iPressPIU),
                   CHF_CONST_INT(iPressPIU),
                   CHF_BOX(a_box) );
}

// Calculate pickup ion source terms with divU
void PITwoEquations::explicitSource(       FArrayBox & a_U,
                                           FArrayBox & a_S,
                                     const FArrayBox & a_W,
                                  const BaseFab<int> & a_REG,
                                     const Real      & a_dt,
                                     const FArrayBox & a_divU,
                                     const int       & a_level,
                                     const Box       & a_box )
{
  int iRhoPIU = m_consStateInterval.begin();
  int iRhoPIW = m_primStateInterval.begin();

  FORT_SOURCE_PICKUP_IONS_REG( CHF_FRA(a_S),
                               CHF_CONST_FRA(a_W),
                               CHF_CONST_FIA1(a_REG,0),
                               CHF_CONST_FRA1(a_divU,0),
                               CHF_CONST_REAL(a_dt),
                               CHF_CONST_INT(iRhoPIU),
                               CHF_CONST_INT(iRhoPIW),
                               CHF_CONST_INT(a_level),
                               CHF_BOX(a_box));

  int iPressPIU  = iRhoPIU + 1;

  FORT_ADDSOURCES( CHF_FRA(a_U),
                   CHF_CONST_FRA(a_S),
                   CHF_CONST_INT(iPressPIU),
                   CHF_CONST_INT(iPressPIU),
                   CHF_BOX(a_box) );
}

//                      Compute pickup ion model specific dt in all patch points
Real PITwoEquations::computeDt( const FArrayBox& a_U,
                                      FArrayBox& a_dt,
                                const int      & a_level,
                                const Box&       a_box )
{
  int iRhoPIU = m_consStateInterval.begin();

  FORT_MINDT_SPHERICAL_PI( CHF_FRA1(a_dt,0),
                           CHF_CONST_FRA(a_U),
                           CHF_CONST_INT(a_level),
                           CHF_CONST_INT(iRhoPIU),
                           CHF_BOX(a_box) );
}
