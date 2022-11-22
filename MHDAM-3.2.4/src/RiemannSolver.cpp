#include "RiemannSolver.H"
#include "AUSMPlusRS.H"
#include "Roe8WavesRS.H"
#include "Roe7WavesRS.H"
#include "RusanovRS.H"
#include "VanLeerRS.H"
#include "VFRoeRS.H"
#include "RoeHHTRS.H"
#include "RiemannHD.H"
#include "DednerRS.H"

// Static method for construction new objects
RiemannSolver * RiemannSolver::make( int                         iFlux,
                                     PhysProblem::ePhysicalModel eModel,
                                     Real                        dEntFix )
{
  RiemannSolver * pRS = NULL;

  if( eModel == PhysProblem::PP_EulerPM )
  {
    switch( iFlux ){
    case 2 : {
      pRS = new AUSMPlusRS();
      break;
    }
    case 3 : {
      pRS = new RiemannHD();
      break;
    }
    case 4 : {
      pRS = new RoeHHTRS();
      break;
    }
    default :
      pRS = new VanLeerRS();
      break;
    }
  }
  else
  {
    switch( iFlux ){
    case 1 : {                   // "Roe" flux with density arithmetic averaging
      Roe8WavesRS * pRoe = new Roe8WavesRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 0, 0 );
      pRS  = pRoe;
      break;
    }
    case 2 : {
                                       // Rusanov flux with arithmetic averaging
      RusanovRS * pRUS = new RusanovRS();
      pRUS->setParameters( 1.0e-9, 0 );
      pRS  = pRUS;
      break;
    }
    case 3 : {
                                   // "Roe" flux with density weighted averaging
      Roe8WavesRS * pRoe = new Roe8WavesRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 1, 0 );
      pRS  = pRoe;
      break;
    }
    case 4 : {
                                 // Rusanov flux with density weighted averaging
      RusanovRS * pRUS = new RusanovRS();
      pRUS->setParameters( 1.0e-9, 1 );
      pRS  = pRUS;
      break;
    }
    case 5 : {
                                         // VFRoe flux with arithmetic averaging
      VFRoeRS * pRoe = new VFRoeRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 0 );
      pRS  = pRoe;
      break;
    }    
    case 11 : {                   // "Roe" flux with density arithmetic averaging
      Roe7WavesRS * pRoe = new Roe7WavesRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 0, 0 );
      pRS  = pRoe;
      break;
    }    
    case 13 : {
                                   // "Roe" flux with density weighted averaging
      Roe7WavesRS * pRoe = new Roe7WavesRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 1, 0 );
      pRS  = pRoe;
      break;
    }        
    case 101 : {
                         // "Roe" flux with density weighted averaging  + Dedner
      DednerRS * pDRS = new DednerRS();
      pDRS->setParameters( dEntFix, 1.0e-9, 0, 0 );
      pRS  = pDRS;
      break;
    }
    case 103 : {
                         // "Roe" flux with density weighted averaging  + Dedner
      DednerRS * pDRS = new DednerRS();
      pDRS->setParameters( dEntFix, 1.0e-9, 1, 0 );
      pRS  = pDRS;
      break;
    }
    default :
      Roe8WavesRS * pRoe = new Roe8WavesRS();
      pRoe->setParameters( dEntFix, 1.0e-9, 0, 0 );
      pRS  = pRoe;
      break;
    }
  }

  return pRS;
}
