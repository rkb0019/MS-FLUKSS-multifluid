#include <time.h>

#include "BoxIterator.H"
#include "IntVectSet.H"

#include "HelioRealBCProblem.H"
#include "HeliosphericF_F.H"

#include "MHDAMDefs.H"

#include "LGintegrator.H"
#include "DebugF_F.H"
#include "EqSysMHDMF.H"
#include "HelioRealBCProblemF_F.H"
#include "EosCommon.H"
#include "HeliosphericPlasmaBCF_F.H"
#include "TecplotIO.H"
#include "SourceCalculator.H"

#include "TMBreechEtAl2008.H"
#include "HeliosphericTurbF_F.H"
#include "PITwoEquations.H"
#include "HeliosphericPIF_F.H"


// Null constructor
HelioRealBCProblem::HelioRealBCProblem()
 : HeliosphericProblem()
{  
  
}

HelioRealBCProblem::~HelioRealBCProblem()
{  

}

void HelioRealBCProblem::input( ParmParse & parser, int verbosity )
{
  HeliosphericProblem::input( parser, verbosity );
  
  m_startBC = 0.0;
  m_startBCPhys = 2005;
  m_iOldOMNIFile = 0;
  parser.query( "startGMIR", m_startBC   );
  parser.query( "startBC",   m_startBC   );
  parser.query( "startBCPhys", m_startBCPhys);    
  
  if (m_subproblem == HPBC_OMNIPROBLEM)   
  {
      parser.query("OMNIData", m_OMNIDataFile);    
      parser.query("iOldOMNIFile", m_iOldOMNIFile);
  }
    
  if (m_subproblem == HPBC_WSOPROBLEM)   
  {
    parser.query("startBCPhys", m_startBCPhys);    
  }
  
  if (m_subproblem == HPBC_CIRPROBLEM)   
  {
    parser.query("CIRfile", m_h5File);           
    if (parser.contains( "h5input" ))
      parser.query("h5input", m_h5File);           
    
  }
  
  if (m_subproblem == HPBC_HDF5INPUT)   
  {
    if (parser.contains("h5input"))
      parser.query("h5input", m_h5File);     
    else
      MayDay::Error("\"h5input\" parameter is not found in the input file");
      
    int inputRotating = 1;
    parser.query("h5input_rotating_frame", inputRotating);        
    m_bDataRotatingFrame = (inputRotating == 1);
  }
  
  if (parser.contains("Vprobe") || parser.contains("Uprobe") || parser.contains("earth_probe"))
    MayDay::Error("Vprobe, Uprobe or earth_probe parameters are not valid any more");
        
}

// Factory method - this object is its own factory:
//   Return a pointer to a new PhysProblem object with m_isDefined = false (i.e.,
//   its define() must be called before it is used) and m_isFortranCommonSet
//   set to value of m_isFortranCommonset in the current (factory) object.
PhysProblem* HelioRealBCProblem::new_PhysProblem()
{
  HelioRealBCProblem* retval = new HelioRealBCProblem();
  
  retval->copy_PhysProblem(this);  

  return static_cast<PhysProblem*>(retval);

}

/// Copy method 
//     Copy all data from a_PP to this instance.   
void HelioRealBCProblem::copy_PhysProblem(const PhysProblem* a_PP)
{
  const HelioRealBCProblem* PP = dynamic_cast<const HelioRealBCProblem*>(a_PP);
  
  if (PP == NULL) MayDay::Error("HelioRealBCProblem::copy_PhysProblem. Wrong argument");
  
  HeliosphericProblem::copy_PhysProblem(a_PP);
  
  this->m_startBC  = PP->m_startBC;  
  this->m_OMNIDataFile = PP->m_OMNIDataFile;
  this->m_startBCPhys = PP->m_startBCPhys;  
  this->m_h5File = PP->m_h5File;
  this->m_bDataRotatingFrame = PP->m_bDataRotatingFrame;
    
}

void HelioRealBCProblem::define(const ProblemDomain& a_domain,                                  
                                const int            a_level)
{
  HeliosphericProblem::define(a_domain,a_level);
  
  if (m_subproblem == HPBC_GMIRPROBLEM)
    readGMIRData();
  if (m_subproblem == HPBC_SUESSPROBLEM)
    readThetaDistrib();
  if( m_subproblem == HPBC_OMNIPROBLEM )
  {
    TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
    if( pTurbMod != NULL )
    {
      readOMNIDataTurb();
    } else {
      readOMNIData();
    }
  }
  if (m_subproblem == HPBC_WSOPROBLEM)
    readWSOData();
  if (m_subproblem == HPBC_CIRPROBLEM)
    readCIRData();
  if (m_subproblem == HPBC_HDF5INPUT)
    prepareHDF5InputProblem();
            
}

void HelioRealBCProblem::defineMesh(const ProblemDomain & a_prob_domain,
                                     const Vector<Real>  & a_domainBox)
{  
  if (m_csh->coordinateSystem() != CoordinateSystemHandler::CS_Spherical) return;    
  
  HeliosphericProblem::defineMesh(a_prob_domain,a_domainBox); return;      
      
}

void HelioRealBCProblem::fillGhostCells(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  if (a_dir == 0)
    HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
  
  
  if (m_subproblem == HPBC_GMIRPROBLEM)
    fillGhostCellsGMIR(a_W, a_U, a_dir, a_time);      
  if (m_subproblem == HPBC_SUESSPROBLEM)
    fillGhostCellsSuess(a_W, a_U, a_dir, a_time);      
  if (m_subproblem == HPBC_OMNIPROBLEM)
    fillGhostCellsOMNI(a_W, a_U, a_dir, a_time);      
  if (m_subproblem == HPBC_WSOPROBLEM)
    fillGhostCellsWSO(a_W, a_U, a_dir, a_time);      
  if (m_subproblem == HPBC_CIRPROBLEM)
    fillGhostCellsCIR(a_W, a_U, a_dir, a_time);      
  if (m_subproblem == HPBC_HDF5INPUT)
    fillGhostCellsH5Input(a_W, a_U, a_dir, a_time);      
}


                                                             // Fill ghost cells
void HelioRealBCProblem::fillGhostCellsOdstrcil(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
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
      
      if (a_dir == 0)      
      FORT_ODSTRCILGSSPHERICAL(
          CHF_FRA(a_W),                     
          CHF_CONST_INT(sign),
          CHF_CONST_FRA(m_OdstrcilData),          
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
            
      if (a_dir == 2)
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
      
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
      FORT_ODSTRCILGSSPHERICAL(
          CHF_FRA(a_W),                     
          CHF_CONST_INT(sign),
          CHF_CONST_FRA(m_OdstrcilData),          
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
      
      if (a_dir == 2)
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }    

  }
  
  
}

void HelioRealBCProblem::fillGhostCellsGMIR(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
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
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
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
      
      if (a_dir == 0)      
      FORT_GMIRGSSPHERICAL(
          CHF_FRA(a_W),                     
          CHF_CONST_FRA(a_U),    
          CHF_CONST_REAL(m_startBC),                 
          CHF_CONST_R1D(m_GMIRData.time,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.N,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.V,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.T,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.B,m_GMIRData.nPoints),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(sign),
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
            
      if (a_dir == 2)
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
      
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
      FORT_GMIRGSSPHERICAL(
          CHF_FRA(a_W),                     
          CHF_CONST_FRA(a_U),   
          CHF_CONST_REAL(m_startBC),
          CHF_CONST_R1D(m_GMIRData.time,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.N,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.V,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.T,m_GMIRData.nPoints),
          CHF_CONST_R1D(m_GMIRData.B,m_GMIRData.nPoints),
          CHF_CONST_INT(iRhoN),
          CHF_CONST_INT(fluids),
          CHF_CONST_INT(sign),
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
      
      if (a_dir == 2)
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }    

  }
  
  
}

void HelioRealBCProblem::fillGhostCellsSuess(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iHCSb  = (m_eqSys->numTrackingSurfaces()>1 ? iHCS+1 : -1);

    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  

  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain;  
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = abs(indW-indD);    
           
        
    if( indW != indD )
    {
      /*int sign         = 1;
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD + 1, nGS );
            
      Box srcBox = WBox;
      Box desBox = WBox;
      int i_srcBox = indD;
      int i_desBox = indD + 1;
      
      if (a_dir == 0)      
      {      
        for (int i = 0; i < boundaryBox.size(a_dir); i++) 
        {
          srcBox.setRange( a_dir, i_srcBox, 1 );
          desBox.setRange( a_dir, i_desBox, 1 );
          a_W.copy(a_W, srcBox, 0, desBox, 0, a_W.nComp());        
          
          i_srcBox--;
          i_desBox++;
        }
      }*/

      
      
     
      
            
      if (a_dir == 2)      
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
            
      
    }
       
    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if ( indW != indD )
    {
      
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
      
      if (a_dir == 0)
      {
      
        Real Tref   = eos_AU/(m_lismV*3600.0*24.0*365.25);
        Real timeCycle   = (a_time-m_startBC)*Tref; // years
        
        if (timeCycle<0.0) timeCycle = 0.0;
        //if (timeCycle>m_ThetaDistrib[m_nThetaDistrib-1].time) 
        //  timeCycle = m_ThetaDistrib[m_nThetaDistrib-2].time+0.5*(m_ThetaDistrib[m_nThetaDistrib-1].time-m_ThetaDistrib[m_nThetaDistrib-2].time);
        
        Real tilt = m_sunTILT;
        if (m_WSOData.time!=NULL)
        {
          Real * tilt_iter = std::lower_bound(m_WSOData.time, m_WSOData.time + m_WSOData.nPoints, timeCycle);
          if (tilt_iter!=m_WSOData.time + m_WSOData.nPoints)
          {
            int tilt_ind = tilt_iter - m_WSOData.time;
            CH_assert((tilt_ind>=0)&&(tilt_ind<m_WSOData.nPoints));
            if (tilt_ind == m_WSOData.nPoints - 1)         
              tilt = m_WSOData.Tilt[tilt_ind];
            else {
              Real time1 = m_WSOData.time[tilt_ind];
              Real time2 = m_WSOData.time[tilt_ind+1];
              Real tilt1 = m_WSOData.Tilt[tilt_ind];
              Real tilt2 = m_WSOData.Tilt[tilt_ind+1];
              tilt = (tilt2-tilt1)/(time2-time1)*(timeCycle-time1)+tilt1;
            }
          }
        }
        
        if (procID() == 0) pout() << "tilt: " << tilt << endl;
        
        int timeind = (int)(timeCycle/(m_ThetaDistrib[1].time-m_ThetaDistrib[0].time));        
        
        int kdomain = m_domain.domainBox().size(2);
        
        if (procID() == 0) pout() << "ind = " << timeind << ", time = " << timeCycle*365.25 << " " << "days " << endl;
        
        if ((timeind>=0) && ((timeind-2)<m_nThetaDistrib))
        {   
          if (! ((m_ThetaDistrib[timeind].time<=timeCycle)&&(timeCycle<=m_ThetaDistrib[timeind+1].time)))
          pout() << 
              "CH_ASSERT: m_ThetaDistrib[" << timeind << "].time = " << m_ThetaDistrib[timeind].time << 
              " <=  " << timeCycle << " <=  " << 
              "m_ThetaDistrib[" << timeind+1 << "].time = " << m_ThetaDistrib[timeind+1].time <<  endl;
            
          CH_assert((m_ThetaDistrib[timeind].time<=timeCycle)&&(timeCycle<=m_ThetaDistrib[timeind+1].time));
          
          FORT_SUESSTHETABC(
              CHF_FRA(a_W),                    
              CHF_CONST_R1D(m_ThetaDistrib[timeind].N,kdomain),
              CHF_CONST_R1D(m_ThetaDistrib[timeind].V,kdomain),
              CHF_CONST_R1D(m_ThetaDistrib[timeind].T,kdomain),     
              CHF_CONST_REAL(m_ThetaDistrib[timeind].time),     
              CHF_CONST_R1D(m_ThetaDistrib[timeind+1].N,kdomain),
              CHF_CONST_R1D(m_ThetaDistrib[timeind+1].V,kdomain),
              CHF_CONST_R1D(m_ThetaDistrib[timeind+1].T,kdomain),     
              CHF_CONST_REAL(m_ThetaDistrib[timeind+1].time),     
              CHF_CONST_REAL(timeCycle),
              CHF_CONST_REAL(a_time),
              CHF_CONST_REAL(tilt),
              CHF_CONST_INT(iHCS),          
              CHF_CONST_INT(iHCSb),          
              CHF_CONST_INT(m_level),    
              CHF_BOX(boundaryBox));
        } else
        {
          MayDay::Abort("Stop cycle");
          FORT_HELIOGSPLASMASPHDEFAULT(
           CHF_FRA(a_W),     
           CHF_CONST_INT(iHCS),
           CHF_CONST_INT(iHCS),
           CHF_CONST_INT(m_level),
           CHF_CONST_REAL(a_time),
           CHF_BOX(boundaryBox));
        }
      }
      if (a_dir == 2)  
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }
  
  }
}


void HelioRealBCProblem::fillGhostCellsOMNI(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;

  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  
  int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iGMIR  = (m_eqSys->numTrackingSurfaces()>1 ? m_eqSys->lvlsStateInterval().begin()+1 : -1);

  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
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
            
      Box srcBox = WBox;
      Box desBox = WBox;
      int i_srcBox = indD;
      int i_desBox = indD + 1;
      
      srcBox.setRange( a_dir, i_srcBox, 1 );
      
      if (a_dir == 0)      
      {      
        for (int i = 0; i < boundaryBox.size(a_dir); i++) 
        {          
          desBox.setRange( a_dir, i_desBox, 1 );
          a_W.copy(a_W, srcBox, 0, desBox, 0, WNUM);        
                    
          i_desBox++;
        }

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
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);
      
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
  
      if( a_dir == 0 )
      {
        FORT_OMNIGSSPHERICAL(
          CHF_FRA(a_W),
          CHF_CONST_REAL(m_startBC),
          CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.N,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.V,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.T,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Br,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Bp,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Bt,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.B,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Tilt,m_OMNIData.nPoints),
          CHF_CONST_INT(iHCS),
          CHF_CONST_INT(iGMIR),
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
              FORT_OMNIGSSPHERICAL_TM( CHF_FRA(a_W),
                                       CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
                                       CHF_CONST_R1D(m_OMNIData.Z2,m_OMNIData.nPoints),
                                       CHF_CONST_INT(iZ2),
                                       CHF_CONST_REAL(a_time),
                                       CHF_CONST_REAL(m_startBC),
                                       CHF_BOX(boundaryBox) );
            }
          }

          PickupIons * pPickUp = m_eqSys->getPickupIons();
          if( pPickUp != NULL )
          {
            if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
            {
              int iRhoPI  = pPickUp->primInterval().begin();
              FORT_OMNIGSSPHERICAL_PI( CHF_FRA(a_W),
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
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }
  }
}


void HelioRealBCProblem::fillGhostCellsWSO(       FArrayBox&      a_W,
                                    const FArrayBox&      a_U,
                                    const int&            a_dir,
                                    const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);  
  int iRhoN  = eqSys->densityIndexPrim(1);  
  int fluids = nFluids();
  int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iHCSb  = (m_eqSys->numTrackingSurfaces()>1 ? iHCS+1 : -1);
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  

  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain;  
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = abs(indW-indD);    
           
        
    if( indW != indD )
    {      
                                   
      if (a_dir == 2)      
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
            
      
    }
       
    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if ( indW != indD )
    {
      
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
            
      
      if (a_dir == 0)
      FORT_WSOGSSPHERICAL(
          CHF_FRA(a_W),            
          CHF_CONST_REAL(m_startBC),                             
          CHF_CONST_REAL(m_startBCPhys),                             
          CHF_CONST_R1D(m_WSOData.time,m_WSOData.nPoints),          
          CHF_CONST_R1D(m_WSOData.Tilt,m_WSOData.nPoints),          
          CHF_CONST_INT(iHCS),          
          CHF_CONST_INT(iHCSb),          
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
      
      if (a_dir == 2)  
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }
  
  }
}

void HelioRealBCProblem::fillGhostCellsCIR(           FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  Box WBox = a_W.box();
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);  
  int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iHCSb  = (m_eqSys->numTrackingSurfaces()>1 ? iHCS+1 : -1);
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  

  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain.domainBox();    
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = abs(indW-indD);    
           
        
    if( indW != indD )
    {      
                                   
      if (a_dir == 2)      
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
            
      
    }
       
    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if ( indW != indD )
    {
      
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
      
      Real time0 = 0.0;
                  
      if (a_dir == 0)
      FORT_CIRGSSPHERICAL(
          CHF_FRA(a_W),            
          CHF_CONST_FRA(m_h5Data1),            
          CHF_CONST_REAL(time0),                                       
          CHF_CONST_REAL(m_dataR0),                                       
          CHF_CONST_INT(iHCS),          
          CHF_CONST_INT(iHCSb),          
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(a_time),
          CHF_BOX(boundaryBox));
      
      if (a_dir == 2)  
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }
  
  }
  
}

void HelioRealBCProblem::fillGhostCellsH5Input(           FArrayBox&      a_W,
                               const FArrayBox&      a_U,
                               const int&            a_dir,
                               const Real&           a_time )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
    

  Box WBox = a_W.box();
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);  
  int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iHCSb  = (m_eqSys->numTrackingSurfaces()>1 ? iHCS+1 : -1);
    
                         // See if this chops off the high side of the input box
  Box tmp  = WBox;
  tmp     &= m_domain;  

  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Cartesian)     
  {
    
  }
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical)
  {           
    int indW,indD,nGS;    
    
    Box WBoxDomain  = WBox;
    WBoxDomain     &= m_domain.domainBox();    
    
    // Outer boundary
    indW = WBox.bigEnd( a_dir );
    indD = WBoxDomain.bigEnd( a_dir );    
    
    nGS  = abs(indW-indD);    
           
        
    if( indW != indD )
    {                                         
      if (a_dir == 2)      
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);                              
    }
       
    // Inner boundary  
    indW = WBox.smallEnd( a_dir );
    indD = WBoxDomain.smallEnd( a_dir );
    
    nGS  = abs(indW-indD);
    
    if ( indW != indD )
    {            
      Box boundaryBox  = WBoxDomain;
      boundaryBox.setRange( a_dir, indD-nGS , nGS );
            
      
      if (a_dir == 0)
      {
        if (m_level > 0)
          MayDay::Abort("HelioRealBCProblem::fillGhostCellsH5Input can not work when m_level > 0");
        
        Real Tref       = eos_AU/(m_lismV*3600.0*24.0);
        Real timeDays  = (a_time-m_startBC)*Tref; // days
        
        bool cirBC      = false;
        bool needRescan = false;
        
        if ((timeDays<0.0) || (timeDays < m_firstDataSetTime))
        {
          cirBC = true; 
          if (procID() == 0)
            pout() << "cir bc" << std::endl;
          
          CH_assert(m_iData1<=0);
          if (m_iData1<0)
          {   
            m_iData1 = 0;            
            readH5DataSet(m_h5fd, m_iData1, m_h5Data1, m_time1);                       
          }
        } else
        if (timeDays>=m_lastDataSetTime)
        {
          cirBC = true;
          if (m_iData2<m_nDataSets-1)
          {            
            m_iData2 = m_nDataSets-1;
            readH5DataSet(m_h5fd, m_iData2, m_h5Data2, m_time2);                       
          }
        } else
        {
          if (m_iData2 < 0) needRescan = true;     
  
          if ((timeDays>m_time2) && (!needRescan))
          {                                  
            CH_assert(m_iData2<m_nDataSets-1);
            Real t;
            // Reading next data set
            FArrayBox fab(m_h5Data2.box(), m_h5Data2.nComp()); 
            readH5DataSet(m_h5fd, m_iData2+1, fab, t);                    
            if (timeDays<=t)
            {
              m_h5Data1.copy(m_h5Data2);
              m_iData1  = m_iData2;
              m_time1   = m_time2;
              m_h5Data2.copy(fab);
              m_iData2++;
              m_time2   = t;
            } else 
            {
              needRescan = true;
            }
          }                      
        }
        
        if (needRescan)
        {
          Real * datasets_time = new Real[m_nDataSets];
          hid_t dataset; herr_t status;  
    
#ifdef H516
  dataset = H5Dopen(m_h5fd, "datasets_time");  
#else
  dataset = H5Dopen2(m_h5fd, "datasets_time",H5P_DEFAULT);
#endif
 
          status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, datasets_time);      
          H5Dclose(dataset);     
          
          CH_assert(timeDays <= datasets_time[m_nDataSets-1]);
          
          if (timeDays < datasets_time[0])
          {
            char buf[40];
            std::string err("Current time is less than time for the first dataset in ");
            err+=m_h5File;
            sprintf(buf,"(%g < %g)",timeDays,datasets_time[0]);
            err+=std::string(buf);            
            MayDay::Error(err.c_str());
          }          
          Real * res = std::lower_bound(datasets_time, datasets_time + m_nDataSets,timeDays);
          int   iRes = res - datasets_time;
          CH_assert((iRes>=0)&&(iRes<m_nDataSets));
          if (iRes == 0)
          {
            m_iData1 = 0;
            m_iData2 = 1;
          } else
          {
            CH_assert((datasets_time[iRes-1]<=timeDays)&&(timeDays<=datasets_time[iRes]));
            m_iData1 = iRes - 1;
            m_iData2 = iRes;
          }
          readH5DataSet(m_h5fd, m_iData1, m_h5Data1, m_time1);                       
          readH5DataSet(m_h5fd, m_iData2, m_h5Data2, m_time2);                               
        }
        
        if (cirBC == true)
        {
          FArrayBox & cirData = (((timeDays<0.0) || (timeDays < m_firstDataSetTime)) ? m_h5Data1 : m_h5Data2);
          //Real        cirTime = (timeDays<0.0 ? m_time1   : m_time2);
		      //cirTime/=Tref;
          Real cirTime = m_startBC;
          if (m_bDataRotatingFrame == false)  cirTime = a_time;
          FORT_CIRGSSPHERICAL(
            CHF_FRA(a_W),            
            CHF_CONST_FRA(cirData),            
            CHF_CONST_REAL(cirTime),                                       
            CHF_CONST_REAL(m_dataR0),                                       
            CHF_CONST_INT(iHCS),          
            CHF_CONST_INT(iHCSb),          
            CHF_CONST_INT(m_level),    
            CHF_CONST_REAL(a_time),
            CHF_BOX(boundaryBox));
        } else
        {        
          CH_assert((m_time1<=timeDays)&&(timeDays<=m_time2));
          
          Real c1 = 1.0 - (timeDays - m_time1)/(m_time2-m_time1);
          Real c2 = 1.0 - (m_time2 - timeDays)/(m_time2-m_time1);
          
          Box b(m_h5Data1.box());
          //b&=boundaryBox; // add limitation theta!!!
          b.setSmall(2,boundaryBox.smallEnd(2));
          b.setBig(  2,boundaryBox.bigEnd(2));
          
          FArrayBox fab(b,m_h5Data1.nComp());
          fab.copy(m_h5Data1);
          fab.mult(c1);
          fab.plus(m_h5Data2,c2);
          
          if ((procID()==0) && (1))
          { 
            pout() << "m_iData1=" << m_iData1 << " c1=" << c1 << ", m_iData2=" << m_iData2 << " c2=" << c2 << std::endl;
          }
                    
          
          if ((procID()==0) && (0))
          {             
            FArrayBox fab_tmp(b,m_h5Data1.nComp());
            fab_tmp.copy(fab);
            
            NodeFArrayBox coords(m_h5Data1.box(),3);          
            m_csh->getNodeCoordsCartesian(coords, m_level);         
            fab_tmp.mult(m_lismN, WRHO);
            fab_tmp.mult(m_lismV/1e+5, WVELX, 3);  
            fab_tmp.mult(m_Bref*m_Bref/1e-12, WPRES);  
            fab_tmp.mult(m_Bref/1e-6 , WBX  , 3);      
          
            static FILE* tfile = OpenTecplotFile("fab.dat","VARIABLES=\"X\" \"Y\" \"Z\" \"rho\" \"vr\" \"vphi\" \"vt\" \"pres\" \"br\" \"bphi\" \"bt\"");    
            WriteFArrayBoxToTecplotFile(tfile, fab_tmp, fab_tmp.box(), Interval(0, fab_tmp.nComp()-1), coords.getFab()); 
            
            static FILE* tfile1 = OpenTecplotFile("m_h5Data1.dat","VARIABLES=\"X\" \"Y\" \"Z\" \"rho\" \"vr\" \"vphi\" \"vt\" \"pres\" \"br\" \"bphi\" \"bt\"");
            if (tfile1 != NULL)
            {
              WriteFArrayBoxToTecplotFile(tfile1, m_h5Data1, m_h5Data1.box(), Interval(0, m_h5Data1.nComp()-1), coords.getFab());
              CloseTecplotFile(tfile1);      
              tfile1 = NULL;
            }
            
            static FILE* tfile2 = OpenTecplotFile("m_h5Data2.dat","VARIABLES=\"X\" \"Y\" \"Z\" \"rho\" \"vr\" \"vphi\" \"vt\" \"pres\" \"br\" \"bphi\" \"bt\"");
            if (tfile2 != NULL)
            {
              WriteFArrayBoxToTecplotFile(tfile2, m_h5Data2, m_h5Data2.box(), Interval(0, m_h5Data2.nComp()-1), coords.getFab());
              CloseTecplotFile(tfile2);      
              tfile2 = NULL;
            }

            
            //CloseTecplotFile(tfile);      
          }
          
          if (m_bDataRotatingFrame == true)                    
            FORT_HDF5INPUTGCROTATING(
              CHF_FRA(a_W),            
              CHF_CONST_FRA(fab),                      
              CHF_CONST_REAL(m_startBC),                                       
              CHF_CONST_REAL(m_dataR0),                                       
              CHF_CONST_INT(iHCS),          
              CHF_CONST_INT(iHCSb),          
              CHF_CONST_INT(m_level),              
              CHF_CONST_REAL(a_time),                                       
              CHF_BOX(boundaryBox));
          else
            FORT_HDF5INPUTGC(
              CHF_FRA(a_W),            
              CHF_CONST_FRA(fab),                                    
              CHF_CONST_REAL(m_dataR0),                                       
              CHF_CONST_INT(iHCS),          
              CHF_CONST_INT(iHCSb),          
              CHF_CONST_INT(m_level),              
              CHF_CONST_REAL(a_time),                                       
              CHF_BOX(boundaryBox));
            
            
        }
          
      }
      
      if (a_dir == 2)  
        HeliosphericProblem::fillGhostCells( a_W, a_U, a_dir, a_time);            
    }
  
  }
  
}



                                                          // Set boundary fluxes
void HelioRealBCProblem::fluxBC(       FArrayBox&      a_F,
                                        FArrayBox&      a_Bn,
                                  const FArrayBox&      a_WMinus,
                                  const FArrayBox&      a_WPlus,
                                  const int&            a_dir,
                                  const Side::LoHiSide& a_side,
                                  const Real&           a_time)
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);
  
  if (m_csh->coordinateSystem() == CoordinateSystemHandler::CS_Spherical) 
  {
    HeliosphericProblem::fluxBC( a_F, a_Bn, a_WMinus, a_WPlus, a_dir, a_side,  a_time);
    return;
  }
  
}


                                                    // Set up initial conditions
void HelioRealBCProblem::initialize( LevelData<FArrayBox>& a_U )
{
  CH_assert(m_isFortranCommonSet == true);
  CH_assert(m_isDefined == true);

  DataIterator dit = a_U.boxLayout().dataIterator();
  
  EqSysMHDMF * eqSys = static_cast<EqSysMHDMF*>(m_eqSys);
  int iRhoN  = eqSys->densityIndexCons(1);
  int iHCS = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  int iHCSb = (m_eqSys->numTrackingSurfaces()>1 ? iHCS+1 : -1);
  int fluids = nFluids();
  
  HeliosphericProblem::initialize( a_U );

                                          // Iterator of all grids in this level
  for (dit.begin(); dit.ok(); ++dit)
  {
                                                     // Storage for current grid
    FArrayBox& U = a_U[dit()];

                                                          // Box of current grid
    Box uBox = U.box();
    uBox &= m_domain.domainBox();
    
    FArrayBox W   (uBox,eqSys->numPrimitives());      
    FArrayBox Ubuf(uBox,eqSys->numStates());      
    
    m_eqSys->stateToPrim(W,U,uBox);
    m_csh->transCartesianVectToCurv(W,uBox,m_level);    
    Ubuf.copy(U,uBox);
    
                                        // Set up initial condition in this grid
    if (m_subproblem == HPBC_GMIRPROBLEM)
    FORT_HELIOGMIRINITSPHERICAL( CHF_CONST_FRA(U), 
                      CHF_CONST_INT(iRhoN ),                     
                      CHF_CONST_INT(fluids),
                      CHF_CONST_INT(m_level),
                      CHF_BOX(uBox) );
                      
    if (m_subproblem == HPBC_SUESSPROBLEM)
    {            
                  
      Real time = 0.0;
      
      FORT_HELIOGSPLASMASPHDEFAULT(
         CHF_FRA(W),     
         CHF_CONST_INT(iHCS),
         CHF_CONST_INT(iHCSb),
         CHF_CONST_INT(m_level),
         CHF_CONST_REAL(time),
         CHF_BOX(uBox));
         
      Real tilt = m_sunTILT;
      
      int kdomain = m_domain.domainBox().size(2);
      FORT_SUESSTHETABC(
            CHF_FRA(W),                    
            CHF_CONST_R1D(m_ThetaDistrib[0].N,kdomain),
            CHF_CONST_R1D(m_ThetaDistrib[0].V,kdomain),
            CHF_CONST_R1D(m_ThetaDistrib[0].T,kdomain),     
            CHF_CONST_REAL(m_ThetaDistrib[0].time),     
            CHF_CONST_R1D(m_ThetaDistrib[1].N,kdomain),
            CHF_CONST_R1D(m_ThetaDistrib[1].V,kdomain),
            CHF_CONST_R1D(m_ThetaDistrib[1].T,kdomain),     
            CHF_CONST_REAL(m_ThetaDistrib[1].time),     
            CHF_CONST_REAL(time),
            CHF_CONST_REAL(time),
            CHF_CONST_REAL(tilt),
            CHF_CONST_INT(iHCS),          
            CHF_CONST_INT(iHCSb),          
            CHF_CONST_INT(m_level),    
            CHF_BOX(uBox));
            
      m_csh->transCurvVectToCartesian(W,uBox,m_level);
      
      m_eqSys->primToState(Ubuf,W,uBox);
      //U.copy(Ubuf,uBox,URHO,uBox,URHO,UNUM);// Copy only plasma data 
      U.copy(Ubuf,uBox);
                          
    }
    if (m_subproblem == HPBC_OMNIPROBLEM)
    {            

      int iHCS   = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
      int iGMIR  = (m_eqSys->numTrackingSurfaces()>1 ? m_eqSys->lvlsStateInterval().begin()+1 : -1);   
      
      Real time = 0.0;
                        
      FORT_OMNIGSSPHERICAL(
          CHF_FRA(W),
          CHF_CONST_REAL(time),
          CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.N,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.V,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.T,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Br,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Bp,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Bt,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.B,m_OMNIData.nPoints),
          CHF_CONST_R1D(m_OMNIData.Tilt,m_OMNIData.nPoints),
          CHF_CONST_INT(iHCS),
          CHF_CONST_INT(iGMIR),    
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(time),
          CHF_BOX(uBox));

      m_csh->transCurvVectToCartesian(W,uBox,m_level);

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          if( pTurbMod->modelID() == TurbulenceModel::TM_BREECH_ET_AL_2008 )
          {
            int iZ2  = pTurbMod->primInterval().begin();
            FORT_OMNIGSSPHERICAL_TM( CHF_FRA(W),
                                     CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
                                     CHF_CONST_R1D(m_OMNIData.N,m_OMNIData.nPoints),
                                     CHF_CONST_INT(iZ2),
                                     CHF_CONST_REAL(time),
                                     CHF_CONST_REAL(time),
                                     CHF_BOX(uBox) );
          }
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          if( pPickUp->modelID() == PickupIons::PI_TWO_EQNS )
          {
            int iRhoPI  = pPickUp->primInterval().begin();
            int sign    =-1;
            int dir     = 0;
            FORT_HELIOGSSPHERICAL_PI( CHF_FRA(W),
                                      CHF_CONST_FRA(U),
                                      CHF_CONST_INT(sign),
                                      CHF_CONST_INT(dir),
                                      CHF_CONST_INT(iRhoPI),
                                      CHF_CONST_INT(m_level),
                                      CHF_CONST_REAL(time),
                                      CHF_BOX(uBox) );
          }
        }
      }

      m_eqSys->primToState(Ubuf,W,uBox);
      U.copy(Ubuf,uBox,URHO,uBox,URHO,UNUM);// Copy only plasma data 

      if( m_eqSys != NULL )
      {
        TurbulenceModel * pTurbMod = m_eqSys->getTurbulenceModel();
        if( pTurbMod != NULL )
        {
          int iRhoZ2 = pTurbMod->consInterval().begin();
          int iNumTM = pTurbMod->numConservative();
          U.copy( Ubuf, uBox, iRhoZ2, uBox, iRhoZ2, iNumTM );
        }

        PickupIons * pPickUp = m_eqSys->getPickupIons();
        if( pPickUp != NULL )
        {
          int iRhoPI = pPickUp->consInterval().begin();
          int iNumPI = pPickUp->numConservative();
          U.copy( Ubuf, uBox, iRhoPI, uBox, iRhoPI, iNumPI );
        }
      }
    }
    
    if (m_subproblem == HPBC_HDF5INPUT)
    {              
      m_iData1 = 0;
      readH5DataSet(m_h5fd, m_iData1, m_h5Data1, m_time1);
      
      m_iData2 = 1;
      readH5DataSet(m_h5fd, m_iData2, m_h5Data2, m_time2);                    
    }
        
    if ((m_subproblem == HPBC_CIRPROBLEM) || (m_subproblem == HPBC_HDF5INPUT))
    {
        Real time = 0.0;
        
        FORT_CIRGSSPHERICAL(
          CHF_FRA(W),            
          CHF_CONST_FRA(m_h5Data1),            
          CHF_CONST_REAL(time),
          CHF_CONST_REAL(m_dataR0),                                       
          CHF_CONST_INT(iHCS),          
          CHF_CONST_INT(iHCSb),          
          CHF_CONST_INT(m_level),    
          CHF_CONST_REAL(time),
          CHF_BOX(uBox));
          
        m_csh->transCurvVectToCartesian(W,uBox,m_level);      
        m_eqSys->primToState(Ubuf,W,uBox);
        
        // Plasma data
        U.copy(Ubuf,uBox,URHO,uBox,URHO,UNUM);
        // Levelset data
        if (m_eqSys->numTrackingSurfaces()>0)
        U.copy(Ubuf,uBox,iHCS,uBox,iHCS,m_eqSys->numTrackingSurfaces());                    
    }    
                                 
  }
}
                                              // Problem specific postprocessing
void HelioRealBCProblem::postprocessing(       FArrayBox & a_U,
                                          const FArrayBox & a_W,
                                          const Real      & a_dt,
                                          const Real      & a_time,                                          
                                          const Box       & a_box       )
{  

  
}

//                              Creates tagged cells for dynamic mesh refinement
void HelioRealBCProblem::tagCells( const FArrayBox&  a_U,
                              const Box&        a_box,
                              IntVectSet& a_tags )
{  
  //if (m_subproblem == HPBC_CIRPROBLEM)
  //  tagCellsCIR(a_U,a_box,a_tags);
  
  tagCellsCME(a_U,a_box,a_tags);
    
}

void HelioRealBCProblem::tagCellsCME(const FArrayBox&  a_U,
                    const Box&        a_box,
                          IntVectSet& a_tags)
{
  BoxIterator bit( a_box ); RealVect cv_coord;
  
  Real cme_tr, rho, ux, uy, uz;
  
  int iCME = (m_eqSys->numTrackingSurfaces()>0 ? m_eqSys->lvlsStateInterval().begin() : -1);
  if (iCME < 0) return;

  for( bit.begin(); bit.ok(); ++bit )
  {
    const IntVect& iv = bit();
    
    cme_tr = a_U.get(iv,iCME);
        
    m_csh->getCellCenter(cv_coord,iv,m_level);    
    if (cv_coord[0]>0.18)
    {         
      if (cme_tr > 0) a_tags |= iv;            
      else
      {
        rho    = a_U.get(iv,URHO);
        ux     = 1e-5*m_lismV*a_U.get(iv,UMOMX)/rho;
        uy     = 1e-5*m_lismV*a_U.get(iv,UMOMY)/rho;
        uz     = 1e-5*m_lismV*a_U.get(iv,UMOMZ)/rho;
        
        if ((ux*ux + uy*uy + uz*uz)>500*500) a_tags |= iv;
      }
    }
  }    
  
}

// Tags CIR region 
void HelioRealBCProblem::tagCellsCIR(const FArrayBox&  a_U,
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
    
    Real theta1,theta2;
    // SC22 min
    theta1 = (90.0-(45.0 + 5.0))*d_PI_180;
    theta2 = (90.0+(45.0 + 5.0))*d_PI_180;
    
    // SC23 min
    theta1 = (90.0-(57.0))*d_PI_180;
    theta2 = (90.0+(57.0))*d_PI_180;
    
    Real th_axis = 17.0;
    Real thetaReg2Min = (th_axis)*d_PI_180;
    Real thetaReg2Max = (180.0-th_axis)*d_PI_180;
        
        
    if (m_level == 0)
    {            
      Real phi    = 9.0;
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;
              
      if ((reg == 3) && (cv_coord[0]>20.0) &&
          (cv_coord[2]>theta1) && (cv_coord[2]<theta2) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      if ((reg == 3) && (cv_coord[0]>65.0) &&
          (cv_coord[2]>theta1-5.0*d_PI_180) && (cv_coord[2]<theta2+5.0*d_PI_180) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      if ((reg == 2) && 
          (cv_coord[2]>thetaReg2Min) && (cv_coord[2]<thetaReg2Max) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
              
    }          
    
    
    if (m_level == 1)
    {        
    
      Real phi    = 50.0;      
      Real phi1   = (90.0  - phi)*d_PI_180;
      Real phi2   = (270.0 + phi)*d_PI_180;

                          
      if ((reg == 3) && (cv_coord[0]>65.0) &&
          (cv_coord[2]>theta1-8.0*d_PI_180) && (cv_coord[2]<theta2+8.0*d_PI_180) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
      
      if ((reg == 2) && 
          (cv_coord[2]>thetaReg2Min) && (cv_coord[2]<thetaReg2Max) && 
          ((cv_coord[1]<phi1) || (cv_coord[1]>phi2)) ) 
      {
        a_tags |= iv;
        continue;
      }
              
    }          

  }
}



//                             Check geometrical limitations for grid adaptation
void HelioRealBCProblem::lockedCellsRegrid( BaseFab<int> & a_flag,
                                       const FArrayBox&  a_U,
                                       const Box&     a_box)
{
  a_flag.setVal(0);
        
}                            

//                            Return boundary condition flags for all boundaries
void HelioRealBCProblem::getBCFlags( eBoundaryConditions leftBC,
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
PhysProblem::eBoundaryConditions HelioRealBCProblem::getBCFlags( int a_dir, Side::LoHiSide a_sd )
{
  if( a_dir == 0 )
  {
    return ( a_sd == Side::Lo ) ? BC_Continuous : BC_Fixed;
  } else {
    return BC_Continuous;
  }
}

void HelioRealBCProblem::readOdstrcilData()
{
	int j,k,var;
	FILE* data_file=fopen("outbnd.dat","r");
  
  // Nlat - k-index (shirota)
  // Nlon - j-index (dolgota)
  int Nlat, Nlon; float fbuf;  
	fscanf(data_file,"%d %d",&Nlat,&Nlon);
  
  CH_assert(m_domain.size(1)==Nlon);
  CH_assert(m_domain.size(2)==Nlat);
  
  m_OdstrcilData.define(Box(
    IntVect(D_DECL(0,0,0)),
    IntVect(D_DECL(0,Nlon-1,Nlat-1))),8);
  
  for (k=0;k<Nlat;k++) 
    fscanf(data_file,"%f",&fbuf);
  for (j=0;j<Nlon;j++) 
    fscanf(data_file,"%f",&fbuf);
    
  
  for (var=0; var<8; var++)
  for (j=0; j<Nlon; j++)
  for (k=0; k<Nlat; k++)
  {
    fscanf(data_file,"%f",&fbuf);
    //m_OdstrcilData.set(IntVect(D_DECL(0,(Nlon/2+j)%Nlon,k)),var,(Real)fbuf);            
    m_OdstrcilData.set(IntVect(D_DECL(0,j,k)),var,(Real)fbuf);            
  }
  fclose(data_file);	  
        
  FORT_ODSTRCILPREPARE(
     CHF_FRA(m_OdstrcilData),
     CHF_BOX(m_OdstrcilData.box()));
        
}

void HelioRealBCProblem::readGMIRData()
{	
	FILE* data_file=fopen("GMIRdata.txt","r");
  
  if (data_file==NULL) MayDay::Abort("GMIRData.txt not found");
  
  const int NBuf = 200;  
  char buffer[NBuf];
  
  int i;
  
  double time, N,V,T,B, time0 = 0.0;
  
  fgets(buffer, NBuf-1, data_file);
  
  fscanf(data_file,"%d",&m_GMIRData.nPoints);
  
  m_GMIRData.time = new Real[m_GMIRData.nPoints];
  m_GMIRData.N = new Real[m_GMIRData.nPoints];
  m_GMIRData.V = new Real[m_GMIRData.nPoints];
  m_GMIRData.B = new Real[m_GMIRData.nPoints];
  m_GMIRData.T = new Real[m_GMIRData.nPoints];
  
  for (i=0;i<m_GMIRData.nPoints;i++)
  {
    fscanf(data_file,"%lf %lf %lf %lf %lf",&time, &B, &N, &T, &V);
    
    if (i==0) time0 = time;
    
    m_GMIRData.time[i] = time;
    m_GMIRData.B[i]    = B;
    m_GMIRData.N[i]    = N;
    m_GMIRData.T[i]    = T;
    m_GMIRData.V[i]    = V;
    
    m_GMIRData.B[i]*=10.0;  // transfer to mG
    m_GMIRData.B[i]*=1e-6;  // transfer to G
    m_GMIRData.V[i]*=1.1;  
    m_GMIRData.time[i]-=time0;        
  }
  
  fclose(data_file);  
        
}

void HelioRealBCProblem::readThetaDistrib()
{
  FILE* Nf=fopen("N10smooth.dat","r");
  if (Nf==NULL) MayDay::Abort("N10smooth.dat not found");  
  FILE* Tf=fopen("Tarray10.dat","r");
  if (Tf==NULL) MayDay::Abort("Tarray10.dat not found");
  FILE* Vf=fopen("vel10smooth.dat","r");
  if (Vf==NULL) MayDay::Abort("vel10smooth.dat not found");
        
  m_nThetaDistrib = 240;
  const int kread = 181;
  Real RealBuf[kread];
    
  
  int dir = 2;
  
  IntVect iv_off(IntVect::Zero);
  iv_off[dir]=1;
    
  Box cBox(m_domain.domainBox().smallEnd()*iv_off, 
           m_domain.domainBox().bigEnd()  *iv_off);    
  
  FArrayBox thetac(cBox,1);
        
  m_csh->getCellCenters(thetac,thetac.box(),dir,m_level);
            
  int i,k;
  int kdomain = m_domain.domainBox().size(2);
  
  m_startBCPhys = 1991+15.0/365.0;
  
  Real dt = (1.0/12.0); // years
  Real time = 0.0;
  
  m_ThetaDistrib = new ThetaDistrib[m_nThetaDistrib];
  for (i=0;i<m_nThetaDistrib;i++)
  {
    m_ThetaDistrib[i].time = time;
    m_ThetaDistrib[i].N    = new Real[kdomain];
    m_ThetaDistrib[i].V    = new Real[kdomain];
    m_ThetaDistrib[i].T    = new Real[kdomain];
    for (k=0;k<kread;k++)
    {
      fscanf(Nf,"%lf",&RealBuf[k]);            
    }
    FORT_THETAINTERPOLATEDATA(                    
          CHF_R1D(m_ThetaDistrib[i].N,kdomain),
          CHF_CONST_R1D(RealBuf,kread),
          CHF_CONST_FRA1(thetac,0));                  
    for (k=0;k<kread;k++)
    {
      fscanf(Tf,"%lf",&RealBuf[k]);            
    }
    FORT_THETAINTERPOLATEDATA(                    
          CHF_R1D(m_ThetaDistrib[i].T,kdomain),
          CHF_CONST_R1D(RealBuf,kread),
          CHF_CONST_FRA1(thetac,0));                  
    for (k=0;k<kread;k++)
    {
      fscanf(Vf,"%lf",&RealBuf[k]);            
    } 
    FORT_THETAINTERPOLATEDATA(                    
          CHF_R1D(m_ThetaDistrib[i].V,kdomain),
          CHF_CONST_R1D(RealBuf,kread),
          CHF_CONST_FRA1(thetac,0));                     
    time+=dt;  
  }
 
  
  fclose(Nf);  
  fclose(Tf);  
  fclose(Vf);    

  readWSOData(); 
}

void HelioRealBCProblem::readOMNIData()
{	
	FILE* data_file=fopen(m_OMNIDataFile.c_str(),"r");  
  
  if (data_file==NULL) 
  {
    std::string err(m_OMNIDataFile+" not found");      
    MayDay::Abort(err.c_str());
  }
      
  int i;
  
  double time, N,V,T,Br,Bp,Bt,B;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double year, day, hour;

//  const int iOldOMNIFile = 0;
    
  //fscanf(data_file,"%d",&m_OMNIData.nPoints);
  m_OMNIData.nPoints = 0;
  
  while (!feof(data_file))
  {
    if( m_iOldOMNIFile ==  1 )
    {
      fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour, &tmp1, &tmp2, 
              &Br, &Bp, &Bt, &B, &V, &tmp4, &tmp5, &N, &T );
    } else
    {
      fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour,
              &Br, &Bp, &Bt, &B, &V, &N, &T, &tmp1 );
    }
    if ((T<9999998.0) && (N<999.0) && (V<9998.0)) m_OMNIData.nPoints++;
  }    
  
  fseek(data_file,0,SEEK_SET);
  
  m_OMNIData.time = new Real[m_OMNIData.nPoints];
  m_OMNIData.N    = new Real[m_OMNIData.nPoints];
  m_OMNIData.V    = new Real[m_OMNIData.nPoints];
  m_OMNIData.Br   = new Real[m_OMNIData.nPoints];
  m_OMNIData.Bp   = new Real[m_OMNIData.nPoints];
  m_OMNIData.Bt   = new Real[m_OMNIData.nPoints];
  m_OMNIData.B    = new Real[m_OMNIData.nPoints];
  m_OMNIData.Tilt = new Real[m_OMNIData.nPoints];
  m_OMNIData.T    = new Real[m_OMNIData.nPoints];
  
  i = 0;
  while (!feof(data_file))  
  {
    if( m_iOldOMNIFile ==  1 )
    {
      fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour, &tmp1, &tmp2, 
              &Br, &Bp, &Bt, &B, &V, &tmp4, &tmp5, &N, &T );
    } else
    {
      fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour,
              &Br, &Bp, &Bt, &B, &V, &N, &T, &tmp1 );
    }
      
    if  (  !((T<9999998.0) && (N<999.0) && (V<9998.0)) ) continue;

    time = year + (day-1)/365 + hour/(365*24.0); 
    
    if (i>0) 
      CH_assert(time>m_OMNIData.time[i-1]);
            
    m_OMNIData.time[i] = time;    
    m_OMNIData.N[i]    = N;
    m_OMNIData.T[i]    = T;
    m_OMNIData.V[i]    = V;
    m_OMNIData.Br[i]   = Br;
    m_OMNIData.Bp[i]   = Bp;
    m_OMNIData.Bt[i]   =-Bt;
    m_OMNIData.B[i]    = B;
    
    m_OMNIData.Br[i]*=10.0;  // transfer to mG
//    m_OMNIData.Br[i]*=1e-6;  // transfer to G
    m_OMNIData.Bp[i]*=10.0;  
//    m_OMNIData.Bp[i]*=1e-6;  
    m_OMNIData.Bt[i]*=10.0;  
//    m_OMNIData.Bt[i]*=1e-6;          

    m_OMNIData.B[i]*=10.0;  // transfer to mG
    
    i++;
  }
  
  fclose(data_file);  
  
  int nHCSpoints;
  FILE* HCS_file=fopen("WSOdata.txt","r");
  if (HCS_file==NULL) MayDay::Abort("WSOdata.txt not found");
  fscanf(HCS_file,"%d",&nHCSpoints);
  
  Real * WSOTilt = new Real[nHCSpoints];
  Real * WSOTime = new Real[nHCSpoints];
  Real R_av,R_n,R_s,L_av,L_n,L_s;

  for (i=0;i<nHCSpoints;i++)
  {
    fscanf(HCS_file,"%lf %lf %lf %lf %lf %lf %lf", &time, &R_av, &R_n, &R_s, &L_av, &L_n, &L_s);
    WSOTime[i] = time;    
    WSOTilt[i] = L_av;    
  }        
  fclose(HCS_file);
  
  FORT_OMNIINTERPOLATETILT(                    
      CHF_R1D(m_OMNIData.Tilt,m_OMNIData.nPoints),
      CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
      CHF_CONST_R1D(WSOTilt,nHCSpoints),
      CHF_CONST_R1D(WSOTime,nHCSpoints));
          
  delete[] WSOTilt;
  delete[] WSOTime;    
  
  for (i=1;i<m_OMNIData.nPoints;i++) m_OMNIData.time[i] -= m_OMNIData.time[0];  
  
  m_startBCPhys = m_OMNIData.time[0];
  
  m_OMNIData.time[0] = 0.0;
  
  if (procID() == 0) remove("fort.66");    
          
}

void HelioRealBCProblem::readOMNIDataTurb()
{	
	FILE* data_file=fopen(m_OMNIDataFile.c_str(),"r");

  if( data_file == NULL )
  {
    std::string err(m_OMNIDataFile+" not found");
    MayDay::Abort(err.c_str());
  }

  int i;

  double time, N, V, T, Br, Bp, Bt, B, Z2;
  double tmp1, tmp2, tmp3, tmp4, tmp5;
  double year, day, hour;

  m_OMNIData.nPoints = 0;

  while( !feof( data_file ) )
  {
    fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour,
            &Br, &Bp, &Bt, &B, &V, &N, &T, &Z2 );
    if( (T < 9999998.0) && (N < 999.0) && (V < 9998.0) ) m_OMNIData.nPoints++;
  }

  fseek(data_file,0,SEEK_SET);

  m_OMNIData.time = new Real[m_OMNIData.nPoints];
  m_OMNIData.N    = new Real[m_OMNIData.nPoints];
  m_OMNIData.V    = new Real[m_OMNIData.nPoints];
  m_OMNIData.Br   = new Real[m_OMNIData.nPoints];
  m_OMNIData.Bp   = new Real[m_OMNIData.nPoints];
  m_OMNIData.Bt   = new Real[m_OMNIData.nPoints];
  m_OMNIData.B    = new Real[m_OMNIData.nPoints];
  m_OMNIData.Tilt = new Real[m_OMNIData.nPoints];
  m_OMNIData.T    = new Real[m_OMNIData.nPoints];
  m_OMNIData.Z2   = new Real[m_OMNIData.nPoints];

  i      = 0;
  while (!feof(data_file))  
  {
    fscanf( data_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &year, &day, &hour,
            &Br, &Bp, &Bt, &B, &V, &N, &T, &Z2 );

    if( (T >= 9999998.0) || (N >= 998.0) || (V >= 9998.0) ) continue;

    int iYear  = (int)year;
    if( iYear%4 == 0 )
    {
      time = year + (day - 1.0 + hour/24.0)/366.0;
    } else {
      time = year + (day - 1.0 + hour/24.0)/365.0;
    }

    if( i > 0 )
      CH_assert(time>m_OMNIData.time[i-1]);
//                                                                transfer to mG
    Br  *=10.0;
    Bp  *=10.0;
    Bt  *=10.0;
    B   *=10.0;

    m_OMNIData.time[i] = time;
    m_OMNIData.N[i]    = N;
    m_OMNIData.T[i]    = T;
    m_OMNIData.V[i]    = V;
    m_OMNIData.Br[i]   = Br;
    m_OMNIData.Bp[i]   = Bp;
    m_OMNIData.Bt[i]   =-Bt;
    m_OMNIData.B[i]    = B;
    m_OMNIData.Z2[i]   = Z2;

    i++;
  }

  fclose(data_file);

  //                                                                 read WSO data
  int nHCSpoints;
  FILE* HCS_file=fopen("WSOdata.txt","r");
  if( HCS_file == NULL ) MayDay::Abort("WSOdata.txt not found");

  fscanf(HCS_file,"%d",&nHCSpoints);

  Real * WSOTilt = new Real[nHCSpoints];
  Real * WSOTime = new Real[nHCSpoints];

  Real R_av, R_n, R_s, L_av, L_n, L_s;

  for( i=0; i<nHCSpoints; i++ )
  {
    fscanf(HCS_file,"%lf %lf %lf %lf %lf %lf %lf", &time, &R_av, &R_n, &R_s, &L_av, &L_n, &L_s);
    WSOTime[i] = time;
    WSOTilt[i] = L_av;
  }

  fclose(HCS_file);

  FORT_OMNIINTERPOLATETILT(
      CHF_R1D(m_OMNIData.Tilt,m_OMNIData.nPoints),
      CHF_CONST_R1D(m_OMNIData.time,m_OMNIData.nPoints),
      CHF_CONST_R1D(WSOTilt,nHCSpoints),
      CHF_CONST_R1D(WSOTime,nHCSpoints));

  delete[] WSOTilt;
  delete[] WSOTime;

  for (i=1;i<m_OMNIData.nPoints;i++) m_OMNIData.time[i] -= m_OMNIData.time[0];

  m_startBCPhys = m_OMNIData.time[0];

  m_OMNIData.time[0] = 0.0;
}

void HelioRealBCProblem::readWSOData()
{	      
              
  FILE* HCS_file=fopen("WSOdata.txt","r");
  if (HCS_file==NULL) MayDay::Abort("WSOdata.txt not found");
  fscanf(HCS_file,"%d",&m_WSOData.nPoints);
  
  m_WSOData.time = new Real[m_WSOData.nPoints];
  m_WSOData.Tilt = new Real[m_WSOData.nPoints];
  
  int i;
  Real time, R_av,R_n,R_s,L_av,L_n,L_s;

  for (i=0;i<m_WSOData.nPoints;i++)
  {
    fscanf(HCS_file,"%lf %lf %lf %lf %lf %lf %lf", &time, &R_av, &R_n, &R_s, &L_av, &L_n, &L_s);
    m_WSOData.time[i] = time;    
    m_WSOData.Tilt[i] = R_av;    
  }        
  fclose(HCS_file);    
              
  for (i=0;i<m_WSOData.nPoints;i++) m_WSOData.time[i] -= m_startBCPhys;      
  
  if (procID() == 0) remove("fort.66");
    
}


void HelioRealBCProblem::readCIRData()
{	      
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();  
  
  if (CoordinateSystem != CoordinateSystemHandler::CS_Spherical) return;
  if (m_level>0) return;

#ifdef CH_USE_HDF5
  hid_t attr,h5f,dataset; herr_t status;

  h5f = H5Fopen(m_h5File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
  
  m_nDataSets = 1;
  
#ifdef H516
  attr = H5Aopen_name(h5f,"r0");
#else
  attr = H5Aopen_by_name(h5f, ".", "r0", H5P_DEFAULT,H5P_DEFAULT);
#endif

  H5Aread(attr, H5T_NATIVE_DOUBLE, &m_dataR0);
  H5Aclose(attr);
  
  int domain[2];  

#ifdef H516  
  attr = H5Aopen_name(h5f, "domain"); 
#else
  attr = H5Aopen_by_name(h5f, ".", "domain", H5P_DEFAULT,H5P_DEFAULT);
#endif

  H5Aread(attr, H5T_NATIVE_INT, domain);   
  H5Aclose(attr);      
  
  const Box & dBox = m_domain.domainBox();   
  if ((domain[0]!=dBox.size(1)) || (domain[1]!=dBox.size(2)) )
    MayDay::Error("HelioRealBCProblem::readCIRData - problem domain does not match data in h5 file");       
                                                   
  int nComp;  

#ifdef H516
  attr = H5Aopen_name(h5f, "num_components"); 
#else
  attr = H5Aopen_by_name(h5f, ".", "num_components", H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_INT, &nComp);   
  H5Aclose(attr);
  if (nComp<UNUM+m_eqSys->numTrackingSurfaces())
    MayDay::Error("HelioRealBCProblem::readCIRData - not enough components in h5 file");       
      
        
  Box b( IntVect(D_DECL(-1, dBox.smallEnd()[1], dBox.smallEnd()[2])), 
         IntVect(D_DECL(-1, dBox.bigEnd()[1],   dBox.bigEnd()[2])));   
  m_h5Data1.define(b,nComp);
      
  m_iData1 = 0;
  readH5DataSet(h5f, m_iData1, m_h5Data1, m_time1);
  
  H5Fclose(h5f);       
                                   
#endif
  
  
  if ((procID()==0) && (1))
  {    
    NodeFArrayBox coords(m_h5Data1.box(),3);          
    m_csh->getNodeCoordsCartesian(coords, m_level);   
    
    FArrayBox tmpFAB(m_h5Data1.box(),m_h5Data1.nComp());
    tmpFAB.copy(m_h5Data1);

    tmpFAB.mult(m_lismN, WRHO);
    tmpFAB.mult(m_lismV, WVELX, 3);  
    tmpFAB.mult(m_Bref*m_Bref, WPRES);  
    tmpFAB.mult(m_Bref , WBX  , 3);      
      

    FILE* tfile = OpenTecplotFile("readCIR.dat","VARIABLES=\"X\" \"Y\" \"Z\" \"rho\" \"vr\" \"vphi\" \"vt\" \"pres\" \"br\" \"bphi\" \"bt\"");    
    WriteFArrayBoxToTecplotFile(tfile, tmpFAB, tmpFAB.box(), Interval(0, tmpFAB.nComp()-1), coords.getFab()); 
    CloseTecplotFile(tfile);      
  }
    
}

void HelioRealBCProblem::prepareHDF5InputProblem()
{
  CoordinateSystemHandler::eCoordinateSystem CoordinateSystem = m_csh->coordinateSystem();  
  
  if (CoordinateSystem != CoordinateSystemHandler::CS_Spherical) return;
  if (m_level>0) return;

#ifdef CH_USE_HDF5
  hid_t attr,dataset; herr_t status;

  m_h5fd = H5Fopen(m_h5File.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);  
  
#ifdef H516
  attr = H5Aopen_name(m_h5fd,"r0");
#else
  attr = H5Aopen_by_name(m_h5fd, ".", "r0", H5P_DEFAULT,H5P_DEFAULT);
#endif

  H5Aread(attr, H5T_NATIVE_DOUBLE, &m_dataR0);
  H5Aclose(attr);
  
#ifdef H516
  attr = H5Aopen_name(m_h5fd,"num_datasets");
#else
  attr = H5Aopen_by_name(m_h5fd, ".", "num_datasets", H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_INT, &m_nDataSets);
  H5Aclose(attr);    
  
  int domain[2];    

#ifdef H516
  attr = H5Aopen_name(m_h5fd, "domain"); 
#else
  attr = H5Aopen_by_name(m_h5fd, ".", "domain", H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_INT, domain);   
  H5Aclose(attr);      
  
  const Box & dBox = m_domain.domainBox();   
  if ((domain[0]!=dBox.size(1)) || (domain[1]!=dBox.size(2)) )
    MayDay::Error("HelioRealBCProblem::readCIRData - problem domain does not match data in h5 file");       
                                                   
  int nComp;  

#ifdef H516
  attr = H5Aopen_name(m_h5fd, "num_components"); 
#else
  attr = H5Aopen_by_name(m_h5fd, ".", "num_components", H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_INT, &nComp);   
  H5Aclose(attr);
  if (nComp<UNUM+m_eqSys->numTrackingSurfaces())
    MayDay::Error("HelioRealBCProblem::readCIRData - not enough components in h5 file");       
      
        
  Box b( IntVect(D_DECL(-1, dBox.smallEnd()[1], dBox.smallEnd()[2])), 
         IntVect(D_DECL(-1, dBox.bigEnd()[1],   dBox.bigEnd()[2])));   
  m_h5Data1.define(b,nComp);
  m_h5Data2.define(b,nComp);
    
      
  m_iData1 = -1;
  m_iData2 = -1;    
  
  char dataspace_name[20];
    
  sprintf(dataspace_name,"data%d",0);

#ifdef H516
  dataset = H5Dopen(m_h5fd, dataspace_name);     
  attr = H5Aopen_name(dataset,"time");
#else
  dataset = H5Dopen2(m_h5fd, dataspace_name,H5P_DEFAULT);
  attr = H5Aopen_by_name(dataset,".", "time",H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_DOUBLE, &m_firstDataSetTime);
  H5Aclose(attr);  
  H5Dclose(dataset);     
  
  sprintf(dataspace_name,"data%d",m_nDataSets-1);

#ifdef H516
  dataset = H5Dopen(m_h5fd, dataspace_name);     
  attr = H5Aopen_name(dataset,"time");
#else
  dataset = H5Dopen2(m_h5fd, dataspace_name,H5P_DEFAULT);
  attr = H5Aopen_by_name(dataset,".", "time",H5P_DEFAULT,H5P_DEFAULT);
#endif 

  H5Aread(attr, H5T_NATIVE_DOUBLE, &m_lastDataSetTime);
  H5Aclose(attr);  
  H5Dclose(dataset);     
                                       
#endif
  if (m_nDataSets<2)
  {
    MayDay::Error("Number of data sets for the HDF5 Input subproblem should be greater then 2");
  }
}

void HelioRealBCProblem::readH5DataSet(hid_t a_h5f, int a_iDataSet, FArrayBox & a_data, Real & a_time)
{
#ifdef CH_USE_HDF5
  hid_t attr,dataset; herr_t status;  

  char dataspace_name[20];
  sprintf(dataspace_name,"data%d",a_iDataSet);
    
#ifdef H516
  dataset = H5Dopen(a_h5f, dataspace_name);   
#else
  dataset = H5Dopen2(a_h5f, dataspace_name,H5P_DEFAULT);
#endif 

  status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, a_data.dataPtr());      

#ifdef H516
  attr = H5Aopen_name(dataset,"time");
#else
  attr = H5Aopen_by_name(dataset,".","time",H5P_DEFAULT,H5P_DEFAULT);
#endif

  H5Aread(attr, H5T_NATIVE_DOUBLE, &a_time);
  H5Aclose(attr);
  
  H5Dclose(dataset);     
        
#endif

  if ((procID() == 0) && (0))
  {
    NodeFArrayBox coords(a_data.box(),3);          
    m_csh->getNodeCoordsCartesian(coords, m_level);         
    NodeFArrayBox coords_cv(a_data.box(),3);          
    m_csh->getNodeCoords(coords_cv, m_level);         
        
    FILE* tfile = OpenTecplotFile(strcat(dataspace_name,"input.dat"),a_data.nComp()+SpaceDim);    
    WriteFArrayBoxToTecplotFile(tfile, a_data, a_data.box(), Interval(0, a_data.nComp()-1), coords.getFab(), coords_cv.getFab());
    CloseTecplotFile(tfile);      
  }


  a_data.mult(1.0/m_lismN, WRHO);
  a_data.mult(1.0/m_lismV, WVELX, 3);  
  a_data.mult(1.0/(m_Bref*m_Bref), WPRES);  
  a_data.mult(1.0/m_Bref , WBX  , 3);      
 
}

// For a tracking surface 'a_s' returns first component of the velocity field for the level set method
int HelioRealBCProblem::lsIndexField(int a_s)
{
  return WVELX;
}

time_t getPosixTime(Real a_time)
{
  int year     = (int)(floor(a_time));
  int daysYear = (year % 4 == 0 ? 366 : 365);
    
  tm time_tm;
  bzero(&time_tm,sizeof(tm));
  
  time_tm.tm_year = year - 1900;
  time_tm.tm_sec  = (int)((a_time - year)*daysYear*24.0*3600.0);
  time_tm.tm_mday = 1;
  
  time_t posix_time = mktime(&time_tm);
  
  return posix_time;  
}

Real HelioRealBCProblem::getPhysTime(Real a_time)
{      
  Real   timeRef      = eos_AU/m_lismV; // sec    
  time_t startBC      = getPosixTime(m_startBCPhys);
  Real   timeSec      = (a_time-m_startBC)*timeRef; // time in sec
  time_t timeSec_t    = (time_t)(floor(timeSec));
  
  time_t curTime      = startBC + timeSec_t;
  
  tm curTime_tm;
  //gmtime_r(&curTime,&curTime_tm);
  localtime_r(&curTime,&curTime_tm);
  
  
  int daysYear = (curTime_tm.tm_year % 4 == 0 ? 366 : 365);
  
  Real PhysTime = curTime_tm.tm_year + 1900.0 + 
                  (Real)(curTime_tm.tm_yday)/(Real)(daysYear) +
                  (Real)(curTime_tm.tm_hour)/(Real)(daysYear*24.0) +
                  (Real)(curTime_tm.tm_min )/(Real)(daysYear*24.0*60.0) +
                  (Real)(curTime_tm.tm_sec+timeSec-timeSec_t)/(Real)(daysYear*24.0*3600.0);
                        
  return PhysTime;                 
    
}

Real HelioRealBCProblem::probeFilesDimensiolessDt(Real a_dt)
{  
  // For heliospheric problems timeDt is given in hours, convert it to dimensionless units
  
  Real timeRef    = eos_AU/m_lismV; // sec    
  Real dimlessDt  = a_dt*3600.0/timeRef;
  return dimlessDt;
  
/*
  Real dt = m_probeDt;   
  
  dt /= 24.0*365.0; // dt in years
  Real time = getPhysTime(a_time);
  if ((time - m_last1dprobe) > dt)
  {
    m_last1dprobe = time;
    return true;
  } else return false;
  */
}

void HelioRealBCProblem::probeFilesHeader(int & a_numVars, std::string& a_header)
{  
  int i;
  
  a_numVars = m_eqSys->numPrimitives() + 
              numPlotVars() +
              getSourceCalculator()->numPlotVars();
  
  a_header = "VARIABLES=\"time\" \"Distance\" \"lon\" \"lat\"";  
              
  Vector<std::string> varNames;
  varNames.reserve(a_numVars);        
                        
  varNames = m_eqSys->primitiveNames();
  for (i = 0; i < (int)varNames.size(); ++i)
  {
    a_header += " \"";
    a_header += varNames[i];
    a_header += "\"";
  }
  varNames = plotNames();
  for (i = 0; i < (int)varNames.size(); ++i)
  {
    a_header += " \"";
    a_header += varNames[i];
    a_header += "\"";
  }
  varNames = getSourceCalculator()->plotNames();
  for (i = 0; i < (int)varNames.size(); ++i)
  {
    a_header += " \"";
    a_header += varNames[i];
    a_header += "\"";
  }  
  a_header += "\nZONE\nDT=(";
  for (i = 0; i < 4 + a_numVars; ++i) a_header += "DOUBLE ";
  a_header += ")";
}

/// Returns curvilinear coordinates of probes and corresponding file names  
/*void HelioRealBCProblem::probe1DPoints(Vector<RealVect> & a_points,                                                                              
                                       Real a_time)
{
  // Do not forget to modify needToWrite1DProbe
  
  std::map<Real,RealVect>::const_iterator vpos_iter;
  Real time = getPhysTime(a_time);
  RealVect pAU;    
  
  if ((m_subproblem == HPBC_OMNIPROBLEM) && (0))
  {        
    pAU[1] = 0.01;
    pAU[2] = d_PI_2;    
    
    pAU[0] = 15.0;
    a_points.push_back(pAU);          
    
    pAU[0] = 40.0;
    a_points.push_back(pAU);        
    
    pAU[0] = 60.0;
    a_points.push_back(pAU);        
  }
  
  if ((m_subproblem == HPBC_CIRPROBLEM) && (0))
  {  
    
    
    //TSb3.dat
    pAU[0] = 104.0;
    pAU[1] = 0.01;
    pAU[2] = 1.12;             
    a_points.push_back(pAU);
    
    
  }
  
  if (m_Vprobe == true)
  {
    vpos_iter = m_V1pos.upper_bound(time);
    pAU[0] = vpos_iter->second[0];
    pAU[1] = vpos_iter->second[1];
    pAU[2] = vpos_iter->second[2];         
    a_points.push_back(pAU);    
    
    vpos_iter = m_V2pos.upper_bound(time);
    pAU[0] = vpos_iter->second[0];
    pAU[1] = vpos_iter->second[1];
    pAU[2] = vpos_iter->second[2];         
    a_points.push_back(pAU);    
  }
  
  if (m_Uprobe == true)
  {
    vpos_iter = m_Ulyssespos.upper_bound(time);
    pAU[0] = vpos_iter->second[0];
    pAU[1] = vpos_iter->second[1];
    pAU[2] = vpos_iter->second[2];         
    a_points.push_back(pAU);    
  }
  
  if (m_Earthprobe == true)
  {
    vpos_iter = m_Earthpos.upper_bound(time);
    pAU[0] = vpos_iter->second[0];
    pAU[1] = vpos_iter->second[1];
    pAU[2] = vpos_iter->second[2];         
    a_points.push_back(pAU);    
  }
                  
}                   

void HelioRealBCProblem::readVoyagersPos()
{
  Real time,lat,lon,vecHGI[3],vec[3]; 
  
  FILE* V_file=fopen("trajV1.dat","r");
  if (V_file==NULL) MayDay::Warning("trajV1.dat not found");
  while (!feof(V_file))
  {
    Real vdata[5];
    
    fscanf(V_file,"%lf %lf %lf %lf %lf", &vdata[0], &vdata[1], &vdata[2], &vdata[3], &vdata[4]);           
           
    time = vdata[0] + (vdata[1]-1)/365.25; 
    
    if (time<m_startBCPhys) continue;
    
    lat  = vdata[3]/180.0*d_PI;
    lon  = vdata[4]/180.0*d_PI;
           
    vecHGI[0] = cos(lat)*cos(lon);
    vecHGI[1] = cos(lat)*sin(lon);
    vecHGI[2] = sin(lat);
    
    vec[0]    = -0.9998430685*vecHGI[0]  + 0.01771548383*vecHGI[1];
    vec[1]    = -0.01771548383*vecHGI[0] - 0.9998430685*vecHGI[1];
    vec[2]    = vecHGI[2];    
    
    lon = acos(vec[0]/(sqrt(vec[0]*vec[0]+vec[1]*vec[1])));   
    if (vec[1]<0.0) lon = d_2PI-lon;
                           
    m_V1pos[time][0] = vdata[2];
    m_V1pos[time][1] = lon;
    m_V1pos[time][2] = d_PI_2-lat;
  }
  fclose(V_file);  
  
  V_file=fopen("trajV2.dat","r");
  if (V_file==NULL) MayDay::Warning("trajV2.dat not found");
  while (!feof(V_file))
  {
    Real vdata[5];
    
    fscanf(V_file,"%lf %lf %lf %lf %lf", &vdata[0], &vdata[1], &vdata[2], &vdata[3], &vdata[4]);           
           
    time = vdata[0] + (vdata[1]-1)/365.25; 
    
    if (time<m_startBCPhys) continue;
           
    lat  = vdata[3]/180.0*d_PI;
    lon  = vdata[4]/180.0*d_PI;
           
    vecHGI[0] = cos(lat)*cos(lon);
    vecHGI[1] = cos(lat)*sin(lon);
    vecHGI[2] = sin(lat);
    
    vec[0]    = -0.9998430685*vecHGI[0]  + 0.01771548383*vecHGI[1];
    vec[1]    = -0.01771548383*vecHGI[0] - 0.9998430685*vecHGI[1];
    vec[2]    = vecHGI[2];    
    
    lon = acos(vec[0]/(sqrt(vec[0]*vec[0]+vec[1]*vec[1])));   
    if (vec[1]<0.0) lon = d_2PI-lon;
                           
    m_V2pos[time][0] = vdata[2];
    m_V2pos[time][1] = lon;
    m_V2pos[time][2] = d_PI_2-lat;
  }
  fclose(V_file);  
  
}

void HelioRealBCProblem::readUlyssesPos()
{
  Real time,lat,lon,vecHGI[3],vec[3],r; 
  
//  FILE* V_file=fopen("ulysses_daily.dat","r");
//  if (V_file==NULL) MayDay::Warning("ulysses_daily.dat not found");
  
  FILE* V_file=fopen("trajU.dat","r");  
  if (V_file==NULL) MayDay::Warning("trajU.dat not found");
  
  while (!feof(V_file))
  {
    //Real vdata[17];
    //fscanf(V_file,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", 
    //       &vdata[0], &vdata[1], &vdata[2], &vdata[3],
    //       &vdata[4], &vdata[5], &vdata[6], &vdata[7],
    //       &vdata[8], &vdata[9], &vdata[10],&vdata[11],
    //       &vdata[12],&vdata[13],&vdata[14],&vdata[15],&vdata[16],&vdata[17]);
    
    Real vdata[5];    
    fscanf(V_file,"%lf %lf %lf %lf %lf", &vdata[0], &vdata[1], &vdata[2], &vdata[3], &vdata[4]);           
           
    time = vdata[0] + (vdata[1]-1)/365.25; 
    
    if (time<m_startBCPhys) continue;
    
    //r    = vdata[3];
    //lat  = vdata[4]/180.0*d_PI;
    //lon  = vdata[5]/180.0*d_PI;
    
    r    = vdata[2];
    lat  = vdata[3]/180.0*d_PI;
    lon  = vdata[4]/180.0*d_PI;
           
    vecHGI[0] = cos(lat)*cos(lon);
    vecHGI[1] = cos(lat)*sin(lon);
    vecHGI[2] = sin(lat);
    
    vec[0]    = -0.9998430685*vecHGI[0]  + 0.01771548383*vecHGI[1];
    vec[1]    = -0.01771548383*vecHGI[0] - 0.9998430685*vecHGI[1];
    vec[2]    = vecHGI[2];    
    
    lon = acos(vec[0]/(sqrt(vec[0]*vec[0]+vec[1]*vec[1])));   
    if (vec[1]<0.0) lon = d_2PI-lon;
                           
    m_Ulyssespos[time][0] = r;
    m_Ulyssespos[time][1] = lon;
    m_Ulyssespos[time][2] = d_PI_2-lat;
  }
  fclose(V_file);  
}

void HelioRealBCProblem::readEarthPos()
{
  Real time,lat,lon,vecHGI[3],vec[3],r; 
    
  FILE* V_file=fopen("trajEarth.dat","r");  
  if (V_file==NULL) MayDay::Warning("trajU.dat not found");
  
  while (!feof(V_file))
  {    
    
    Real vdata[5];   int year; 
    fscanf(V_file,"%i %lf %lf %lf %lf", &year, &vdata[1], &vdata[2], &vdata[3], &vdata[4]);           
           
    int daysYear = (year%4 == 0 ? 366 : 365);
    time = year + (vdata[1]-1)/daysYear; 
    
    if (time<m_startBCPhys) continue;
        
    
    r    = vdata[2];
    lat  = vdata[3]/180.0*d_PI;
    lon  = vdata[4]/180.0*d_PI;
           
    vecHGI[0] = cos(lat)*cos(lon);
    vecHGI[1] = cos(lat)*sin(lon);
    vecHGI[2] = sin(lat);
    
    vec[0]    = -0.9998430685*vecHGI[0]  + 0.01771548383*vecHGI[1];
    vec[1]    = -0.01771548383*vecHGI[0] - 0.9998430685*vecHGI[1];
    vec[2]    = vecHGI[2];    
    
    lon = acos(vec[0]/(sqrt(vec[0]*vec[0]+vec[1]*vec[1])));   
    if (vec[1]<0.0) lon = d_2PI-lon;
                           
    m_Earthpos[time][0] = r;
    m_Earthpos[time][1] = lon;
    m_Earthpos[time][2] = d_PI_2-lat;
  }
  fclose(V_file);  
}
*/




/// Creates "DATASETAUXDATA" fields for tecplot
/**
*/
void HelioRealBCProblem::auxDataTecplot(std::string & a_str,
                                        Real          a_time,
                                        int           a_type)
{    
  Real physTime = getPhysTime(a_time);
  
  //physTime = (physTime - m_startBCPhys)*366.0*24.0;
  
  //if (procID() == 0) pout() << "HelioRealBCProblem::auxDataTecplot physTime: " << physTime << endl;
  
  char buf[50];
  if (a_type == 0)
      sprintf(buf,"DATASETAUXDATA PhysTime = \"%.15g\"",physTime);
  if (a_type == 1)  
    sprintf(buf,       "AUXDATA PhysTime = \"%.15g\"",physTime);
  a_str+=buf;
}

