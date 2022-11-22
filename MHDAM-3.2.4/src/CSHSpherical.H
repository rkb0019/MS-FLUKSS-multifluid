#ifndef _CSHSPHERICAL_H_
#define _CSHSPHERICAL_H_

#include "CSHandler.H"
#include <RealVect.H> 

class CSHSpherical : public CoordinateSystemHandler
{
public:
    enum eRspacing {
                        R_Const        = 0,
                        R_UniformCells = 1,
                        R_Exp          = 2
                   };

                            
  CSHSpherical();  
  virtual ~CSHSpherical();
  
  void define(eCoordinateSystem   a_CS,
              int                 a_max_level,
              const Vector<int>&  a_ref_ratios,
              const Box&          a_prob_domain,
              RealVect            a_domainLength,
              const Vector<Real>& a_domainBox,
              const EquationSystem & a_eqSys
             );
             
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );
  
  // Set grid spacing in a_dir direction
  virtual void setGridSpacing(const FArrayBox & a_du,
                                          int   a_dir);
    
    // Returns true if grid spacing is constant in this direction  
  virtual bool constStep(int a_dir);
  
  // Returns grid spacing if it is constant in this direction  
  virtual Real dx(int a_dir, int a_level);
    
  
  // Returns grid spacing. Works for all types of directions.
  /**
    a_box can contain cells outside of the problem domain.
    a_dx must be virtually 1D array in a_dir direction, 
    i.e. for all dirs except of a_dir the follwing should be fulfilled
         a_dx.box().smallEnd(dir) == 0
         a_dx.box().bigEnd(dir) == 0  
  */
  virtual void dx(FArrayBox & a_dx,
                  const Box & a_box,
                  int         a_dir,
                  int         a_level);
  
  // Returns cell volumes  
  virtual void getCellVolumes(FArrayBox & a_vol,
                              const Box & a_box,
                              int         a_level);
  
  // Get cell centers in a_dir direction
  /**
    a_box can contain cells outside of the problem domain.
    a_dx must be virtually 1D array in a_dir direction, 
    i.e. for all dirs except of a_dir the follwing should be fulfilled
         a_dx.box().smallEnd(dir) == 0
         a_dx.box().bigEnd(dir) == 0  
  */
  virtual void getCellCenters(FArrayBox & a_centers,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level);
                              
  virtual void getCellCenter( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level);
  
  // Get curvilinear coordiantes of nodes
  virtual void getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level);
                              
  // Get curvilinear coordiantes of nodes
  virtual void getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level);
  
  // Get cartesian coordiantes of nodes
  virtual void getNodeCoordsCartesian( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level);
                                
                              
  // Get cartesian coordiantes of a node
  virtual void getNodeCoordsCartesian( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level);
                              
  // Get the closest cell for a point given in curvilinear coordinates
  virtual void getClosestCell(IntVect  & a_iv,
                        const RealVect & a_coords,
                        int              a_level);
                        
  // Get the closest index for a given curvilinear coordinate
  virtual void getClosestIndex(int  & a_i,
                           const Real & a_coord,
                           int          a_dir,                           
                           int          a_level) ;
                              
  // Transform curvilinear coordiantes to cartesian
  virtual void transCurvCoordsToCartesian( RealVect & a_cs,
                              const RealVect & a_crv);
                              
  // Transform cartesian coordiantes to curvilinear
  virtual void transCartesianCoordsToCurv( RealVect & a_crv,
                              const RealVect & a_cs);
  
  // Get cells areas
  virtual void getAreas(FArrayBox & a_areas,
                        const Box & a_box,
                        int         a_dir,
                        int         a_level);
                        
                        
  // Modifies conservative variables in a_interp to  values
  // that will be used PiecewiseLinearFillPatchMHDAM for filling ghost cells
  virtual void PLFP_ConsToInterpVars(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level);
    
  // Reverse operation to PLFP_ConsToInterpVars
  virtual void PLFP_InterpVarsToCons(
                        FArrayBox & a_interp,                   
                  const Box       & a_box,
                  int               a_level);
                  
                  
  // Transform Cartesian components of a vector into curvilinear
  virtual void transCartesianVectToCurv(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level);
                  
  // Transform curvilinear components of a vector into Cartesian
  virtual void transCurvVectToCartesian(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level);
                  
  // Transform fluxes calculated by Riemann problem class for using them in UpdateState
  virtual void transFluxesForUpdateState(                        
                        FArrayBox & a_F,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level);
                  
  virtual void scalingFactor(
                        FArrayBox & a_scale, 
                  const Box       & a_box,                        
                        int         a_level);
                        
   
  // Retutn whether fluxes should be scaled in LevelFluxRegister Class (parameter m_scaleFineFluxes)
  virtual bool scaleFineFluxes();
  
  // Calculates fluxes for using them in FluxRegistry.incrementCoarse
  virtual void modifyFluxesFRincrementCoarse(
                        FArrayBox & fluxC, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level);

  // Calculates fluxes for using them in FluxRegistry.incrementFine
  virtual void modifyFluxesFRincrementFine(
                        FArrayBox & fluxF, 
                  const FArrayBox & flux, 
                        int         idir,                  
                  int               a_level);
                  
  // square root of orthogonal components (g_ii) of the metric tensor
  virtual void sqrtMetricTensorCoeff(
                        FArrayBox & a_g,
                  const Box       & a_box,
                  int               a_dir,
                  int               a_level);

#ifdef CH_USE_HDF5
  // Reads coordinate specific data  
  virtual void readGeomInfo(hid_t a_fileID, bool a_dirs[CH_SPACEDIM]);

  // Writes coordinate specific data
  virtual void writeGeomInfo(hid_t a_fileID);
#endif

  // Boundary conditions through the z-axis. 
  /** We assume that a_U ghost cells are already filled with the conservative variables
      shifted on PI in phi direction and all vectors are transformed to curvilinear coordinate system.    */
  void zaxisBC(        FArrayBox  & a_W,
                 const Box        & a_box,
                       int          a_sign);
  
   
protected:

  struct SphericalLevel
  {
    Real m_dudvdw_const[CH_SPACEDIM]; // dr, dphi, dtheta if constStep[dir] == true
    
    FArrayBox m_dudvdw[CH_SPACEDIM];
    FArrayBox m_uvwc[CH_SPACEDIM];
    FArrayBox m_uvwn[CH_SPACEDIM];
    
    Box m_prob_domain;
            
  };
  SphericalLevel * m_gridData;
  
  bool m_constStep[CH_SPACEDIM];
      

  // By default we provide possibility to have nonuniform mesh in 'r' direction.
  // if other direction should be non-uniform use PhysProblem::defineMesh()
  eRspacing m_Rspacing;  
  Real      m_expAlpha;
  
};

#endif
