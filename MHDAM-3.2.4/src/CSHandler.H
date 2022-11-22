#ifndef _COORDINATESYSTEMHANDLER_H_
#define _COORDINATESYSTEMHANDLER_H_

#include <FArrayBox.H>
#include <Vector.H> 
#include <Box.H> 
#include <RealVect.H> 
#include <NodeFArrayBox.H> 
#include <ParmParse.H>
#include <CH_HDF5.H>


class EquationSystem;

class CoordinateSystemHandler
{
public:
  enum eCoordinateSystem {  
                        CS_Cartesian    = 0, // Cartesian 2D, 3D
                        CS_Axisymmetric = 1, // Cartesian axially-symmetrical, vers 1.
                        CS_Cylindrical  = 2, // Cartesian axially-symmetrical, vers 2.
                        CS_Polar        = 3, // Polar grid
                        CS_PolarAxisym  = 4, // Polar axially-symmetrical, vers 1.
                        CS_Spherical    = 5  // 3D spherical
                        };
                        
  
  
  CoordinateSystemHandler();  
  virtual ~CoordinateSystemHandler();
  
  void define(eCoordinateSystem  a_CS,
              int                a_max_level,
              const Vector<int>& a_ref_ratios,
              const Box&         a_prob_domain,
              const EquationSystem & a_eqSys);
              
  /// Input parameters
  /**
   */
  virtual void input( ParmParse & parser, int verbosity );


  
  eCoordinateSystem coordinateSystem( void ) const;    
                        
  // Set grid spacing in a_dir direction
  virtual void setGridSpacing(const FArrayBox & a_du,
                                          int   a_dir) = 0;
  
  // Returns true if grid spacing is constant in this direction  
  virtual bool constStep(int a_dir) = 0;
  
  // Returns grid spacing if it is constant in this direction  
  virtual Real dx(int a_dir, int a_level) = 0;
  
  const RealVect &  smallEndCurv () const;  
  const RealVect &  bigEndCurv () const;

    
  
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
                  int         a_level) = 0;
  
  // Returns cell volumes  
  virtual void getCellVolumes(FArrayBox & a_vol,
                              const Box & a_box,
                              int         a_level) = 0;
  
  // Get cell centers in a_dir direction in curvilinear coordinate system
  /**
    a_box can contain cells outside of the problem domain.
    a_centers must be virtually 1D array in a_dir direction, 
    i.e. for all dirs except of a_dir the follwing should be fulfilled
         a_dx.box().smallEnd(dir) == 0
         a_dx.box().bigEnd(dir) == 0  
  */
  virtual void getCellCenters(FArrayBox & a_centers,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level) = 0;
  
  // Get curvilinear coordiantes of a cell center
  virtual void getCellCenter( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level) = 0;
  
  // Get curvilinear coordiantes of nodes
  virtual void getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,
                              int         a_dir,
                              int         a_level) = 0;
                              
  // Get curvilinear coordiantes of nodes
  virtual void getNodeCoords( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level) = 0;
                              
  // Get curvilinear coordiantes of nodes
  virtual void getNodeCoords( NodeFArrayBox & a_coords,                              
                              int         a_level);
  
  // Get cartesian coordiantes of nodes
  virtual void getNodeCoordsCartesian( FArrayBox & a_coords,
                              const Box & a_box,                              
                              int         a_level) = 0;
  
  // Get cartesian coordiantes of nodes
  virtual void getNodeCoordsCartesian( NodeFArrayBox & a_coords,                              
                              int         a_level);
  
  // Get cartesian coordiantes of a node
  virtual void getNodeCoordsCartesian( RealVect & a_coords,
                              const IntVect & a_iv,                              
                              int             a_level) = 0;
  
  // Get the closest cell for a point given in curvilinear coordinates
  virtual void getClosestCell(IntVect  & a_iv,
                        const RealVect & a_coords,
                        int              a_level) = 0;
                        
  // Get the closest index for a given curvilinear coordinate
  virtual void getClosestIndex(int  & a_i,
                           const Real & a_coord,
                           int          a_dir,                           
                           int          a_level) = 0;
                              
  // Transform curvilinear coordiantes to cartesian
  virtual void transCurvCoordsToCartesian( RealVect & a_cs,
                              const RealVect & a_crv) = 0;
                              
  // Transform cartesian coordiantes to curvilinear
  virtual void transCartesianCoordsToCurv( RealVect & a_crv,
                              const RealVect & a_cs) = 0;
  
  // Get cells areas. This is not physical areas but coefficients fluxes are multiplyed by.
  virtual void getAreas(FArrayBox & a_areas,
                        const Box & a_box,
                        int         a_dir,
                        int         a_level) = 0;
                              
  
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
                           
  // Transform Cartesian components of vectors into curvilinear
  // a_vectors contains first components of vectors needed to be transformed
  virtual void transCartesianVectToCurv(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level);
                  
  // Transform Cartesian components of vectors into curvilinear  
  virtual void transCartesianVectToCurv(
                        FArrayBox & a_U,                   
                  const Box       & a_box,
                  int               a_level);
                  
  // Transform curvilinear components of vectors into Cartesian
  // a_vectors contains first components of vectors needed to be transformed
  virtual void transCurvVectToCartesian(
                        FArrayBox & a_U, 
                  const int       * a_vectors,
                        int         a_numvectors,
                  const Box       & a_box,
                  int               a_level);
                  
  // Transform curvilinear components of vectors into Cartesian  
  virtual void transCurvVectToCartesian(
                        FArrayBox & a_U,                   
                  const Box       & a_box,
                  int               a_level);
                    
  
  virtual void multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level);                  
                  
  virtual void multiplyFluxesByArea(                        
                        FArrayBox & a_F,
                  const Interval  & a_comps,
                  const FArrayBox & a_areas,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level);                  
  
  // Transform fluxes calculated by Riemann problem class for using them in UpdateState
  virtual void transFluxesForUpdateState(                        
                        FArrayBox & a_F,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level);                  
    
  
  /**
    This methods returns a coefficent 
    used in updateState to update conservative variables.
    U = U - dt * scale * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
    
    Usually scale is a inversion of cell volumes
    For cartesian coordinate system this parameter can be equal to dx
    
    The same coefficient is used for refluxing
  */
  virtual void scalingFactor(
                        FArrayBox & a_scale, 
                  const Box       & a_box,                        
                        int         a_level) = 0;
  
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
                  int               a_level) = 0;
  
  
#ifdef CH_USE_HDF5
  // Reads coordinate specific data  
  virtual void readGeomInfo(hid_t a_fileID);

  // Reads coordinate specific data  
  virtual void readGeomInfo(hid_t a_fileID, bool a_dirs[CH_SPACEDIM]);

  // Writes coordinate specific data
  virtual void writeGeomInfo(hid_t a_fileID);
#endif

protected:
  eCoordinateSystem m_CoordinateSystem;
  
  Vector<int> m_ref_ratios;
  
  int m_max_level;
  
  // Problem domain of the base level
  Box         m_prob_domain;
  
  RealVect m_loCurv;
  RealVect m_hiCurv;  
  
  // First components of vectors in the array of conservative/primitive variables
  int *m_vectors;
  
  // Number of vectors 
  int m_nVectors;    
  
  // First components of velocity vectors in the array of conservative/primitive variables
  int *m_velVectors;
  
  // Number of velocity vectors 
  int m_nVelVectors;    
  
  int m_verbosity;



};

#endif
