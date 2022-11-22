#ifndef _CSHCARTESIAN_H_
#define _CSHCARTESIAN_H_

#include "CSHandler.H"
#include <RealVect.H> 

class CSHCartesian : public CoordinateSystemHandler
{
public:
                            
  CSHCartesian();  
  virtual ~CSHCartesian();
  
  void define(eCoordinateSystem  a_CS,
              int                a_max_level,
              const Vector<int>& a_ref_ratios,
              const Box&         a_prob_domain,
              RealVect           a_domainLength,              
              const EquationSystem & a_eqSys);
              
  // Set grid spacing in a_dir direction
  virtual void setGridSpacing(const FArrayBox & a_du,
                                          int   a_dir);
              
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
  /**
    Strictly speaking, this methods should return a coefficent 
    used in updateState to update conservative variables.
    U = U - dt/a_vol * (a_F(a_dir,HiSide) - a_F(a_dir,LoSide))
    
    For cartesian coordinate system this parameter can be equal to dx
  */
  virtual void getCellVolumes(FArrayBox & a_vol,
                              const Box & a_box,
                              int         a_level);
  
  // Get cell centers in a_dir direction in curvilinear coordinate system
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
                           int          a_level);
                        
  // Transform curvilinear coordiantes to cartesian                            
  virtual void transCurvCoordsToCartesian( RealVect & a_cs,
                              const RealVect & a_crv);
                              
  // Transform cartesian coordiantes to curvilinear
  virtual void transCartesianCoordsToCurv( RealVect & a_crv,
                              const RealVect & a_cs);
                              
                              
  virtual void getAreas(FArrayBox & a_areas,
                        const Box & a_box,
                        int         a_dir,
                        int         a_level);
                              
                                       
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
                  
  virtual void scalingFactor(
                        FArrayBox & a_scale, 
                  const Box       & a_box,                        
                        int         a_level);
                        
  // square root of orthogonal components (g_ii) of the metric tensor
  virtual void sqrtMetricTensorCoeff(
                        FArrayBox & a_g,
                  const Box       & a_box,
                  int               a_dir,                  
                  int               a_level);
                  
  
  
protected:

  int m_max_level;
  RealVect*  m_dxdydz;

};

#endif
