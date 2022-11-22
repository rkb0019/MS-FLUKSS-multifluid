#ifndef _CoarseAverageFace_H_
#define _CoarseAverageFace_H_

#include "REAL.H"
#include "EdgeBox.H"
#include "LevelData.H"
#include "NamespaceHeader.H"


/// replaces edge-centered coarse-level data w/ averaged fine-level data
/** This class replaces edge-centered data on a coarse level of refinement
    with the average of the finer-level data which shares the edge.
    This class is slightly modified version of CoarseAverageFace */
class CoarseAverageEdge
{
public:
  /// Default constructor
  CoarseAverageEdge();

  /// defining constructor
  CoarseAverageEdge( const DisjointBoxLayout& a_fineGrids, 
		                       int                a_nComp,
                           int                a_nRef       );

  /// destructor
  ~CoarseAverageEdge(){};

  /// defines the object
  void define( const DisjointBoxLayout& a_fineGrids, 
	                   int                a_nComp,
                     int                a_nRef       );
  
  ///
  bool isDefined() const;

  /// averages fine-level data to coarse level
  void averageToCoarse(       LevelData<EdgeBox>& a_coarse_data, 
		                    const LevelData<EdgeBox>& a_fine_data,
                              double              a_weight = 1.0 );

  /// averages fine-level data to coarse level using harmonic averaging
  void averageToCoarseHarmonic(       LevelData<EdgeBox>& a_coarse_data, 
                                const LevelData<EdgeBox>& a_fine_data,
                                      double              a_weight = 1.0 );

protected:

  ///
  enum averageType{ arithmetic = 0,
                    harmonic, 
                    NUM_AVERAGE_TYPES
                  };

  bool m_isDefined;

  int m_nRef;

  // work array for coarsening of fine data, same "shape" as fine data
  LevelData<EdgeBox> m_coarsenedFineData;

  /** utility function called by both averageToCoarse and 
      averageCoarseHarmonic (to avoid code duplication)
  */
  void computeAverages(       LevelData<EdgeBox>& a_coarse_data, 
		                    const LevelData<EdgeBox>& a_fine_data,
                              int                 a_averageType,
                              double              a_weight      );

  /// averages entire single grid data from fine->coarse
  void averageGridData(       EdgeBox& a_coarsenedFine, 
                        const EdgeBox& fine,
                              int      a_averageType,
                              double   a_weight      ) const;
};

#include "NamespaceFooter.H"
#endif
  



