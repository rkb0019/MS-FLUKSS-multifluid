#include "CoarseAverageEdge.H"
#include "AverageEdgeF_F.H"
#include "NamespaceHeader.H"

// ----------------------------------------------------------
CoarseAverageEdge :: CoarseAverageEdge()
: m_isDefined( false )
{
}

CoarseAverageEdge :: CoarseAverageEdge( const DisjointBoxLayout& a_fineGrids,
                                              int                a_nComp,
                                              int                a_nRef       ) 
: m_isDefined( false )
{
  define( a_fineGrids, a_nComp, a_nRef );
}

                                                           // defines the object
void CoarseAverageEdge :: define( const DisjointBoxLayout& a_fineGrids,
			                                  int                a_nComp,
                                        int                a_nRef        )
{
  m_nRef = a_nRef;
  
  DisjointBoxLayout coarsened_fine_domain;
  coarsen( coarsened_fine_domain, a_fineGrids, m_nRef );
  m_coarsenedFineData.define( coarsened_fine_domain, a_nComp );

  m_isDefined = true;
}

bool CoarseAverageEdge :: isDefined() const
{
  return m_isDefined;
}

                                     // averages fine-level data to coarse level
void CoarseAverageEdge :: averageToCoarse(       LevelData<EdgeBox>& a_coarseData,
				                                   const LevelData<EdgeBox>& a_fineData,
                                                 double              a_weight      )
{
  computeAverages( a_coarseData, a_fineData, arithmetic, a_weight );
}

            // averages fine-level data to coarse level using harmonic averaging
void CoarseAverageEdge :: averageToCoarseHarmonic(       LevelData<EdgeBox>& a_coarseData,
                                                   const LevelData<EdgeBox>& a_fineData,
                                                         double              a_weight      )
{
  computeAverages( a_coarseData, a_fineData, harmonic, a_weight );
}

    // utility function called by both averageToCoarse and averageCoarseHarmonic
void CoarseAverageEdge :: computeAverages(       LevelData<EdgeBox>& a_coarseData, 
                                           const LevelData<EdgeBox>& a_fineData,
                                                 int                 a_averageType,
                                                 double              a_weight      )
{
  CH_assert(isDefined());

  DataIterator dit = a_fineData.dataIterator();
  for( dit.reset(); dit.ok(); ++dit ) 
  {
    EdgeBox& coarsenedFine = m_coarsenedFineData[dit()];
    const EdgeBox& fine = a_fineData[dit()];

                // coarsen from the entire fine grid onto the entire coarse grid
    averageGridData( coarsenedFine, fine, a_averageType, a_weight );
  }

  m_coarsenedFineData.copyTo( m_coarsenedFineData.interval(),
                              a_coarseData, a_coarseData.interval() );
}

                           // averages entire single grid data from fine->coarse
void CoarseAverageEdge :: averageGridData(       EdgeBox& a_coarsenedFine, 
                                           const EdgeBox& a_fine,
                                                 int      a_averageType,
                                                 double   a_weight        ) const
{
  if( SpaceDim ==  3 )
  {
    for( int dir = 0; dir < SpaceDim; dir++ )
    {
      FArrayBox& coarseFab = a_coarsenedFine[dir];
      const FArrayBox& fineFab = a_fine[dir];

      const Box& coarseBox = coarseFab.box();

                                                        // set up refinement box
      IntVect hiVect(D_DECL(0,0,0));
      hiVect.setVal( dir, m_nRef - 1 );
      IntVect loVect(D_DECL(0,0,0));
      Box refBox(loVect, hiVect);    

      if( a_averageType == arithmetic )
      {
        FORT_AVERAGEEDGE( CHF_FRA(coarseFab),
                          CHF_CONST_FRA(fineFab),
                          CHF_BOX(coarseBox),
                          CHF_CONST_INT(m_nRef),
                          CHF_CONST_REAL(a_weight),
                          CHF_BOX(refBox)          );
      }
      else if( a_averageType == harmonic )
      {
        FORT_AVERAGEEDGEHARMONIC( CHF_FRA(coarseFab),
                                  CHF_CONST_FRA(fineFab),
                                  CHF_BOX(coarseBox),
                                  CHF_CONST_INT(m_nRef),
                                  CHF_CONST_REAL(a_weight),
                                  CHF_BOX(refBox)         );
      }
      else 
      {
        MayDay::Error("CoarseAverageEdge::averageGridData -- bad averageType");
      }
    }
  }
  else
  {
    FArrayBox& coarseFab = a_coarsenedFine[0];
    const FArrayBox& fineFab = a_fine[0];

    const Box& coarseBox = coarseFab.box();

    FORT_FINETOCOARSE( CHF_FRA(coarseFab),
                       CHF_CONST_FRA(fineFab),
                       CHF_BOX(coarseBox),
                       CHF_CONST_REAL(a_weight),
                       CHF_CONST_INT(m_nRef)   );
  }
}
  
#include "NamespaceFooter.H"
