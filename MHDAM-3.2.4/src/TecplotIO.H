#ifndef _TECPLOT_IO_H_
#define _TECPLOT_IO_H_

#include <FArrayBox.H> 
#include <Interval.H> 


FILE* OpenTecplotFile(const char* filename,const char* header);
void  WriteFArrayBoxToTecplotFile(FILE* tfile, const FArrayBox& a_fab, const Box& a_box, const Interval& a_comps, Real a_dx);
void  WriteFArrayBoxToTecplotFile(FILE* tfile, const FArrayBox& a_fab, const Box& a_box, const Interval& a_comps, const FArrayBox& a_coords, const FArrayBox& a_coords_cv = FArrayBox());
void  CloseTecplotFile(FILE* tfile);

// Creates standard header for tecplot file.
FILE* OpenTecplotFile(const char* filename,int num_vars);

#endif
