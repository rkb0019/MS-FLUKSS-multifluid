#include "TecplotIO.H"
#include <string>

FILE* extern_tecplotfile1;

FILE* OpenTecplotFile(const char* filename,const char* header)
{
  FILE* plot_file=fopen(filename,"w");  
  fprintf(plot_file,"%s\n",header);
  return plot_file;
}

// Creates standard header for tecplot file.
FILE* OpenTecplotFile(const char* filename,int num_vars)
{
  FILE* plot_file=fopen(filename,"w");  
  if (CH_SPACEDIM == 2) fprintf(plot_file,"VARIABLES=\"X\" \"Y\"");
  if (CH_SPACEDIM == 3) fprintf(plot_file,"VARIABLES=\"X\" \"Y\" \"Z\"");
  
  int i; std::string variables_string;char buffer[100];  
  for (i=0;i<num_vars;i++)
  {
    sprintf(buffer," \"var%i\"",i);
    variables_string+=buffer;
  }
  fprintf(plot_file,"%s\n",variables_string.c_str());
  
  return plot_file;
}


// Attn! This function has never been tested in 3D.
void  WriteFArrayBoxToTecplotFile(FILE* tfile, const FArrayBox& a_fab, const Box& a_box, const Interval& a_comps, Real a_dx)
{
  int D_DECL(i,j,k);    
  int iComp;  
  Real D_DECL(*x_coords=NULL,*y_coords=NULL,*z_coords=NULL);
  
  int max_numbers_in_row = 30;
  
  int written_numbers = 0;
  

  
  D_TERM(  
    int nx=a_box.size(0);,
    int ny=a_box.size(1);,
    int nz=a_box.size(2););
    
  
  IntVect LCorner=a_box.smallEnd();
  IntVect UCorner=a_box.bigEnd();
  
  D_TERM(  
    CH_assert(nx==(UCorner[0]-LCorner[0]+1));,
    CH_assert(ny==(UCorner[1]-LCorner[1]+1));,
    CH_assert(nz==(UCorner[2]-LCorner[2]+1)););
  
  D_TERM(  
    x_coords = new Real[nx+1];,
    y_coords = new Real[ny+1];,
    z_coords = new Real[nz+1];);
  
  D_TERM(  
  for (i=LCorner[0];i<=UCorner[0]+1;i++)  x_coords[i-LCorner[0]] = i*a_dx;,
  for (j=LCorner[1];j<=UCorner[1]+1;j++)  y_coords[j-LCorner[1]] = j*a_dx;,    
  for (k=LCorner[2];k<=UCorner[2]+1;k++)  z_coords[k-LCorner[2]] = k*a_dx;);

#if CH_SPACEDIM == 2   
  fprintf(tfile,"ZONE I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([3-%i]=CELLCENTERED)\n", nx+1, ny+1, a_comps.size()+CH_SPACEDIM);
#endif  
#if CH_SPACEDIM == 3   
  fprintf(tfile,"ZONE I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([4-%i]=CELLCENTERED)\n", nx+1, ny+1, nz+1, a_comps.size()+CH_SPACEDIM);
#endif
  
  written_numbers = 0;
  
#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
#endif
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++) 
  {
    fprintf(tfile,"%.6e ",x_coords[i-LCorner[0]]);
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }
  }

#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
#endif  
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++)
  {
    fprintf(tfile,"%.6e ",y_coords[j-LCorner[1]]);    
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }    
  }

#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++)
  {
    fprintf(tfile,"%.6e ",z_coords[k-LCorner[2]]);    
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }    
  }
#endif    
  
  IntVect aux_IntVect; Real value;
  
  for (iComp = a_comps.begin();iComp <= a_comps.end();++iComp)
#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2];k++)
#endif  
  for (j=LCorner[1];j<=UCorner[1];j++)               
  for (i=LCorner[0];i<=UCorner[0];i++)
  { 
    D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
    value=a_fab.get(aux_IntVect,iComp);      
    fprintf(tfile,"%.6e ",value);        
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }
  }          

    
  D_TERM(
    delete[] x_coords;,
    delete[] y_coords;,  
    delete[] z_coords;);  
}


void  WriteFArrayBoxToTecplotFile(FILE* tfile, const FArrayBox& a_fab, const Box& a_box, const Interval& a_comps, 
                                              const FArrayBox& a_coords, const FArrayBox& a_coords_cv )
{
  int D_DECL(i,j,k);
  Real D_DECL(x,y,z);
    
  int iComp;      
  int max_numbers_in_row = 30;  
  int written_numbers = 0;
  
  D_TERM(  
    int nx=a_box.size(0);,
    int ny=a_box.size(1);,
    int nz=a_box.size(2););
    
  
  IntVect LCorner=a_box.smallEnd();
  IntVect UCorner=a_box.bigEnd();
  
  D_TERM(  
    CH_assert(nx==(UCorner[0]-LCorner[0]+1));,
    CH_assert(ny==(UCorner[1]-LCorner[1]+1));,
    CH_assert(nz==(UCorner[2]-LCorner[2]+1)););
    
#if CH_SPACEDIM == 2   
  if (a_coords_cv.nComp() == 0)
    fprintf(tfile,"ZONE I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([3-%i]=CELLCENTERED)\n", nx+1, ny+1, a_comps.size()+CH_SPACEDIM);    
  else
    fprintf(tfile,"ZONE I=%i,J=%i, DATAPACKING=BLOCK, VARLOCATION=([5-%i]=CELLCENTERED)\n", nx+1, ny+1, a_comps.size()+2*CH_SPACEDIM);
    
#endif  
#if CH_SPACEDIM == 3   
  if (a_coords_cv.nComp() == 0)
    fprintf(tfile,"ZONE I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([4-%i]=CELLCENTERED)\n", nx+1, ny+1, nz+1, a_comps.size()+CH_SPACEDIM);
  else
    fprintf(tfile,"ZONE I=%i,J=%i,K=%i, DATAPACKING=BLOCK, VARLOCATION=([7-%i]=CELLCENTERED)\n", nx+1, ny+1, nz+1, a_comps.size()+2*CH_SPACEDIM);
#endif
  
  written_numbers = 0;
  
#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
#endif
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++) 
  {
    x = a_coords.get(IntVect(D_DECL(i,j,k)),0); 
    fprintf(tfile,"%.6e ",x);
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }
  }

#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
#endif  
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++)
  {
    y = a_coords.get(IntVect(D_DECL(i,j,k)),1); 
    fprintf(tfile,"%.6e ",y);
    
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }    
  }

#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2]+1;k++)
  for (j=LCorner[1];j<=UCorner[1]+1;j++)
  for (i=LCorner[0];i<=UCorner[0]+1;i++)
  {
    z = a_coords.get(IntVect(D_DECL(i,j,k)),2); 
    fprintf(tfile,"%.6e ",z);
        
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }    
  }
#endif    

  // Curvilinear coordinates
  if (a_coords_cv.nComp() > 0)
  {
  
  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    for (i=LCorner[0];i<=UCorner[0]+1;i++) 
    {
      x = a_coords_cv.get(IntVect(D_DECL(i,j,k)),0); 
      fprintf(tfile,"%.6e ",x);
      written_numbers++;
      if (written_numbers > max_numbers_in_row)
      {
        fprintf(tfile,"\n");
        written_numbers = 0;
      }
    }

  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
  #endif  
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    for (i=LCorner[0];i<=UCorner[0]+1;i++)
    {
      y = a_coords_cv.get(IntVect(D_DECL(i,j,k)),1); 
      fprintf(tfile,"%.6e ",y);
      
      written_numbers++;
      if (written_numbers > max_numbers_in_row)
      {
        fprintf(tfile,"\n");
        written_numbers = 0;
      }    
    }

  #if CH_SPACEDIM == 3 
    for (k=LCorner[2];k<=UCorner[2]+1;k++)
    for (j=LCorner[1];j<=UCorner[1]+1;j++)
    for (i=LCorner[0];i<=UCorner[0]+1;i++)
    {
      z = a_coords_cv.get(IntVect(D_DECL(i,j,k)),2); 
      fprintf(tfile,"%.6e ",z);
          
      written_numbers++;
      if (written_numbers > max_numbers_in_row)
      {
        fprintf(tfile,"\n");
        written_numbers = 0;
      }    
    }
  #endif    
  }

  
  IntVect aux_IntVect; Real value;
  
  for (iComp = a_comps.begin();iComp <= a_comps.end();++iComp)
#if CH_SPACEDIM == 3 
  for (k=LCorner[2];k<=UCorner[2];k++)
#endif  
  for (j=LCorner[1];j<=UCorner[1];j++)               
  for (i=LCorner[0];i<=UCorner[0];i++)
  { 
    D_TERM(aux_IntVect[0]=i;,aux_IntVect[1]=j;,aux_IntVect[2]=k;);
    value=a_fab.get(aux_IntVect,iComp);      
    fprintf(tfile,"%.6e ",value);        
    written_numbers++;
    if (written_numbers > max_numbers_in_row)
    {
      fprintf(tfile,"\n");
      written_numbers = 0;
    }
  }          
    
  
}


void  CloseTecplotFile(FILE* tfile)
{
  if (tfile!=NULL) fclose(tfile);
}
