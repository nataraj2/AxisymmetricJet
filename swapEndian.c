/* 
 * swap the endian-ness of an unformatted PLOT3D file
 *
 * initial version Monday, July 16, 2007  7:52:40 -0500
 *
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define S_PREC 4
#define D_PREC 8
#define MAXCHAR 256

void endianSwap(long int, int, int);
size_t safe_fread(long int, int, int, FILE *);
size_t safe_fwrite(long int, int, int, FILE *);
void swapGridFile(char *, char *, char *, char *, char *, char *, char *);
// void swapSolnFile(char *, char *, char *, char *, char *, char *, char *);

/* global integers */
int SWAP_READ = 0;
int SWAP_WRITE = 0;

int main(int argc, char *argv[])
{

  FILE *in;
  char prec[2], ib[2], vf[2], ftype[2], gf[2], filename_in[MAXCHAR];
  char filename_out[MAXCHAR];

  if (argc != 8) {
    fprintf(stderr, "ERROR: Incorrect command line.\n");
    fprintf(stderr, "Usage: %s <precision> <grid format> <volume format> <iblank> <file type> <in file> <out file>\n",argv[0]);
    exit(EXIT_FAILURE);
  }

  /* read command line */
  strncpy(prec,argv[1],1); prec[1] = '\0';
  strncpy(gf,argv[2],1); gf[1] = '\0';
  strncpy(vf,argv[3],1); vf[1] = '\0';
  strncpy(ib,argv[4],1); ib[1] = '\0';
  strncpy(ftype,argv[5],1); ftype[1] = '\0';
  strncpy(filename_in,argv[6],MAXCHAR-1);
  strncpy(filename_out,argv[7],MAXCHAR-1);

  /* print command line */
  fprintf(stdout,"precision: %s\n", prec);
  fprintf(stdout,"grid format: %s\n", gf);
  fprintf(stdout,"volume format: %s\n", vf);
  fprintf(stdout,"iblank present: %s\n", ib);
  fprintf(stdout,"file type: %s\n", ftype);
  fprintf(stdout,"in file; %s\n", filename_in);
  fprintf(stdout,"out file: %s\n", filename_out);
 
  /* query user for swap order */
  fprintf(stdout,">>>> Do we swap on read? (0=no, 1=yes): ");
  fscanf(stdin,"%d",&SWAP_READ);
  if (SWAP_READ == 0) {
    fprintf(stdout,"Swapping on write.\n");
    SWAP_WRITE = 1;
  } else {
    fprintf(stdout,"Swapping on read.\n");
  }

  if ((in = fopen(filename_in,"rb+")) == NULL) {
    fprintf(stdout,"Unable to open file \"%s\".", filename_in);
    exit(EXIT_FAILURE);
  }

  if (strcmp(ftype,"g") == 0) 
     swapGridFile(prec,gf,vf,ib,ftype,filename_in,filename_out);
  // if (strcmp(ftype,"s") == 0) swapSolnFile(prec,gf,vf,ib,ftype,filename_in,filename_out);
 
  return EXIT_SUCCESS;

}

void swapGridFile(char *prec, char *gf, char *vf, char *ib, char *ftype, char *fname_in, char *fname_out)
{
  int record, nzones, i, ND, *N;
  FILE *in, *out;
  double *dX;
  float *fX;

  in = fopen(fname_in,"rb");
  out = fopen(fname_out,"wb");

  /* number of zones, if present */
  nzones = 1;
  if (strcmp(gf,"m") == 0) {
    safe_fread((long int)&record,sizeof(int),1,in);
    safe_fread((long int)&nzones,sizeof(int),1,in);
    safe_fread((long int)&record,sizeof(int),1,in);

    safe_fwrite((long int)&record,sizeof(int),1,out);
    safe_fwrite((long int)&nzones,sizeof(int),1,out);
    safe_fwrite((long int)&record,sizeof(int),1,out);
  }
  fprintf(stdout,"\"%s\" has %d zone(s)\n", fname_in, nzones);

  /* query user for number of dimensions */ 
  fprintf(stdout,">>>> Input number of dimensions: ");
  fscanf(stdin,"%d", &ND);

  /* read the zonal dimensions */
  safe_fread((long int)&record,sizeof(int),1,in);
  N = (int *)malloc((size_t)record);
  safe_fread((long int)N, sizeof(int), (size_t)(record/sizeof(int)), in);
  safe_fread((long int)&record,sizeof(int),1,in);

  safe_fwrite((long int)&record,sizeof(int),1,out);
  safe_fwrite((long int)N, sizeof(int), (size_t)(record/sizeof(int)), out);
  safe_fwrite((long int)&record,sizeof(int),1,out);

  /* read the grid */
  if (strcmp(vf,"w") == 0) {

    for (i = 0; i < nzones; i++) {
      safe_fread((long int)&record,sizeof(int),1,in);
      if (strcmp(prec,"d") == 0) {
        dX = (double *)malloc((size_t)record);
        safe_fread((long int)dX,sizeof(double),(size_t)(record/sizeof(double)),in);
      } else {
        fX = (float *)malloc((size_t)record);
        safe_fread((long int)fX,sizeof(float),(size_t)(record/sizeof(float)),in);
      }
      safe_fread((long int)&record,sizeof(int),1,in);

      safe_fwrite((long int)&record,sizeof(int),1,out);
      if (strcmp(prec,"d") == 0) {
        safe_fwrite((long int)dX,sizeof(double),(size_t)(record/sizeof(double)),out);
        free(dX);
      } else {
        safe_fwrite((long int)fX,sizeof(float),(size_t)(record/sizeof(float)),out);
        free(fX);
      }
      safe_fwrite((long int)&record,sizeof(int),1,out);
    }

  } else {
 
    fprintf(stdout,"ERROR: grid file, planes format not supported");
    exit(EXIT_FAILURE);

  }
    
  fclose(in);
  fclose(out);  

  return;

}


size_t safe_fread(long int addr, int size_el, int num_el, FILE *in)
{
  size_t retval;
  char *a;

  a = (char *)addr;
  
  retval = fread(a, size_el, num_el, in);
 
  if (SWAP_READ == 1) endianSwap((long int)a, num_el, size_el);

  return(retval);
}

size_t safe_fwrite(long int addr, int size_el, int num_el, FILE *in)
{
  size_t retval;
  char *a;

  a = (char *)addr;

  /* swap coming in */ 
  if (SWAP_WRITE == 1) endianSwap((long int)a, num_el, size_el);

  retval = fwrite(a, size_el, num_el, in);

  /* swap again to keep data readable */
  if (SWAP_WRITE == 1) endianSwap((long int)a, num_el, size_el);

  return(retval);
}

/*********************************************************************
 *
 * endianSwap(addr,num_el,size_el): Changes between big- and 
 * little-endian by address-swapping.
 *
 * addr: (long int) address of variable to be swapped.
 * num_el: (int) number of elements in array pointed to by addr
 * size_el: (int) size of individual elements of *addr
 *
 * Written by Daniel J. Bodony (bodony@Stanford.EDU)
 * Copyright (c) 2001
 *
 * WARNING: This is only works on arrays that contiguous in memory.
 *
 *********************************************************************/

#define SWAP(a,b) temp=(a); (a)=(b); (b)=temp;

void endianSwap (long int addr, int num_el, int size_el) 
{
   int i, x;
   char *a;
   char temp;

   a = (char *)addr;

   for (i = 0; i < num_el; i++) {
      a = (char *)(addr + i*size_el);
      for (x = 0; x < size_el/2; x++) {
         SWAP(a[x],a[size_el-x-1]);
      }
   }
}
