/* 
 * Automatically detect the type of PLOT3D file to be read
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define GRID_FILE 0
#define SOLN_FILE 1

void endianSwap(long int, int, int);
size_t safe_fread(long int, int, int, FILE *);

void p3d_detect_(int *len, char *fname, char *prec, int *ND, char *gf, char *vf, char *ib, int *ftype)
{

  FILE *in;
  int i, record, nzones, *N, npts, el_size;

  fname[*(len)*sizeof(char)] = '\0';

  if ((in = fopen(fname,"rb")) == NULL)
    perror((char *)sprintf("Unable to open file %s",fname));

  /* read first record: if record == sizeof(int): multizone, else single zone */
  safe_fread((long int)&record,sizeof(int),1,in);

  /* if we have a multi-zone file */
  if (record == sizeof(int)) {
    sprintf(gf,"%s","m");
    safe_fread((long int)&nzones,sizeof(int),1,in);
    safe_fread((long int)&record,sizeof(int),1,in);
    safe_fread((long int)&record,sizeof(int),1,in);

  } else { /* we have a single-zone file */    
    sprintf(gf,"%s","s");
    nzones = 1;
  }

  /* check number of dimensions: ND = record / (nzones * sizeof(int)) */
  *(ND) = record / (nzones * sizeof(int));

  printf("%d %d %d", record, nzones, sizeof(int));

  /* read the zonal dimensions */
  N = (int *)malloc((size_t)record);
  safe_fread((long int)N, sizeof(int), (size_t)(record/sizeof(int)), in);
  safe_fread((long int)&record,sizeof(int),1,in);

  /* are we a solution file or a grid file */
  safe_fread((long int)&record,sizeof(int),1,in);

  if ((record == (4 * sizeof(float))) || (record == (4 * sizeof(double)))) {
    *(ftype) = SOLN_FILE;
  } else {
    *(ftype) = GRID_FILE;
  }
  fseek(in,(long)(-sizeof(int)),SEEK_CUR);

  if (*(ftype) == GRID_FILE) {

    /* figure out if we're single/double precision and whole or planes */
    npts = 1; for (i = 0; i < *(ND); i++) npts *= N[i]; 
    safe_fread((long int)&record,sizeof(int),1,in);
    if (record == (*(ND) * npts * sizeof(float) + npts * sizeof(int))) {
      sprintf(vf,"%s","w");
      sprintf(prec,"%s","s");
      sprintf(ib,"%s","y");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(double) + npts * sizeof(int))) {
      sprintf(vf,"%s","w");
      sprintf(prec,"%s","d");
      sprintf(ib,"%s","y");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(float))) {
      sprintf(vf,"%s","w");
      sprintf(prec,"%s","s");
      sprintf(ib,"%s","n");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(double))) {
      sprintf(vf,"%s","w");
      sprintf(prec,"%s","d");
      sprintf(ib,"%s","n");
      fclose(in);
      return;
    }

    /* we must be planes */
    npts = 1; for (i = 0; i < *(ND)-1; i++) npts *= N[i]; 
    if (record == (*(ND) * npts * sizeof(float)) + npts * sizeof(int)) {
      sprintf(vf,"%s","p");
      sprintf(prec,"%s","s");
      sprintf(ib,"%s","y");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(double) + npts * sizeof(int))) {
      sprintf(vf,"%s","p");
      sprintf(prec,"%s","d");
      sprintf(ib,"%s","y");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(float))) {
      sprintf(vf,"%s","p");
      sprintf(prec,"%s","s");
      sprintf(ib,"%s","n");
      fclose(in);
      return;
    } else if (record == (*(ND) * npts * sizeof(double))) {
      sprintf(vf,"%s","p");
      sprintf(prec,"%s","d");
      sprintf(ib,"%s","n");
      fclose(in);
      return;
    }


  }

  if (*(ftype) == SOLN_FILE) {

    /* figure out if we're single/double precision */
    sprintf(prec,"%s","u");
    safe_fread((long int)&record,sizeof(int),1,in);
    if (record == ( 4 * sizeof(float) )) {
      sprintf(prec,"%s","s");
      el_size = 4;
    } else if (record == ( 4 * sizeof(double) )) {
      sprintf(prec,"%s","d");
      el_size = 8;
    } else sprintf(prec,"%s","u");
    fseek(in,record,SEEK_CUR);
    safe_fread((long int)&record,sizeof(int),1,in);

    /* figure out if we're planes or whole */
    npts = 1; for (i = 0; i < *(ND); i++) npts *= N[i]; 
    safe_fread((long int)&record,sizeof(int),1,in);
    if (record == ((*(ND)+2) * npts * el_size)) {
      sprintf(vf,"%s","w");
      fclose(in);
      return;
    }

    /* we must be planes */
    npts = 1; for (i = 0; i < *(ND)-1; i++) npts *= N[i]; 
    if (record == ((*(ND)+2) * npts * el_size)) {
      sprintf(vf,"%s","p");
      fclose(in);
      return;
    }

    sprintf(vf,"%s","u");

  }

  fclose(in);
  return;

}

size_t safe_fread(long int addr, int size_el, int num_el, FILE *in)
{
  size_t retval;
  char *a;

  a = (char *)addr;
  
  retval = fread(a, size_el, num_el, in);
 
  #ifdef LE
  endianSwap((long int)a, num_el, size_el);
  #endif

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
