#ifndef output_h
#define output_h

#include <iostream>
#include <fstream>

#include "Environment.h"

struct outputFiles{
  FILE *analysis1;
  FILE *analysis2;
  FILE *analysis3;
  FILE *analysis4;
  FILE *analysis5;
  FILE *analysis6;
  FILE *analysis7;
  FILE *analysis8;
  FILE *fp;
  //FILE *fp,*fp2; // Output files for Rstudio analysis
  //FILE *fp4; // Input file for vessel position
  //  FILE *f1,*f2,*f3,*f4,*f5, *fmr, *fma, *fhrma, *fhama, *frneu, *faneu, *famc, *frmc; // Cell type files
//  FILE *fimmunekill, *ftotalcells;
};

#endif /* output_h */
