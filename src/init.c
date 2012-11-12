/*
 *  Copyright (C) 2009-2010  B. Wang
 *  Unlimited use and distribution (see LICENCE).
 */

#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

/*
void F77_SUB(FFTSupport)(double *DT, int *NDT, double *DLO, double *DHI,
			 double *WINDOW, double *SIG, double *FT,
			 double *SMOOTH, int *NFT);
void llrHLaplace(double *range, double *x, double *sig, double *y,
		 double *scb, double *bw, int *bwyn, int *size, 
		 double *sens, double *level);
void DMLaplace(double *xv, int *size, double *bw, double *sig, double *l);
void LLSLaplace(double *y, double *x, double *sig, int *size,
		double *bw, double *rhat);
void lpgauss(double *y, double *x, double *f, int *size,
	      double *bw, double *bws, int *k, double *sensitivity,
	      double *lx, double *rhat);

void NormNorm1(int *n, double *Rfx, double *s2, double *h1, 
	       double *grid, double *ub);

void bwBoot1(double *h0,int *size,int *type,double *y,double *sig, 
	     int *grid,double *ub);

void F77_SUB(npqlevel)(double *x, int *n, double *x0, double *y0, int *m);

void F77_SUB(gaussbin)(double *x, double *sigma, int *n, double *a, double *b, int *m, 
		     int *trun, double *gcounts);
					    
void F77_SUB(linbin)(double *x, int *n, double *a,
		     double *b, int *m, int *trun, double *gcounts);
*/


static const R_FortranMethodDef FortEntries[] = {
  /*
  {"llrHLaplace", (DL_FUNC) & llrHLaplace, 10},
  {"DMLaplace", (DL_FUNC) & DMLaplace, 5},
  {"LLSLaplace", (DL_FUNC) & LLSLaplace, 6},
  {"lpgauss", (DL_FUNC) & lpgauss, 10},
  {"FFTSupport", (DL_FUNC) &F77_SUB(FFTSupport), 9},
  {"gaussbin", (DL_FUNC) &F77_SUB(gaussbin),  8},
  {"linbin", (DL_FUNC) &F77_SUB(linbin),  7},
  {"npqlevel", (DL_FUNC) &F77_SUB(npqlevel),  5},
  {"NormNorm1", (DL_FUNC) & NormNorm1, 6},
  {"bwBoot1", (DL_FUNC) & bwBoot1, 7},
*/
  {NULL, NULL, 0}
};


void R_init_miRada(DllInfo *dll)
{
  //    R_registerRoutines(dll, NULL, NULL, callMethods, NULL);
  R_registerRoutines(dll, NULL, NULL, FortEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
