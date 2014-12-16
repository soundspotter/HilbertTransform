/* hilbert.h - hilbert transform via allpass filters
 * Michael A. Casey - Dartmouth College, USA
 */

#ifndef __HILBERT_H
#define __HILBERT_H
#include "iirfilter.h"

/* #define __HILBERTTEST__ */ /* Set to compile main() test program */

typedef struct SAMPLETVEC {
  sampleT* data;
  int len;
} vec;

typedef struct HILBERTSTRUCT{
  double fs;
  FILTER** ap_pair;
  int buflen;
  vec* imvec;
  sampleT* y1;
  sampleT* y2;
  sampleT* cpx;
  sampleT phase;
} Hilbert;


/* API */
Hilbert * init_hilbert(int buffer_length, double fs);
void free_hilbert(Hilbert* h);
void analytic(Hilbert* h, sampleT* x);
void freq_shift(Hilbert* h, sampleT* x, double f0);

#endif
