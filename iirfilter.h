/* iirfilter.h */

/* Author: Michael A. Casey
 * Language: C
 * Copyright (C) 1997 Michael A. Casey, MIT Media Lab, All Rights Reserved
 *
 * Implementation of filter opcode for general purpose filtering.
 * This opcode implements the following difference equation:
 *
 * (1)*y(n) = b(0)*x(n) + b(1)*x(n-1) + ... + b(nb)*x(n-nb)
 *                      - a(1)*y(n-1) - ... - a(na)*y(n-na)
 *
 * whose system function is represented by:
 *
 *                               -1              -nb
 *        jw  B(z)   b(0) + b(1)z + .... + b(nb)z
 *     H(e) = ---- = ----------------------------
 *                               -1              -na
 *            A(z)    1   + a(1)z + .... + a(na)z
 *
 *
 * This is the same as scipy.signal's lfilter is a Direct Form II Transposed Filter:
 *                             -1               -nb
 *                 b[0] + b[1]z  + ... + b[nb] z
 *         Y(z) = ---------------------------------- X(z)
 *                             -1               -na
 *                 a[0] + a[1]z  + ... + a[na] z    
 *
 */

#ifndef __filter_h
#define __filter_h
#include <stdint.h>

/*#define __FILTERTEST__*/
#define MAXZEROS 50 /* Allow up to 50th-order digital filters */
#define MAXPOLES 50
#define OK 104
#define CS_KSMPS 4096 /* samples per buffer */

typedef double sampleT;
typedef struct FCOMPLEX {sampleT r,i;} fcomplex;

/* Structures for FILTER opcode */
typedef struct {
  sampleT *out;       /* output signal */
  sampleT *in;        /* input signal */
  int *nb, *na;   /* filter-order input arguments */
  sampleT coeffs[MAXPOLES+MAXZEROS+1]; /* filter-coefficient input arguments */
  sampleT *d1,*d2;    /* These allow ZFILTER to access FILTER routines */

  int numa;         /* i-var p-time storage registers */
  int numb;

  sampleT* delay;     /* delay-line state memory base pointer */
  sampleT* currPos;  /* delay-line current position pointer */ /* >>Was float<< */
  int   ndelay;    /* length of delay line (i.e. filter order) */
} FILTER;

typedef struct {
  sampleT *out;       /* output signal */
  sampleT *in;        /* input signal */
  sampleT *kmagf, *kphsf; /* magnitude and phase pole nudging factors */
  int *nb, *na;   /* filter-order input arguments */
  sampleT coeffs[MAXPOLES+MAXZEROS+1]; /* filter-coefficient input arguments */

  int numa;         /* i-var p-time storage registers */
  int numb;

  sampleT* delay;       /* delay-line state memory base pointer */
  sampleT* currPos;  /* delay-line current position pointer */ /* >>Was float<< */
  int   ndelay;      /* length of delay line (i.e. filter order) */
  fcomplex* roots;       /* pole roots memory for zfilter */
} ZFILTER;

/* API */
int ifilter(FILTER* p);
void free_filter(FILTER* p);
int izfilter(ZFILTER *p);
void free_zfilter(ZFILTER* p);
int afilter(FILTER* p, uint32_t nsmps);
int azfilter(ZFILTER* p, uint32_t nsmps);
int kfilter(FILTER* p);

#endif



