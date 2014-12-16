/* hilbert.c - hilbert transform via allpass filters
 *
 * Ideas from: http://yehar.com/blog/?p=368
 * originally from http://ldesoras.free.fr/prod.html#src_hiir
 * Also see http://www.mathworks.com/help/signal/examples/single-sideband-modulation-via-the-hilbert-transform.html
 *
 * Copyright (C) 2014 Michael A. Casey, Bregman Media Labs, Dartmouth College, USA
 */

/*****************************************************************************************************************************/
/* Here's a quick diagram of the allpass pair:										     */
/* 															     */
/*         ................ filter 1 .................									     */
/*    +--> allpass --> allpass --> allpass --> allpass --> delay --> out1						     */
/*    |															     */
/*   in															     */
/*    |    ................ filter 2 .................									     */
/*    +--> allpass --> allpass --> allpass --> allpass ------------> out2 (+90 deg)					     */
/* 															     */
/* We can use cookbook formulas to convert an allpass section into code. A general IIR recurrence relation:		     */
/* 															     */
/*   out(t) = a0*in(t) + a1*in(t-1) + a2*in(t-2) + ...									     */
/*          + b1*out(t-1) + b2*out(t-2) + ...										     */
/* 															     */
/* results in the transfer function:											     */
/* 															     */
/*          a0 + a1*z^-1 + a2*z^-2 + ...										     */
/*   H(z) = ----------------------------										     */
/*           1 - b1*z^-1 - b2*z^-2 - ...										     */
/* 															     */
/* The allpass section in question has the following transfer function:							     */
/* 															     */
/*           a^2 - z^-2													     */
/*   H(z) = ------------												     */
/*          1 - a^2 z^-2												     */
/* 															     */
/* We want to convert this into the recurrence relation. According to the cookbook formulas and the above transfer function: */
/* 															     */
/*   a0 = a^2, a2 = -1, b2 = a^2, rest of coefficients zero								     */
/* 															     */
/*   =>  out(t) = a^2*in(t) - in(t-2) + a^2*out(t-2)									     */
/* 															     */
/* which simplifies to the one-multiplication allpass section:								     */
/* 															     */
/*   out(t) = a^2*(in(t) + out(t-2)) - in(t-2)										     */
/*****************************************************************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "hilbert.h"

sampleT HilbertCoeffs[2][4] = {{0.6923878, 0.9360654322959, 0.9882295226860, 0.9987488452737},
				{0.4021921162426, 0.8561710882420, 0.9722909545651, 0.9952884791278}};

static vec ** H_sect(sampleT a);
static void free_H_sect(vec** allpass_coeffs);
static vec * imp(size_t N);
static vec* convolve(vec* X, vec* Y);
static vec ** H(vec* A0);
static vec ** H_1(vec* A0);
static vec ** H_2(vec* A0);
static FILTER* init_biquad_filter(vec* B, vec* A);
static void free_allpass_filters(FILTER** ap_pair);
static FILTER ** init_allpass_filters();


/* Utility functions */

static vec* new_vec(int N, int allocate){
  vec* v = (vec*) calloc(1,sizeof(vec));
  if(allocate){
    v->data = (sampleT*) calloc(N,sizeof(sampleT));
  }
  else{
    v->data = NULL; 
  }
  v->len = N;
  return v;
}

static void free_vec(vec* v){
  if(v->data!=NULL){
    free(v->data);
  }
  free(v);
}


/*
 * Impulse signal
 */
static vec * imp(size_t N){
  vec *I = new_vec(N,1);
  I->data[0] = 1.0;
  return I;
} 

/*
 * Full direct convolution of two vectors
 * Returns a new vector of length N1+N2-1
 */
static vec* convolve(vec* X, vec* Y){
  int i,k;
  vec* v = new_vec(X->len+Y->len-1, 1);
  for(i=0; i<X->len; i++){
    for(k=0; k<Y->len; k++){
      v->data[i+k] += X->data[i] * Y->data[k];
    }
  }
  return v;
}



/*************************************************************/
/* Allpass biquadratic section:			   	     */
/* 							     */
/*                  a^2 - z^-2				     */
/*          H(z) = -------------			     */
/* 		    1 - a^2 z^-2			     */
/* 							     */
/* Implemenation as direct Form II Transposed Filter:	     */
/* 							     */
/*              -1                       -nb		     */
/*          b[0] + b[1]z  + ... + b[nb] z		     */
/*  Y(z) = ---------------------------------- X(z)	     */
/*                      -1               -na		     */
/*          a[0] + a[1]z  + ... + a[na] z		     */
/*************************************************************/
static vec ** H_sect(sampleT a){
  vec ** coeffs = (vec**) calloc(2,sizeof(vec*));
  coeffs[0] = new_vec(4,1);
  coeffs[1] = new_vec(4,1);
  coeffs[0]->data[0] = a*a;
  coeffs[0]->data[2] = -1.0;
  coeffs[1]->data[0] = 1.0;
  coeffs[1]->data[2] = -(a*a);
  return coeffs;
}

static void free_H_sect(vec** allpass_coeffs){
 free_vec(allpass_coeffs[0]);free_vec(allpass_coeffs[1]);
 free(allpass_coeffs);
}

/*
 * Allpass hilbert transform sections denominator
 */
static vec ** H(vec* A0){
  vec *B1, *A1, *B2, *A2, *B3, *A3, *B4, *A4, **tmp;
  vec *D2, *C2, *D3, *C3, *outB, *outA;
  sampleT* a0 = A0->data;
  tmp = H_sect(a0[0]); B1=tmp[0]; A1=tmp[1];
  free(tmp); // avoid memory leaks
  tmp = H_sect(a0[1]); B2=tmp[0]; A2=tmp[1];
  free(tmp); // avoid memory leaks
  tmp = H_sect(a0[2]); B3=tmp[0]; A3=tmp[1];
  free(tmp); // avoid memory leaks
  tmp = H_sect(a0[3]); B4=tmp[0]; A4=tmp[1];
  D2 = convolve(B1,B2); C2 = convolve(A1,A2);
  D3 = convolve(D2,B3); C3 = convolve(C2,A3);
  outB = convolve(D3,B4); outA = convolve(C3,A4);
  free_vec(B1);free_vec(A1);
  free_vec(B2);free_vec(A2);
  free_vec(B3);free_vec(A3);
  free_vec(B4);free_vec(A4);
  free_vec(D2);free_vec(C2);
  free_vec(D3);free_vec(C3);
  tmp[0]=outB; tmp[1]=outA;
  return tmp;
}

/*
 * Allpass hilbert transform numerator
 */
static vec ** H_1(vec* A0){
  vec ** tmp = H(A0);
  vec* unit_delay = new_vec(2,1);
  unit_delay->data[1]=1;
  vec *B  = convolve(tmp[0], unit_delay);
  free_vec(unit_delay);
  free_vec(tmp[0]); // avoid memory leaks
  tmp[0] = B;
  return tmp;
}

/*
 * Allpass hilbert transform denominator
 */
static vec ** H_2(vec* A0){
  return H(A0);
}

/*
 * Allpass filters from series biquad filter sections
 */
static FILTER* init_biquad_filter(vec* B, vec* A){
  FILTER* biquad = (FILTER*) calloc(1,sizeof(FILTER));
  biquad->numb = B->len;
  biquad->numa = A->len-1; // Assume A[0]=1 and crop array
  int i;
  for(i=0; i<biquad->numb; i++){
	biquad->coeffs[i] = B->data[i];
  }
  for(i=1; i<biquad->numa; i++){ // Assume A[0]=1 and crop array
    biquad->coeffs[biquad->numb+i-1] = A->data[i];
  }
  ifilter(biquad);
  return biquad;
}

/* 
 * Deallocate allpass filter (multi-section numerator and denominator)
 */
static void free_allpass_filters(FILTER** ap_pair){
  free_filter(ap_pair[0]);
  free_filter(ap_pair[1]);
  free(ap_pair);
}

/*
 * Initialize allpass filter (multi-section numerator and denominator)
 */
static FILTER ** init_allpass_filters(){
  vec ** allpass_coeffs;
  vec* a1 = new_vec(4, 0);
  vec* a2 = new_vec(4, 0);  
  a1->data = HilbertCoeffs[0];
  a2->data = HilbertCoeffs[1];  
  allpass_coeffs = H_1(a1);
  FILTER* allpass1 = init_biquad_filter(allpass_coeffs[0], allpass_coeffs[1]);
  free_H_sect(allpass_coeffs);

  allpass_coeffs = H_2(a2);  
  FILTER* allpass2 = init_biquad_filter(allpass_coeffs[0], allpass_coeffs[1]);
  free_H_sect(allpass_coeffs);

  FILTER** ap_pair = (FILTER**) calloc(2, sizeof(FILTER*));
  ap_pair[0] = allpass1;
  ap_pair[1] = allpass2;
  a1->data = NULL; a2->data = NULL;
  free_vec(a1); free_vec(a2);
  return ap_pair;  
}

/* Hilbert transform constructor */
Hilbert* init_hilbert(int buffer_length, double fs){
  Hilbert* H = (Hilbert*) calloc(1, sizeof(Hilbert));
  H->fs = fs;
  H->buflen = buffer_length;
  H->ap_pair = init_allpass_filters();
  H->imvec = new_vec(H->buflen,1); // Imaginary part of analytic signal
  H->y1 = (sampleT*) calloc(H->buflen,sizeof(sampleT));
  H->y2 = (sampleT*) calloc(H->buflen,sizeof(sampleT));
  H->cpx = (sampleT*) calloc(H->buflen,sizeof(sampleT));
  H->ap_pair[0]->out = H->y1;
  H->ap_pair[1]->out = H->y2;
  H->phase = 0.0;
  if(!(H->ap_pair&&H->buflen&&H->imvec&&H->y1&&H->y2)){
    fprintf(stderr, "Hilbert transformer initialization failed init_hilbert()\n");
    exit(1);
  }
  return H;
}

/* Hilbert transform destructor */
void free_hilbert(Hilbert* H){
  free(H->y1);
  free(H->y2);
  free(H->cpx);
  free_vec(H->imvec);
  free_allpass_filters(H->ap_pair);
  free(H);
}

/*
  Analytic (complex) signal, in-place, given real and zero-imag input.
    H := z -> 0.5*(H_2(z)+I*H_1(z));
*/
void analytic(Hilbert* H, sampleT* x){
  int i;
  if(!H){
    fprintf(stderr, "Hilbert transformer initialization not previously called, analytic()\n");    
    exit(1);
  }
  H->ap_pair[0]->in = x;
  afilter(H->ap_pair[0], H->buflen);
  H->ap_pair[1]->in = x;
  afilter(H->ap_pair[1], H->buflen);  
  for(i=0; i<H->buflen; i++){
    x[i] = H->ap_pair[1]->out[i]*0.5;
    H->cpx[i] = H->ap_pair[0]->out[i]*0.5;
  }
}

/*
    In-place hilbert transformer frequency shifter, by constant offset
    Uses single sideband modulation of input signal to carrier (offset) 
*/
void freq_shift(Hilbert* H, sampleT* x, double f0){
  double ws = 2*M_PI*f0/H->fs; // Carrier freq
  int i;
  analytic(H, x);
  for(i = 0; i<H->buflen; i++){
    x[i] = 2 * (cos(ws*i+H->phase)*x[i] - sin(ws*i+H->phase)*H->cpx[i]);
  }  
  H->phase = fmod(H->phase + ws*H->buflen, 2*M_PI);
}

#ifdef __HILBERTTEST__
int main(int argc, char* argv[]){
  sampleT* in = (sampleT*)calloc(CS_KSMPS, sizeof(sampleT));
  if(!in){
    fprintf(stderr, "Could not allocate in buffer.\n");
    exit(1);
  }
  Hilbert* H = init_hilbert(CS_KSMPS, 44100.0); // initialize hilbert transform freq shifter
  double w0 = 2 * M_PI * 110.0 / H->fs;
  double phase = 0.0;
  // 4 buffers of samples
  int i,j,k;
  for (i=0; i < 4 ; i++){
    for(j=0; j < H->buflen; j++){
      in[j]=cos(w0*j+phase);
    }
    phase = fmod(phase + w0*H->buflen, 2*M_PI);
    freq_shift(H, in, 10.0);
    for(j=0; j < H->buflen; j++){
      fprintf(stdout, "%5.4f ", in[j]);
    }
  }
  fprintf(stdout, "\n");
  free_hilbert(H);
  free(in);
  exit(0);
}

#endif
