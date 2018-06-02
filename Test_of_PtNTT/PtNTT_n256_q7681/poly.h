#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/* 
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] 
 */
typedef struct{
  uint16_t coeffs[KYBER_N];
} poly;

typedef struct{
  uint16_t coeffs[KYBER_N/2];
} poly_half;

void poly_getnoise(poly *r,const unsigned char *seed, unsigned char nonce);

void poly_ntt(poly *r);
void poly_invntt(poly *r);
void poly_ntt_128(poly_half *r);
void poly_invntt_128(poly_half *r);
  
void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);

#endif
