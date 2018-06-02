#ifndef POLY_H
#define POLY_H

#include <stdint.h>
// #include "params.h"
#define NEWHOPE_N 1024

/* 
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] 
 */
typedef struct {
  uint16_t coeffs[NEWHOPE_N];
} poly __attribute__ ((aligned (32)));

typedef struct {
  uint16_t coeffs[NEWHOPE_N/2];
} poly_half __attribute__ ((aligned (32)));

void poly_uniform(poly *a, const unsigned char *seed);
void poly_sample(poly *r, const unsigned char *seed, unsigned char nonce);
void poly_add(poly *r, const poly *a, const poly *b);
void poly_half_sub(poly_half *r, const poly_half *a, const poly_half *b);

void poly_ntt_512(poly_half *r);
void poly_invntt_512(poly_half *r);

void poly_mul_pointwise(poly *r, const poly *a, const poly *b);
void poly_half_mul_pointwise(poly_half *r, const poly_half *a, const poly_half *b);

void poly_ntt_1024(poly_half *r);
void poly_invntt_1024(poly_half *r);



#endif
