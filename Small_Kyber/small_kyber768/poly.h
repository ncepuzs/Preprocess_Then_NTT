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


void poly_compress(unsigned char *r, const poly *a);
void poly_decompress(poly *r, const unsigned char *a);

void poly_tobytes(unsigned char *r, const poly *a);
void poly_frombytes(poly *r, const unsigned char *a);

void poly_frommsg(poly *r, const unsigned char msg[KYBER_SYMBYTES]);
void poly_tomsg(unsigned char msg[KYBER_SYMBYTES], const poly *r);

void poly_getnoise(poly *r,const unsigned char *seed, unsigned char nonce);

void poly_ntt(poly_half *r);
void poly_invntt(poly_half *r);
  
void poly_add(poly *r, const poly *a, const poly *b);
void poly_sub(poly *r, const poly *a, const poly *b);

int print_poly(poly *p);
int poly_equal(poly a,poly b);

void poly_half_add(poly_half *r, const poly_half *a, const poly_half *b);
void poly_half_mul_pointwise(poly_half *r, const poly_half *a, const poly_half *b);
int split_poly(poly *f, poly_half *f0, poly_half *f1);
int shift_poly(poly_half *f, poly_half *g);
int recover_poly(poly *f, poly_half *f0, poly_half *f1);


#endif
