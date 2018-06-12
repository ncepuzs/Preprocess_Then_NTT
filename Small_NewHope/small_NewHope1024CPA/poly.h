#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include "params.h"

/* 
 * Elements of R_q = Z_q[X]/(X^n + 1). Represents polynomial
 * coeffs[0] + X*coeffs[1] + X^2*xoeffs[2] + ... + X^{n-1}*coeffs[n-1] 
 */
typedef struct {
  uint16_t coeffs[NEWHOPE_N];
} poly __attribute__ ((aligned (32)));

typedef struct {
  uint16_t coeffs[NEWHOPE_N/4];
} poly_quarter __attribute__ ((aligned (32)));

void poly_uniform(poly *a, const unsigned char *seed);
void poly_sample(poly *r, const unsigned char *seed, unsigned char nonce);
void poly_add(poly *r, const poly *a, const poly *b);

void poly_ntt(poly_quarter *r);
void poly_invntt(poly_quarter *r);

void poly_frombytes(poly *r, const unsigned char *a);
void poly_tobytes(unsigned char *r, const poly *p);
void poly_compress(unsigned char *r, const poly *p);
void poly_decompress(poly *r, const unsigned char *a);

void poly_frommsg(poly *r, const unsigned char *msg);
void poly_tomsg(unsigned char *msg, const poly *x);
void poly_sub(poly *r, const poly *a, const poly *b);


void poly_quarter_add(poly_quarter *r, const poly_quarter *a, const poly_quarter *b, const poly_quarter *c, const poly_quarter *d);
int split_poly(poly *f, poly_quarter *f00, poly_quarter *f01, poly_quarter *f10, poly_quarter *f11);
int shift_poly(poly_quarter *f, poly_quarter *g);
int recover_poly(poly *f, poly_quarter *f00, poly_quarter *f01, poly_quarter *f10, poly_quarter *f11);
void poly_quarter_add(poly_quarter *r, const poly_quarter *a, const poly_quarter *b, const poly_quarter *c, const poly_quarter *d);
void poly_quarter_mul_pointwise(poly_quarter *r, const poly_quarter *a, const poly_quarter *b);

#endif
