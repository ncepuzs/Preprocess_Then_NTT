#ifndef POLYVEC_H
#define POLYVEC_H

#include "params.h"
#include "poly.h"
#include "pt_ntt.h"

typedef struct{
  poly vec[KYBER_K];
} polyvec;

typedef struct{
  poly_half vec[KYBER_K];
} poly_halfvec;


void polyvec_compress(unsigned char *r, const polyvec *a);
void polyvec_decompress(polyvec *r, const unsigned char *a);

void polyvec_tobytes(unsigned char *r, const polyvec *a);
void polyvec_frombytes(polyvec *r, const unsigned char *a);

// void polyvec_ntt(polyvec *r);
// void polyvec_invntt(polyvec *r);
  
// void polyvec_pointwise_acc(poly *r, const polyvec *a, const polyvec *b);

void polyvec_add(polyvec *r, const polyvec *a, const polyvec *b);


int polyvec_pt_ntt2(polyvec *a);
int polyvec_pt_ntt3(polyvec *s, poly_halfvec *s0, poly_halfvec *s1, poly_halfvec *s2);
int polyvec_pt_ntt_bowtiemultiply(poly *p, polyvec *a, poly_halfvec *s0, poly_halfvec *s1, poly_halfvec *s2);
int polyvec_inv_ptntt(polyvec *a);


#endif
