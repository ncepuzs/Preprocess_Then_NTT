#include <stdio.h>
#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "fips202.h"

/*************************************************
* Name:        poly_getnoise
* 
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter KYBER_ETA
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed 
*              - unsigned char nonce:       one-byte input nonce
**************************************************/
void poly_getnoise(poly *r,const unsigned char *seed, unsigned char nonce)
{
  unsigned char buf[KYBER_ETA*KYBER_N/4];
  unsigned char extseed[KYBER_SYMBYTES+1];
  int i;

  for(i=0;i<KYBER_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[KYBER_SYMBYTES] = nonce;
     
  shake256(buf,KYBER_ETA*KYBER_N/4,extseed,KYBER_SYMBYTES+1);

  cbd(r, buf);
}

/*************************************************
* Name:        poly_ntt
* 
* Description: Computes negacyclic number-theoretic transform (NTT) of
*              a polynomial in place; 
*              inputs assumed to be in normal order, output in bitreversed order
*
* Arguments:   - uint16_t *r: pointer to in/output polynomial
**************************************************/
void poly_ntt(poly *r)
{
  ntt(r->coeffs);
}
void poly_ntt_128(poly_half *r)
{
  ntt_128(r->coeffs);
}

/*************************************************
* Name:        poly_invntt
* 
* Description: Computes inverse of negacyclic number-theoretic transform (NTT) of
*              a polynomial in place; 
*              inputs assumed to be in bitreversed order, output in normal order
*
* Arguments:   - uint16_t *a: pointer to in/output polynomial
**************************************************/
void poly_invntt(poly *r)
{
  invntt(r->coeffs);
}
void poly_invntt_128(poly_half *r)
{
  invntt_128(r->coeffs);
}
 
/*************************************************
* Name:        poly_add
* 
* Description: Add two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/ 
void poly_add(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = barrett_reduce(a->coeffs[i] + b->coeffs[i]);
}
void poly_half_add(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i;
  for(i=0;i<KYBER_N/2;i++)
    r->coeffs[i] = barrett_reduce(a->coeffs[i] + b->coeffs[i]);
}

/*************************************************
* Name:        poly_sub
* 
* Description: Subtract two polynomials
*
* Arguments: - poly *r:       pointer to output polynomial
*            - const poly *a: pointer to first input polynomial
*            - const poly *b: pointer to second input polynomial
**************************************************/ 
void poly_sub(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<KYBER_N;i++)
    r->coeffs[i] = barrett_reduce(a->coeffs[i] + 3*KYBER_Q - b->coeffs[i]);
}
void poly_half_mul_pointwise(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i,j;
  uint16_t t;
  for(j=0;j<KYBER_N/2;j++)
  {
    t = montgomery_reduce(4613* (uint32_t)b->coeffs[j]); // 1674 = 2^{2*18} % q
    r->coeffs[j] = montgomery_reduce(a->coeffs[j] * t);
  }
}


int poly_equal(poly a,poly b)
{
  int i;
  for(i=0;i<KYBER_N;i++)
    if(a.coeffs[i]!=b.coeffs[i])
      return 0;
  return 1;
}
int print_poly(poly *p)
{
    int i;
    for(i=0;i<KYBER_N;i++)
    {
        printf("%d,",p->coeffs[i]);
    }
    printf("\n");
    return 1;
}

void poly_pointwise_mul(poly *r, const poly *a, const poly *b)
{
  int i,j;
  uint16_t t;
  for(j=0;j<KYBER_N;j++)
  {
    t = montgomery_reduce(4613* (uint32_t)b->coeffs[j]); // 1674 = 2^{2*18} % q
    r->coeffs[j] = montgomery_reduce(a->coeffs[j] * t);
  }
}

int split_poly(poly *f, poly_half *f0, poly_half *f1)
{
    int i;
    for(i=0;i<KYBER_N/2;i++)
    {
        f0->coeffs[i]=f->coeffs[2*i];
        f1->coeffs[i]=f->coeffs[2*i+1];
    }
    return 1;
}
int shift_poly(poly_half *f, poly_half *g)
{
    int i;
    f->coeffs[0]=KYBER_Q-(g->coeffs[KYBER_N/2-1]%KYBER_Q);
    for(i=1;i<KYBER_N/2;i++)
    {
        f->coeffs[KYBER_N/2-i]=g->coeffs[KYBER_N/2-i-1];
    }
    return 1;
}
int recover_poly(poly *f, poly_half *f0, poly_half *f1)
{
    int i;
    for(i=0;i<KYBER_N/2;i++)
    {
        f->coeffs[2*i] =   f0->coeffs[i];
        f->coeffs[2*i+1] = f1->coeffs[i];
    }
    return 1;
}