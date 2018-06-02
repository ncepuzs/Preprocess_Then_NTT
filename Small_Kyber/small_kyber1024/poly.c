#include <stdio.h>
#include "poly.h"
#include "ntt.h"
#include "polyvec.h"
#include "reduce.h"
#include "cbd.h"
#include "fips202.h"

/*************************************************
* Name:        poly_compress
* 
* Description: Compression and subsequent serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *a:    pointer to input polynomial
**************************************************/
void poly_compress(unsigned char *r, const poly *a)
{
  uint32_t t[8];
  unsigned int i,j,k=0;

  for(i=0;i<KYBER_N;i+=8)
  {
    for(j=0;j<8;j++)
      t[j] = (((freeze(a->coeffs[i+j]) << 3) + KYBER_Q/2)/KYBER_Q) & 7;

    r[k]   =  t[0]       | (t[1] << 3) | (t[2] << 6);
    r[k+1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
    r[k+2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
    k += 3;
  }
}

/*************************************************
* Name:        poly_decompress
* 
* Description: De-serialization and subsequent decompression of a polynomial; 
*              approximate inverse of poly_compress
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_decompress(poly *r, const unsigned char *a)
{
  unsigned int i;
  for(i=0;i<KYBER_N;i+=8)
  {
    r->coeffs[i+0] =  (((a[0] & 7) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+1] = ((((a[0] >> 3) & 7) * KYBER_Q)+ 4) >> 3;
    r->coeffs[i+2] = ((((a[0] >> 6) | ((a[1] << 2) & 4)) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+3] = ((((a[1] >> 1) & 7) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+4] = ((((a[1] >> 4) & 7) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+5] = ((((a[1] >> 7) | ((a[2] << 1) & 6)) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+6] = ((((a[2] >> 2) & 7) * KYBER_Q) + 4)>> 3;
    r->coeffs[i+7] = ((((a[2] >> 5)) * KYBER_Q) + 4)>> 3;
    a += 3;
  }
}

/*************************************************
* Name:        poly_tobytes
* 
* Description: Serialization of a polynomial
*
* Arguments:   - unsigned char *r: pointer to output byte array
*              - const poly *a:    pointer to input polynomial
**************************************************/
void poly_tobytes(unsigned char *r, const poly *a)
{
  int i,j;
  uint16_t t[2];

  for(i=0;i<KYBER_N/2;i++)
  {
    for(j=0;j<2;j++)
      t[j] = freeze(a->coeffs[2*i+j]);

    r[3*i+ 0] =  t[0]        & 0xff;
    r[3*i+ 1] = (t[0] >>  8) | ((t[1] & 0x0f) << 4);
    r[3*i+ 2] = (t[1] >>  4) & 0xff;
   
  }
}

/*************************************************
* Name:        poly_frombytes
* 
* Description: De-serialization of a polynomial; 
*              inverse of poly_tobytes
*
* Arguments:   - poly *r:                pointer to output polynomial
*              - const unsigned char *a: pointer to input byte array
**************************************************/
void poly_frombytes(poly *r, const unsigned char *a)
{
  int i;
  for(i=0;i<KYBER_N/2;i++)
  {
    r->coeffs[2*i+0] =  a[3*i+ 0]       | (((uint16_t)a[3*i+ 1] & 0x0f) << 8);
    r->coeffs[2*i+1] = (a[3*i+ 1] >> 4) | (((uint16_t)a[3*i+ 2]       ) << 4);
  }
}

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
void poly_ntt(poly_half *r)
{
  ntt(r->coeffs);
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
void poly_invntt(poly_half *r)
{
  invntt(r->coeffs);
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

/*************************************************
* Name:        poly_frommsg
* 
* Description: Convert 32-byte message to polynomial
*
* Arguments:   - poly *r:                  pointer to output polynomial
*              - const unsigned char *msg: pointer to input message
**************************************************/
void poly_frommsg(poly *r, const unsigned char msg[KYBER_SYMBYTES])
{
  uint16_t i,j,mask;

  for(i=0;i<KYBER_SYMBYTES;i++)
  {
    for(j=0;j<8;j++)
    {   
      mask = -((msg[i] >> j)&1);
      r->coeffs[8*i+j] = mask & ((KYBER_Q+1)/2);
    }   
  }
}

/*************************************************
* Name:        poly_tomsg
* 
* Description: Convert polynomial to 32-byte message
*
* Arguments:   - unsigned char *msg: pointer to output message
*              - const poly *a:      pointer to input polynomial
**************************************************/
void poly_tomsg(unsigned char msg[KYBER_SYMBYTES], const poly *a)
{
  uint16_t t;
  int i,j;

  for(i=0;i<KYBER_SYMBYTES;i++)
  {
    msg[i] = 0;
    for(j=0;j<8;j++)
    {
      t = (((freeze(a->coeffs[8*i+j]) << 1) + KYBER_Q/2)/KYBER_Q) & 1;
      msg[i] |= t << j;
    }
  }
}





void poly_half_mul_pointwise(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i,j;
  uint16_t t;
  for(j=0;j<KYBER_N/2;j++)
  {
    t = montgomery_reduce(1674* (uint32_t)b->coeffs[j]); // 1674 = 2^{2*18} % q
    r->coeffs[j] = montgomery_reduce(a->coeffs[j] * t);
  }
}


void poly_half_add(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i;
  for(i=0;i<KYBER_N/2;i++)
    r->coeffs[i] = barrett_reduce(a->coeffs[i] + b->coeffs[i]);
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