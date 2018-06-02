#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "fips202.h"

#define NEWHOPE_N 1024
#define NEWHOPE_Q 12289
#define NEWHOPE_SYMBYTES 32

/*************************************************
* Name:        coeff_freeze
* 
* Description: Fully reduces an integer modulo q in constant time
*
* Arguments:   uint16_t x: input integer to be reduced
*              
* Returns integer in {0,...,q-1} congruent to x modulo q
**************************************************/
static uint16_t coeff_freeze(uint16_t x)
{
  uint16_t m,r;
  int16_t c;
  r = x % NEWHOPE_Q;

  m = r - NEWHOPE_Q;
  c = m;
  c >>= 15;
  r = m ^ ((r^m)&c);

  return r;
}

/*************************************************
* Name:        flipabs
* 
* Description: Computes |(x mod q) - Q/2|
*
* Arguments:   uint16_t x: input coefficient
*              
* Returns |(x mod q) - Q/2|
**************************************************/
static uint16_t flipabs(uint16_t x)
{
  int16_t r,m;
  r = coeff_freeze(x);

  r = r - NEWHOPE_Q/2;
  m = r >> 15;
  return (r + m) ^ m;
}
 
/*************************************************
* Name:        poly_uniform
* 
* Description: Sample a polynomial deterministically from a seed,
*              with output polynomial looking uniformly random
*
* Arguments:   - poly *a:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed
**************************************************/
void poly_uniform(poly *a, const unsigned char *seed)
{
  unsigned int ctr=0;
  uint16_t val;
  uint64_t state[25];
  uint8_t buf[SHAKE128_RATE];
  uint8_t extseed[NEWHOPE_SYMBYTES+1];
  int i,j;

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];

  for(i=0;i<NEWHOPE_N/64;i++) /* generate a in blocks of 64 coefficients */
  {
    ctr = 0;
    extseed[NEWHOPE_SYMBYTES] = i; /* domain-separate the 16 independent calls */
    shake128_absorb(state, extseed, NEWHOPE_SYMBYTES+1);
    while(ctr < 64) /* Very unlikely to run more than once */
    {
      shake128_squeezeblocks(buf,1,state);
      for(j=0;j<SHAKE128_RATE && ctr < 64;j+=2)
      {
        val = (buf[j] | ((uint16_t) buf[j+1] << 8));
        if(val < 5*NEWHOPE_Q)
        {
          a->coeffs[i*64+ctr] = val;
          ctr++;
        }
      }
    }
  }
}

/*************************************************
* Name:        hw
* 
* Description: Compute the Hamming weight of a byte
*
* Arguments:   - unsigned char a: input byte
**************************************************/
static unsigned char hw(unsigned char a)
{
  unsigned char i, r = 0;
  for(i=0;i<8;i++)
    r += (a >> i) & 1;
  return r;
}

/*************************************************
* Name:        poly_sample
* 
* Description: Sample a polynomial deterministically from a seed and a nonce,
*              with output polynomial close to centered binomial distribution
*              with parameter k=8
*
* Arguments:   - poly *r:                   pointer to output polynomial
*              - const unsigned char *seed: pointer to input seed 
*              - unsigned char nonce:       one-byte input nonce
**************************************************/
void poly_sample(poly *r, const unsigned char *seed, unsigned char nonce)
{
  unsigned char buf[128], a, b;
//  uint32_t t, d, a, b, c;
  int i,j;

  unsigned char extseed[NEWHOPE_SYMBYTES+2];

  for(i=0;i<NEWHOPE_SYMBYTES;i++)
    extseed[i] = seed[i];
  extseed[NEWHOPE_SYMBYTES] = nonce;

  for(i=0;i<NEWHOPE_N/64;i++) /* Generate noise in blocks of 64 coefficients */
  {
    extseed[NEWHOPE_SYMBYTES+1] = i;
    shake256(buf,128,extseed,NEWHOPE_SYMBYTES+2);
    for(j=0;j<64;j++)
    {
      a = buf[2*j];
      b = buf[2*j+1];
      r->coeffs[64*i+j] = hw(a) + NEWHOPE_Q - hw(b);
      /*
      t = buf[j] | ((uint32_t)buf[j+1] << 8) | ((uint32_t)buf[j+2] << 16) | ((uint32_t)buf[j+3] << 24);
      d = 0;
      for(k=0;k<8;k++)
        d += (t >> k) & 0x01010101;
      a = d & 0xff;
      b = ((d >>  8) & 0xff);
      c = ((d >> 16) & 0xff);
      d >>= 24;
      r->coeffs[64*i+j/2]   = a + NEWHOPE_Q - b;
      r->coeffs[64*i+j/2+1] = c + NEWHOPE_Q - d;
      */
    }
  }
}

/*************************************************
* Name:        poly_pointwise
* 
* Description: Multiply two polynomials pointwise (i.e., coefficient-wise).
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_mul_pointwise(poly *r, const poly *a, const poly *b)
{
  int i;
  uint16_t t;
  for(i=0;i<NEWHOPE_N;i++)
  {
    t            = montgomery_reduce(3186*b->coeffs[i]); /* t is now in Montgomery domain */
    r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t);  /* r->coeffs[i] is back in normal domain */
  }
}
void poly_half_mul_pointwise(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i;
  uint16_t t;
  for(i=0;i<NEWHOPE_N/2;i++)
  {
    t            = montgomery_reduce(3186*b->coeffs[i]); /* t is now in Montgomery domain */
    r->coeffs[i] = montgomery_reduce(a->coeffs[i] * t);  /* r->coeffs[i] is back in normal domain */
  }
}

/*************************************************
* Name:        poly_add
* 
* Description: Add two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_add(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<NEWHOPE_N;i++)
    r->coeffs[i] = (a->coeffs[i] + b->coeffs[i]) % NEWHOPE_Q;
}
void poly_half_add(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i;
  for(i=0;i<NEWHOPE_N/2;i++)
    r->coeffs[i] = (a->coeffs[i] + b->coeffs[i]) % NEWHOPE_Q;
}
/*************************************************
* Name:        poly_sub
* 
* Description: Subtract two polynomials
*
* Arguments:   - poly *r:       pointer to output polynomial
*              - const poly *a: pointer to first input polynomial
*              - const poly *b: pointer to second input polynomial
**************************************************/
void poly_sub(poly *r, const poly *a, const poly *b)
{
  int i;
  for(i=0;i<NEWHOPE_N;i++)
    r->coeffs[i] = (a->coeffs[i] + 3*NEWHOPE_Q - b->coeffs[i]) % NEWHOPE_Q;
}
void poly_half_sub(poly_half *r, const poly_half *a, const poly_half *b)
{
  int i;
  for(i=0;i<NEWHOPE_N/2;i++)
    r->coeffs[i] = (a->coeffs[i] + 3*NEWHOPE_Q - b->coeffs[i]) % NEWHOPE_Q;
}

/*************************************************
* Name:        poly_ntt
* 
* Description: Forward NTT transform of a polynomial in place
*              Input is assumed to have coefficients in bitreversed order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_ntt_512(poly_half *r)
{
  mul_coefficients_512(r->coeffs, psis_bitrev_montgomery_512);
  ntt_512((uint16_t *)r->coeffs, omegas_bitrev_montgomery_512);
}
void poly_ntt_1024(poly_half *r)
{
  mul_coefficients_1024(r->coeffs, psis_bitrev_montgomery_1024);
  ntt_1024((uint16_t *)r->coeffs, omegas_bitrev_montgomery_1024);
}

/*************************************************
* Name:        poly_invntt
* 
* Description: Inverse NTT transform of a polynomial in place
*              Input is assumed to have coefficients in normal order
*              Output has coefficients in normal order
*
* Arguments:   - poly *r: pointer to in/output polynomial
**************************************************/
void poly_invntt_512(poly_half *r)
{
  bitrev_vector_512(r->coeffs);
  ntt_512((uint16_t *)r->coeffs, omegas_inv_bitrev_montgomery_512);
  mul_coefficients_512(r->coeffs, psis_inv_montgomery_512);
  for(int i=0;i<512;i++)
  {
    if(r->coeffs[i]>=NEWHOPE_Q)
      r->coeffs[i]=r->coeffs[i]%NEWHOPE_Q;
  }
}
void poly_invntt_1024(poly_half *r)
{
  bitrev_vector_1024(r->coeffs);
  ntt_1024((uint16_t *)r->coeffs, omegas_inv_bitrev_montgomery_1024);
  mul_coefficients_1024(r->coeffs, psis_inv_montgomery_1024);
}
/*************************************************
* Name:        poly_invntt
* 
* Description: test if two polynomials equal
*              Input are two polynomials with degree = NEWHOPE_N
*              Output is 0 or 1
*
* Arguments:   - poly a: in/output polynomial
**************************************************/
int poly_equal(poly a,poly b)
{
  int i;
  for(i=0;i<NEWHOPE_N;i++)
    if(a.coeffs[i]!=b.coeffs[i])
      return 0;
  return 1;
}
/*************************************************
* Name:        poly_invntt
* 
* Description: print the value of a polynomial
*              Input is a polynomials with degree = NEWHOPE_N
*              Output the value 
*
* Arguments:   - poly a: the pointer of in/output polynomial
**************************************************/
int print_poly(poly *p)
{
    int i;
    for(i=0;i<NEWHOPE_N;i++)
    {
        printf("%d,",p->coeffs[i]);
    }
    printf("\n");
    return 1;
}

int split_poly(poly *f, poly_half *f0, poly_half *f1)
{
    int i;
    for(i=0;i<NEWHOPE_N/2;i++)
    {
        f0->coeffs[i]=f->coeffs[2*i];
        f1->coeffs[i]=f->coeffs[2*i+1];
    }
    return 1;
}
int shift_poly(poly_half *f, poly_half *g)
{
    int i;
    //coefficients in polynomials of Kyber < q = 3329 which need 2 bytes
    for(i=1;i<NEWHOPE_N/2;i++)
    {
        f->coeffs[NEWHOPE_N/2-i]=g->coeffs[NEWHOPE_N/2-i-1];
    }
    f->coeffs[0]=NEWHOPE_Q-(g->coeffs[NEWHOPE_N/2-1]%NEWHOPE_Q);

    return 1;
}
int recover_poly(poly *f, poly_half *f0, poly_half *f1)
{
    int i;
    for(i=0;i<NEWHOPE_N/2;i++)
    {
        f->coeffs[2*i] =   f0->coeffs[i];
        f->coeffs[2*i+1] = f1->coeffs[i];
    }
    return 1;
}
