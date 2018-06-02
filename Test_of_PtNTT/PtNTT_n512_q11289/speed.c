#include "poly.h"
#include "pt_ntt.h"
#include "cpucycles.h"
#include <stdlib.h>
#include <stdio.h>

#define NTESTS 10000
#define NEWHOPE_SYMBYTES 32
static void gen_a(poly *a, const unsigned char *seed)
{
  poly_uniform(a,seed);
}

static int cmp_llu(const void *a, const void*b)
{
  if(*(unsigned long long *)a < *(unsigned long long *)b) return -1;
  if(*(unsigned long long *)a > *(unsigned long long *)b) return 1;
  return 0;
}

static unsigned long long median(unsigned long long *l, size_t llen)
{
  qsort(l,llen,sizeof(unsigned long long),cmp_llu);

  if(llen%2) return l[llen/2];
  else return (l[llen/2-1]+l[llen/2])/2;
}

static unsigned long long average(unsigned long long *t, size_t tlen)
{
  unsigned long long acc=0;
  size_t i;
  for(i=0;i<tlen;i++)
    acc += t[i];
  return acc/(tlen);
}

static void print_results(const char *s, unsigned long long *t, size_t tlen)
{
  size_t i;
  printf("%s", s);
  for(i=0;i<tlen-1;i++)
  {
    t[i] = t[i+1] - t[i];
  }
  printf("\n");
  printf("median: %llu\n", median(t, tlen));
  printf("average: %llu\n", average(t, tlen-1));
  printf("\n");
}

unsigned long long t[NTESTS];
unsigned char seed[32] = {0};

int main()
{
  poly ahat, ehat, ahat_shat, bhat, shat;
  poly_half f0,f1,g0,g1,g2;
  unsigned char z[2*NEWHOPE_SYMBYTES];
  unsigned char *publicseed = z;
  unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;

  randombytes(z, NEWHOPE_SYMBYTES);
  shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES);

  gen_a(&ahat, publicseed);
  poly_sample(&shat, noiseseed, 0);
  
  int i;
  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ahat, &f0, &f1);
  }

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ahat, &f0, &f1);
  }
  print_results("pt_ntt:           ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ahat, &f0, &f1);
    recover_poly(&ehat, &f0, &f1);
    poly_pt_ntt3(&shat, &g0, &g1, &g2);
    pt_ntt_bowtiemultiply(&ahat_shat, &ehat, &g0, &g1, &g2);
    poly_inv_ptntt(&ahat_shat);
  }
  print_results("mul using pt_ntt:  ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_inv_ptntt(&ahat_shat);
  }
  print_results("inv_pt_ntt:           ", t, NTESTS);
  
  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_ntt_512(&ahat);
  }
  print_results("ntt:           ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    bitrev_vector_512(&ahat);
    bitrev_vector_512(&shat);

    poly_ntt_512(&ahat);
    poly_ntt_512(&shat);
    poly_mul_pointwise(&ahat_shat, &ahat, &shat);
    poly_invntt_512(&ahat_shat);
  }
  print_results("mul using ntt:        ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_invntt_512(&ahat_shat);
  }
  print_results("inv_ntt:           ", t, NTESTS);

  return 0; 
}
