#include "poly.h"
#include "cpucycles.h"
#include <stdlib.h>
#include <stdio.h>

#define NTESTS 10000

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
  unsigned char buf[64];
  unsigned char *publicseed = buf;
  unsigned char *noiseseed = buf+KYBER_SYMBYTES;
  int i;
  unsigned char nonce=0;

  randombytes(buf, 32);
  sha3_512(buf, buf, 32);

  poly ab1, ab2, ap, bp;
  poly_getnoise(&ap, noiseseed, nonce++);
  poly_getnoise(&bp, noiseseed, nonce++);
  poly_half f0,f1,g0,g1,g2;

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ap, &f0, &f1);
  }

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ap, &f0, &f1);
  }
  print_results("pt_ntt:           ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   poly_pt_ntt2(&bp, &g0, &g1, &g2);
  // }
  // print_results("pt_ntt3:           ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   poly_ntt_128(&f0);
  // }
  // print_results("poly_ntt_128:           ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   split_poly(&ap,&f1,&f0);
  // }
  // print_results("split:           ", t, NTESTS);
  
  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt2(&ap, &f0, &f1);
    recover_poly(&ab1, &f0, &f1);
    poly_pt_ntt3(&bp, &g0, &g1, &g2);
    pt_ntt_bowtiemultiply(&ab1, &ab1, &g0, &g1, &g2);
    poly_inv_ptntt(&ab1);
  }
  print_results("mul using pt_ntt:  ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_inv_ptntt(&ab1, &f0, &f1);
  }
  print_results("inv_pt_ntt:           ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_ntt(&ap);
  }
  print_results("ntt:           ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    ntt(&ap);
    ntt(&bp);
    poly_pointwise_mul(&ab2, &ap, &bp);
    invntt(&ab2);
  }
  print_results("mul using ntt:        ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    invntt(&ab2);
  }
  print_results("inv_ntt:        ", t, NTESTS);

  return 0;
}
