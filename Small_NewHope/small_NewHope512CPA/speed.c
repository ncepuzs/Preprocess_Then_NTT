#include "poly.h"
#include "pt_ntt.h"
#include "cpucycles.h"
#include <stdlib.h>
#include <stdio.h>
#include "params.h"
#include "api.h"

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
  poly_quarter  poly_f00, poly_f01, poly_f10, poly_f11, poly_g00, poly_g01, poly_g10, poly_g11, poly_g01_s, poly_g10_s, poly_g11_s;
  unsigned char z[2*NEWHOPE_SYMBYTES];
  unsigned char *publicseed = z;
  unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;
  unsigned char *r[2*NEWHOPE_POLYBYTES];
  unsigned char sk[NEWHOPE_POLYBYTES],pk[NEWHOPE_POLYBYTES],ct[CRYPTO_CIPHERTEXTBYTES], ss[CRYPTO_BYTES], ss1[CRYPTO_BYTES];
  randombytes(z, NEWHOPE_SYMBYTES);
  shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES);

  gen_a(&ahat, publicseed);
  poly_sample(&shat, noiseseed, 0);
  
  int i;
  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    gen_a(&ahat, publicseed);
  }
  print_results("gen_a:            ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_sample(&shat, noiseseed, 0);
  }
  print_results("gen_s:            ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt4(&ahat, &poly_f00, &poly_f01, &poly_f10, &poly_f11);
  }
  print_results("ptntt:            ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_inv_ptntt(&ahat);
  }
  print_results("inv_ptntt:         ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    poly_pt_ntt4(&ahat, &poly_f00, &poly_f01, &poly_f10, &poly_f11);
    recover_poly(&ahat, &poly_f00, &poly_f01, &poly_f10, &poly_f11);
    poly_pt_ntt7(&shat, &poly_g00, &poly_g01, &poly_g10, &poly_g11, &poly_g01_s, &poly_g10_s, &poly_g11_s);
    pt_ntt_bowtiemultiply(&ahat, &ahat, &poly_g00, &poly_g01, &poly_g10, &poly_g11, &poly_g01_s, &poly_g10_s, &poly_g11_s);
    poly_inv_ptntt(&ahat);
  }
  print_results("total_ptntt:      ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_kem_keypair(pk,sk);
  }
  print_results("kepair:           ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_kem_enc(ct, ss, pk);
  }
  print_results("enc:              ", t, NTESTS);

  for(i=0; i<NTESTS; i++)
  {
    t[i] = cpucycles();
    crypto_kem_dec(ss1, ct, sk);
  }
  print_results("dec:              ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   decode_c(&ahat, &shat, ct);
  // }
  // print_results("decode_c:        ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   poly_half_tobytes(&s0, &s1, &s2, sk);

  // }
  // print_results("tobytes:        ", t, NTESTS);  

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   poly_half_frombytes(&s0, &s1, &s2, sk);

  // }
  // print_results("frombytes:        ", t, NTESTS);

  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   pt_ntt_bowtiemultiply(&bhat, &ahat, &s0, &s1, &s2);
  // }
  // print_results("mul:        ", t, NTESTS);
  
  // for(i=0; i<NTESTS; i++)
  // {
  //   t[i] = cpucycles();
  //   poly_tomsg(ss1, &bhat);
  // }
  // print_results("tomsg:        ", t, NTESTS);

  return 0; 
}

