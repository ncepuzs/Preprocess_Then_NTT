#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "precomp.c"
#include <stdlib.h>
#include "pre_ntt.h"

#define NEWHOPE_Q 12289
#define NEWHOPE_N 512
#define NEWHOPE_SYMBYTES 32

uint16_t positive_mod(int32_t r)
{
    if(r<0)
        while(r<0)
        {
            r=r+12289;
        }
    return r;
}
void normal_poly_mul(poly *r,poly *a,poly *b)
{
    int i,j1,j2;
    for(i=0;i<NEWHOPE_N;i++)
            r->coeffs[i]=0;
    
    for(j1=0;j1<NEWHOPE_N;j1++)
        for(j2=0;j2<NEWHOPE_N;j2++)
        {
            if(j1+j2>=NEWHOPE_N)
            {
                i=j1+j2-NEWHOPE_N;
                int32_t temp=r->coeffs[i]-(a->coeffs[j1]*b->coeffs[j2]);
                temp=positive_mod(temp);
                r->coeffs[i]=(temp)%12289;
            } 
            else
            {
                i=j1+j2;
                r->coeffs[i]=(r->coeffs[i]+(a->coeffs[j1]*b->coeffs[j2]))%12289;
            }
        }
}


int main()
{
    int i;
    poly poly_f, poly_g, poly_p1, poly_p2, poly_p3;

    unsigned char z[2*NEWHOPE_SYMBYTES];
    unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;
    randombytes(z, NEWHOPE_SYMBYTES);
    shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES);

    poly_sample(&poly_f, noiseseed, 0);
    poly_sample(&poly_g, noiseseed, 0);
 
    normal_poly_mul(&poly_p1, &poly_f, &poly_g);

    pre_ntt(&poly_p2, &poly_f, &poly_g);

    bitrev_vector_512(&poly_f);
    bitrev_vector_512(&poly_g);

    poly_ntt_512(&poly_f);
    poly_ntt_512(&poly_g);
    poly_mul_pointwise(&poly_p3, &poly_f, &poly_g);
    poly_invntt_512(&poly_p3);
    
    if(poly_equal(poly_p1, poly_p2)&&poly_equal(poly_p1, poly_p3))
    {
        printf("The result of pre_ntt is correct!\n");
    }
   
    else 
    {
        printf("The result is wrong!\n");
        printf("-----------\n");
        print_poly(&poly_p1);
        printf("--------------\n");
        print_poly(&poly_p2);
        printf("--------------\n");
        print_poly(&poly_p3);
    }

    return 1;
}
