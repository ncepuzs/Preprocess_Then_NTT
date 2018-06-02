#include "poly.h"
#include "ntt.h"
#include "reduce.h"
#include "precomp.c"
#include <stdlib.h>
#include "pt_ntt.h"

#define NEWHOPE_Q 12289
#define NEWHOPE_N 1024
#define NEWHOPE_SYMBYTES 32

uint16_t positive_mod(int32_t r)
{
    if(r<0)
        while(r<0)
        {
            r=r+NEWHOPE_Q;
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
                r->coeffs[i]=(temp)%NEWHOPE_Q;
            } 
            else
            {
                i=j1+j2;
                r->coeffs[i]=(r->coeffs[i]+(a->coeffs[j1]*b->coeffs[j2]))%NEWHOPE_Q;
            }
        }
}


int main()
{
    int i;
    poly poly_f, poly_g, poly_p1, poly_p2, poly_p3,poly_fprime;
    poly_half  poly_f0, poly_f1, poly_g0, poly_g1, poly_g2;

    unsigned char z[2*NEWHOPE_SYMBYTES];
    unsigned char *noiseseed = z+NEWHOPE_SYMBYTES;
    randombytes(z, NEWHOPE_SYMBYTES);
    shake256(z, 2*NEWHOPE_SYMBYTES, z, NEWHOPE_SYMBYTES);

    poly_sample(&poly_f, noiseseed, 0);
    poly_sample(&poly_g, noiseseed, 0);

    //multiply f and g normally
    normal_poly_mul(&poly_p1, &poly_f, &poly_g);

    //apply pt_ntt to multiply f and g
    poly_pt_ntt2(&poly_f, &poly_f0, &poly_f1);
    recover_poly(&poly_fprime, &poly_f0, &poly_f1);
    poly_pt_ntt3(&poly_g, &poly_g0, &poly_g1, &poly_g2);
    pt_ntt_bowtiemultiply(&poly_p2, &poly_fprime, &poly_g0, &poly_g1, &poly_g2);
    poly_inv_ptntt(&poly_p2);
    
    //apply ntt to multiply f and g
    bitrev_vector_1024(&poly_f);
    bitrev_vector_1024(&poly_g);

    poly_ntt_1024(&poly_f);
    poly_ntt_1024(&poly_g);
    poly_mul_pointwise(&poly_p3, &poly_f, &poly_g);
    poly_invntt_1024(&poly_p3);
    
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
