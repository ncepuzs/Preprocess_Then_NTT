#include "poly.h"
#include "randombytes.h"
#include "params.h"
#include "pt_ntt.h"

uint16_t posi_mod(int32_t r)
{
    if(r<0)
        while(r<0)
        {
            r=r+KYBER_Q;
        }
    return r;
}

void normal_poly_mul(poly *r,poly *a,poly *b)
{
    int i,j1,j2;
    for(i=0;i<KYBER_N;i++)
            r->coeffs[i]=0;
    
    for(j1=0;j1<KYBER_N;j1++)
        for(j2=0;j2<KYBER_N;j2++)
        {
            if(j1+j2>=KYBER_N)
            {
                i=j1+j2-KYBER_N;
                int32_t temp=r->coeffs[i]-(a->coeffs[j1]*b->coeffs[j2]);
                temp=posi_mod(temp);
                r->coeffs[i]=(temp)%KYBER_Q;
            } 
            else
            {
                i=j1+j2;
                r->coeffs[i]=(r->coeffs[i]+(a->coeffs[j1]*b->coeffs[j2]))%KYBER_Q;
            }
        }
}

int main(void)
{
    unsigned char buf[64];
    unsigned char *publicseed = buf;
    unsigned char *noiseseed = buf+KYBER_SYMBYTES;
    int i;
    unsigned char nonce=0;

    randombytes(buf, 32);
    sha3_512(buf, buf, 32);

    poly poly_f, poly_g, poly_p1, poly_p2, poly_p3, f_prime;
    poly_half f0,f1, g0,g1,g2;

    poly_getnoise(&poly_f, noiseseed, nonce++);
    poly_getnoise(&poly_g, noiseseed, nonce++);

    //multiply f and g normally
    normal_poly_mul(&poly_p1, &poly_f, &poly_g);

    //apply pt_ntt to multiply f and g
    poly_pt_ntt2(&poly_f, &f0, &f1);
    recover_poly(&f_prime, &f0, &f1);
    poly_pt_ntt3(&poly_g, &g0, &g1, &g2);
    pt_ntt_bowtiemultiply(&poly_p2, &f_prime, &g0, &g1, &g2);
    poly_inv_ptntt(&poly_p2);
    
    // ntt(&poly_f);
    // ntt(&poly_g);
    // poly_pointwise_mul(&poly_p3, &poly_f, &poly_g);
    // invntt(&poly_p3);


    if(poly_equal(poly_p1, poly_p2))
    {
        printf("The result of pre_ntt is correct!\n");
    }
    else 
    {
        printf("The result of pre_ntt is wrong!\n");
        printf("\n");
        print_poly(&poly_p1);
        printf("\n");
        print_poly(&poly_p2);

    }

  return 0;
}
