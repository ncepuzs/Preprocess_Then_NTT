#include "pre_ntt.h"
#include "poly.h"


#define NEWHOPE_N 512

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
int shift_poly(poly_half *f)
{
    int i;
    //coefficients in polynomials of Kyber < q = 3329 which need 2 bytes
    uint16_t the_last_coe;
    the_last_coe=f->coeffs[NEWHOPE_N/2-1]%7681;
    for(i=1;i<NEWHOPE_N/2;i++)
    {
        f->coeffs[NEWHOPE_N/2-i]=f->coeffs[NEWHOPE_N/2-i-1];
    }
    f->coeffs[0]=7681-the_last_coe;

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
int pre_ntt(poly *p, poly *f, poly *g)
{
    poly_half f0,f1,g0,g1,p0,p1,p2,p3,r1,r2;
    
    split_poly(f, &f0, &f1);
    split_poly(g, &g0, &g1);

    bitrev_vector_256(&f0);
    bitrev_vector_256(&f1);
    bitrev_vector_256(&g0);
    bitrev_vector_256(&g1);

    poly_ntt_256(&f0);
    poly_ntt_256(&f1);
    poly_ntt_256(&g0);
    poly_ntt_256(&g1);

    //pointwise_mul
    poly_half_mul_pointwise(&p0,&f0,&g0);
    poly_half_mul_pointwise(&p1,&f0,&g1);
    poly_half_mul_pointwise(&p2,&f1,&g0);
    poly_half_mul_pointwise(&p3,&f1,&g1);

    poly_invntt_256(&p0);
    // poly_invntt_512(&p1);
    // poly_invntt_512(&p2);
    poly_invntt_256(&p3);

    shift_poly(&p3);

    poly_half_add(&r1, &p0, &p3);
    poly_half_add(&r2, &p1, &p2);
    poly_invntt_256(&r2);

    recover_poly(p, &r1, &r2);

    return 1;
}