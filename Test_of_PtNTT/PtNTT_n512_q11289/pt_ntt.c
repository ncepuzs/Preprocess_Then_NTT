#include "pt_ntt.h"
#include "poly.h"


#define NEWHOPE_N 512
#define NEWHOPE_Q 12289

void poly_pt_ntt2(poly *p, poly_half *p0, poly_half *p1)
{
    split_poly(p, p0, p1);

    bitrev_vector_256(p0);
    bitrev_vector_256(p1);

    poly_ntt_256(p0);
    poly_ntt_256(p1);
}
void poly_pt_ntt3(poly *p, poly_half *p0, poly_half *p1, poly_half *p2)
{
    split_poly(p, p0, p1);
    shift_poly(p2,p1);

    bitrev_vector_256(p0);
    bitrev_vector_256(p1);
    bitrev_vector_256(p2);

    poly_ntt_256(p0);
    poly_ntt_256(p1);
    poly_ntt_256(p2);
}

void pt_ntt_bowtiemultiply(poly *b, poly *a, poly_half *s0, poly_half *s1, poly_half *s2)
{
    poly_half a0, a1, b0, temp, b1;
    split_poly(a, &a0, &a1);
    
    poly_half_mul_pointwise(&b0, &a0, s0);
    poly_half_mul_pointwise(&temp, &a1, s2);
    poly_half_add(&b0, &temp, &b0);

    poly_half_mul_pointwise(&b1, &a1, s0);
    poly_half_mul_pointwise(&temp, &a0, s1);
    poly_half_add(&b1, &b1, &temp);

    recover_poly(b, &b0, &b1);
} 

void poly_inv_ptntt(poly *b)
{
    poly_half b0,b1;
    split_poly(b, &b0, &b1);
    poly_invntt_256(&b0);
    poly_invntt_256(&b1);

    recover_poly(b,&b0,&b1);
}

