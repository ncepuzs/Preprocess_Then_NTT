#include "pt_ntt.h"
#include "poly.h"
#include "ntt.h"

int poly_pt_ntt2(poly *a, poly_half *a0, poly_half *a1)
{
    split_poly(a, a0, a1);

    poly_ntt(a0);
    poly_ntt(a1);

    return 1;
}
int poly_pt_ntt3(poly *s, poly_half *s0, poly_half *s1, poly_half *s2)
{
    split_poly(s, s0, s1);
    shift_poly(s2, s1);

    poly_ntt(s0);
    poly_ntt(s1);
    poly_ntt(s2);

    return 1;

}
int pt_ntt_bowtiemultiply(poly *p, poly *a, poly_half *s0, poly_half *s1, poly_half *s2)
{
    poly_half a0, a1,r1,r2,r3,r4,p0,p1;
    split_poly(a, &a0, &a1);

    //even elements
    poly_half_mul_pointwise(&r1, &a0, s0);
    poly_half_mul_pointwise(&r2, &a1, s2);
    poly_half_add(&p0, &r1, &r2);

    //odd elements
    poly_half_mul_pointwise(&r3, &a1, s0);
    poly_half_mul_pointwise(&r4, &a0, s1);
    poly_half_add(&p1, &r3, &r4);

    recover_poly(p, &p0, &p1);

    return 1;
}
int poly_inv_ptntt(poly *a)
{
    poly_half a0, a1;
    split_poly(a, &a0, &a1);
    poly_invntt(&a0);
    poly_invntt(&a1);
    recover_poly(a, &a0, &a1);

    return 1;
}
