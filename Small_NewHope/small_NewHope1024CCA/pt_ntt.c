#
#include "pt_ntt.h"
#include "poly.h"
#include "ntt.h"
void poly_pt_ntt4(poly *p)
{
    poly_quarter p00, p01, p10, p11;
    split_poly(p, &p00, &p01, &p10, &p11);

    bitrev_vector_256(&p00);
    bitrev_vector_256(&p10);
    bitrev_vector_256(&p01);
    bitrev_vector_256(&p11);

    poly_ntt(&p00);
    poly_ntt(&p10);
    poly_ntt(&p01);
    poly_ntt(&p11);
    recover_poly(p, &p00, &p01, &p10, &p11);
}
void poly_pt_ntt7(poly *p, poly_quarter *p00, poly_quarter *p01, poly_quarter *p10, poly_quarter *p11, poly_quarter *p01_s, poly_quarter *p10_s, poly_quarter *p11_s)
{
    split_poly(p, p00, p01, p10, p11);
    shift_poly(p10_s,p10);
    shift_poly(p01_s,p01);
    shift_poly(p11_s,p11);

    bitrev_vector_256(p00);
    bitrev_vector_256(p10);
    bitrev_vector_256(p01);
    bitrev_vector_256(p11);

    poly_ntt(p00);
    poly_ntt(p10);
    poly_ntt(p01);
    poly_ntt(p11);

    bitrev_vector_256(p01_s);
    bitrev_vector_256(p10_s);
    bitrev_vector_256(p11_s);

    poly_ntt(p01_s);
    poly_ntt(p10_s);
    poly_ntt(p11_s);
}

void pt_ntt_bowtiemultiply(poly *b, poly *f, poly_quarter *g00, poly_quarter *g01, poly_quarter *g10, poly_quarter *g11, poly_quarter *g01_s, poly_quarter *g10_s, poly_quarter *g11_s)
{
    poly_quarter f00, f01, f10, f11, b00, b10, b01, b11, temp0, temp1, temp2, temp3;
    split_poly(f, &f00, &f01, &f10, &f11);
    // split_poly(g, &g00, &g01, &g10, &g11);
    
    poly_quarter_mul_pointwise(&temp0, &f00, g00);
    poly_quarter_mul_pointwise(&temp1, &f01, g01_s);
    poly_quarter_mul_pointwise(&temp2, &f10, g11_s);
    poly_quarter_mul_pointwise(&temp3, &f11, g10_s);
    poly_quarter_add(&b00, &temp0, &temp1, &temp2, &temp3);

    poly_quarter_mul_pointwise(&temp0, &f00, g10);
    poly_quarter_mul_pointwise(&temp1, &f01, g11_s);
    poly_quarter_mul_pointwise(&temp2, &f10, g00);
    poly_quarter_mul_pointwise(&temp3, &f11, g01_s);
    poly_quarter_add(&b10, &temp0, &temp1, &temp2, &temp3);

    poly_quarter_mul_pointwise(&temp0, &f00, g01);
    poly_quarter_mul_pointwise(&temp1, &f01, g00);
    poly_quarter_mul_pointwise(&temp2, &f10, g10);
    poly_quarter_mul_pointwise(&temp3, &f11, g11_s);
    poly_quarter_add(&b01, &temp0, &temp1, &temp2, &temp3);

    poly_quarter_mul_pointwise(&temp0, &f00, g11);
    poly_quarter_mul_pointwise(&temp1, &f01, g10);
    poly_quarter_mul_pointwise(&temp2, &f10, g01);
    poly_quarter_mul_pointwise(&temp3, &f11, g00);
    poly_quarter_add(&b11, &temp0, &temp1, &temp2, &temp3);

    recover_poly(b, &b00, &b01, &b10, &b11);
} 

void poly_inv_ptntt(poly *b)
{
    poly_quarter b00,b01,b10,b11;
    split_poly(b, &b00, &b01, &b10, &b11);
    poly_invntt(&b00);
    poly_invntt(&b10);
    poly_invntt(&b01);
    poly_invntt(&b11);

    recover_poly(b,&b00, &b01, &b10, &b11);
}

