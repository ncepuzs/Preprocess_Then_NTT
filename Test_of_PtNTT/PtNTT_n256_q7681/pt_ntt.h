#include "poly.h"
#include "params.h"

int poly_pt_ntt2(poly *a, poly_half *a0, poly_half *a1);
int poly_pt_ntt3(poly *s, poly_half *s0, poly_half *s1, poly_half *s2);
int pt_ntt_bowtiemultiply(poly *p, poly *a, poly_half *s0, poly_half *s1, poly_half *s2);
int poly_inv_ptntt(poly *a);