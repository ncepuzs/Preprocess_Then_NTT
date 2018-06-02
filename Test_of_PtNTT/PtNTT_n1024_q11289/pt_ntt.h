#include "poly.h"

void poly_pt_ntt2(poly *p, poly_half *p0, poly_half *p1);
void poly_pt_ntt3(poly *p, poly_half *p0, poly_half *p1, poly_half *p2);
void pt_ntt_bowtiemultiply(poly *b, poly *a, poly_half *s0, poly_half *s1, poly_half *s2);
void poly_inv_ptntt(poly *b);
