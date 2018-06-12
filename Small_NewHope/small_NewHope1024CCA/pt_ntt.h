#ifndef PT_NTT_H
#define PT_NTT_H

#include "poly.h"

void poly_pt_ntt4(poly *p);
void poly_pt_ntt7(poly *p, poly_quarter *p00, poly_quarter *p10, poly_quarter *p01, poly_quarter *p11, poly_quarter *p10_s, poly_quarter *p01_s, poly_quarter *p11_s);
void pt_ntt_bowtiemultiply(poly *b, poly *f, poly_quarter *g00, poly_quarter *g01, poly_quarter *g10, poly_quarter *g11, poly_quarter *g01_s, poly_quarter *g10_s, poly_quarter *g11_s);
void poly_inv_ptntt(poly *b);


#endif