#ifndef NTT_H
#define NTT_H

#include "inttypes.h"

extern uint16_t omegas_bitrev_montgomery_128[];
extern uint16_t omegas_inv_bitrev_montgomery_128[];

extern uint16_t psis_bitrev_montgomery_128[];
extern uint16_t psis_inv_montgomery_128[];

void bitrev_vector_128(uint16_t* poly);
void mul_coefficients_128(uint16_t* poly, const uint16_t* factors);
void ntt_128(uint16_t* poly, const uint16_t* omegas);

#endif
