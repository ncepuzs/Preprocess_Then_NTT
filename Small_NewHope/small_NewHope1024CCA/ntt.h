#ifndef NTT_H
#define NTT_H

#include "inttypes.h"

extern uint16_t omegas_bitrev_montgomery[];
extern uint16_t omegas_inv_bitrev_montgomery[];

extern uint16_t psis_bitrev_montgomery[];
extern uint16_t psis_inv_montgomery[];


void ntt_256(uint16_t * a, const uint16_t* omega);
void bitrev_vector_256(uint16_t* poly);
void mul_coefficients_256(uint16_t* poly, const uint16_t* factors);

#endif
