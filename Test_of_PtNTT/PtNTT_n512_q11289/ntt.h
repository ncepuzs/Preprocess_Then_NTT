#ifndef NTT_H
#define NTT_H

#include "inttypes.h"

extern uint16_t omegas_bitrev_montgomery_512[];
extern uint16_t omegas_inv_bitrev_montgomery_512[];

extern uint16_t psis_bitrev_montgomery_512[];
extern uint16_t psis_inv_montgomery_512[];

extern uint16_t omegas_bitrev_montgomery_256[];
extern uint16_t omegas_inv_bitrev_montgomery_256[];

extern uint16_t psis_bitrev_montgomery_256[];
extern uint16_t psis_inv_montgomery_256[];

void bitrev_vector_512(uint16_t* poly);
void bitrev_vector_256(uint16_t* poly);
void mul_coefficients_512(uint16_t* poly, const uint16_t* factors);
void mul_coefficients_256(uint16_t* poly, const uint16_t* factors);
void ntt_512(uint16_t* poly, const uint16_t* omegas);
void ntt_256(uint16_t* poly, const uint16_t* omegas);

#endif
