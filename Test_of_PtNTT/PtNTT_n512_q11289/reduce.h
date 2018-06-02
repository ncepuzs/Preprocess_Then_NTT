#ifndef REDUCE_H
#define REDUCE_H

#include <stdint.h>

uint16_t montgomery_reduce(uint32_t a);
uint16_t montgomery_reduce_256(uint32_t a);

#endif
