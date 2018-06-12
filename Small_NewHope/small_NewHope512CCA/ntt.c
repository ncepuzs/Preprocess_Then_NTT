#include "inttypes.h"
#include "ntt.h"
#include "params.h"
#include "reduce.h"

/************************************************************
* Name:        bitrev_table
*
* Description: Contains bit-reversed 10-bit indices to be used to re-order 
*              polynomials before number theoratic transform 
************************************************************/
static uint16_t bitrev_table_128[128] = {
0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104,24, 88, 56, 120,
4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121, 
5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123, 
7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127
};


/*************************************************
* Name:        bitrev_vector
* 
* Description: Permutes coefficients of a polynomial into bitreversed order
*
* Arguments:   - uint16_t* poly: pointer to in/output polynomial
**************************************************/
void bitrev_vector_128(uint16_t* poly)
{
    unsigned int i,r;
    uint16_t tmp;

    for(i = 0; i < 128; i++)
    {
        r = bitrev_table_128[i];
        if (i < r)
        {
          tmp = poly[i];
          poly[i] = poly[r];
          poly[r] = tmp;
        }
    }
}
/*************************************************
* Name:        mul_coefficients
* 
* Description: Performs pointwise (coefficient-wise) multiplication
*              of two polynomials
* Arguments:   - uint16_t* poly:          pointer to in/output polynomial
*              - const uint16_t* factors: pointer to input polynomial, coefficients 
*                                         are assumed to be in Montgomery representation
**************************************************/
void mul_coefficients_128(uint16_t* poly, const uint16_t* factors)
{
    unsigned int i;

    for(i = 0; i < 128; i++)
      poly[i] = montgomery_reduce((poly[i] * factors[i]));
}

void /*************************************************
* Name:        ntt_256
* 
* Description: Computes number-theoretic transform (NTT) of
*              a polynomial in place; inputs assumed to be in
*              bitreversed order, output in normal order
*
* Arguments:   - uint16_t * a:          pointer to in/output polynomial
*              - const uint16_t* omega: pointer to input powers of root of unity omega;
*                                       assumed to be in Montgomery domain
**************************************************/
ntt_128(uint16_t * a, const uint16_t* omega)
{
  int i, start, j, jTwiddle, distance;
  uint16_t temp, W;


  for(i=0;i<7;i+=2)
  {
    // Even level
    distance = (1<<i);
    for(start = 0; start < distance;start++)
    {
      jTwiddle = 0;
      for(j=start;j<127;j+=2*distance)
      {
        W = omega[jTwiddle++];
        temp = a[j];
        a[j] = (temp + a[j + distance]); // Omit reduction (be lazy)
        a[j + distance] = montgomery_reduce((W * ((uint32_t)temp + 3*7681 - a[j + distance])));
      }
    }
    if(i+1<7){
      // Odd level
      distance <<= 1;
      for(start = 0; start < distance;start++)
      {
        jTwiddle = 0;
        for(j=start;j<127;j+=2*distance)
        {
          W = omega[jTwiddle++];
          temp = a[j];
          a[j] = (temp + a[j + distance]) % 7681;
          a[j + distance] = montgomery_reduce((W * ((uint32_t)temp + 3*7681 - a[j + distance])));
        }
      }
    }
  }
}
