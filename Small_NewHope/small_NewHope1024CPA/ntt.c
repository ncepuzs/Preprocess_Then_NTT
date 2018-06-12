#include "inttypes.h"
#include "ntt.h"
#include "params.h"
#include "reduce.h"

#if (NEWHOPE_N == 512)
/************************************************************
* Name:        bitrev_table
*
* Description: Contains bit-reversed 9-bit indices to be used to re-order 
*              polynomials before number theoratic transform 
************************************************************/
static uint16_t bitrev_table[128] = {
0, 64, 32, 96, 16, 80, 48, 112, 8, 72, 40, 104,24, 88, 56, 120,
4, 68, 36, 100, 20, 84, 52, 116, 12, 76, 44, 108, 28, 92, 60, 124,
2, 66, 34, 98, 18, 82, 50, 114, 10, 74, 42, 106, 26, 90, 58, 122,
6, 70, 38, 102, 22, 86, 54, 118, 14, 78, 46, 110, 30, 94, 62, 126,
1, 65, 33, 97, 17, 81, 49, 113, 9, 73, 41, 105, 25, 89, 57, 121, 
5, 69, 37, 101, 21, 85, 53, 117, 13, 77, 45, 109, 29, 93, 61, 125,
3, 67, 35, 99, 19, 83, 51, 115, 11, 75, 43, 107, 27, 91, 59, 123, 
7, 71, 39, 103, 23, 87, 55, 119, 15, 79, 47, 111, 31, 95, 63, 127
};

#elif (NEWHOPE_N == 1024)

/************************************************************
* Name:        bitrev_table
*
* Description: Contains bit-reversed 9-bit indices to be used to re-order 
*              polynomials before number theoratic transform 
************************************************************/
static uint16_t bitrev_table_256 [256] = {
0,128,64,192,32,160,96,224,16,144,80,208,48,176,112,240, 
     8,136,72,200,40,168,104,232,24,152,88,216,56,184,120,248, 
     4,132,68,196,36,164,100,228,20,148,84,212,52,180,116,244, 
     12,140,76,204,44,172,108,236,28,156,92,220,60,188,124,252, 
     2,130,66,194,34,162,98,226,18,146,82,210,50,178,114,242, 
     10,138,74,202,42,170,106,234,26,154,90,218,58,186,122,250, 
     6,134,70,198,38,166,102,230,22,150,86,214,54,182,118,246, 
     14,142,78,206,46,174,110,238,30,158,94,222,62,190,126,254, 
     1,129,65,193,33,161,97,225,17,145,81,209,49,177,113,241, 
     9,137,73,201,41,169,105,233,25,153,89,217,57,185,121,249, 
     5,133,69,197,37,165,101,229,21,149,85,213,53,181,117,245, 
     13,141,77,205,45,173,109,237,29,157,93,221,61,189,125,253, 
     3,131,67,195,35,163,99,227,19,147,83,211,51,179,115,243, 
     11,139,75,203,43,171,107,235,27,155,91,219,59,187,123,251, 
     7,135,71,199,39,167,103,231,23,151,87,215,55,183,119,247, 
     15,143,79,207,47,175,111,239,31,159,95,223,63,191,127,255};


#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif


/*************************************************
* Name:        bitrev_vector
* 
* Description: Permutes coefficients of a polynomial into bitreversed order
*
* Arguments:   - uint16_t* poly: pointer to in/output polynomial
**************************************************/
void bitrev_vector_256(uint16_t* poly)
{
    unsigned int i,r;
    uint16_t tmp;

    for(i = 0; i < 256; i++)
    {
        r = bitrev_table_256[i];
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
void mul_coefficients_256(uint16_t* poly, const uint16_t* factors)
{
    unsigned int i;

    for(i = 0; i < 256; i++)
      poly[i] = montgomery_reduce((poly[i] * factors[i]));
}


#if (NEWHOPE_N == 512)

/*************************************************
* Name:        ntt
* 
* Description: Computes number-theoretic transform (NTT) of
*              a polynomial in place; inputs assumed to be in
*              bitreversed order, output in normal order
*
* Arguments:   - uint16_t * a:          pointer to in/output polynomial
*              - const uint16_t* omega: pointer to input powers of root of unity omega;
*                                       assumed to be in Montgomery domain
**************************************************/
void ntt(uint16_t * a, const uint16_t* omega)
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

#elif (NEWHOPE_N == 1024)

/*************************************************
* Name:        ntt
* 
* Description: Computes number-theoretic transform (NTT) of
*              a polynomial in place; inputs assumed to be in
*              bitreversed order, output in normal order
*
* Arguments:   - uint16_t * a:          pointer to in/output polynomial
*              - const uint16_t* omega: pointer to input powers of root of unity omega;
*                                       assumed to be in Montgomery domain
**************************************************/
void ntt_256(uint16_t * a, const uint16_t* omega)
{
  int i, start, j, jTwiddle, distance;
  uint16_t temp, W;


  for(i=0;i<8;i+=2)
  {
    // Even level
    distance = (1<<i);
    for(start = 0; start < distance;start++)
    {
      jTwiddle = 0;
      for(j=start;j<255;j+=2*distance)
      {
        W = omega[jTwiddle++];
        temp = a[j];
        a[j] = (temp + a[j + distance]); // Omit reduction (be lazy)
        a[j + distance] = montgomery_reduce((W * ((uint32_t)temp + 3*7681 - a[j + distance])));
      }
    }

    // Odd level
    distance <<= 1;
    for(start = 0; start < distance;start++)
    {
      jTwiddle = 0;
      for(j=start;j<255;j+=2*distance)
      {
        W = omega[jTwiddle++];
        temp = a[j];
        a[j] = (temp + a[j + distance]) % 7681;
        a[j + distance] = montgomery_reduce((W * ((uint32_t)temp + 3*7681 - a[j + distance])));
      }
    }
  }
}



#else
#error "NEWHOPE_N must be either 512 or 1024"
#endif
