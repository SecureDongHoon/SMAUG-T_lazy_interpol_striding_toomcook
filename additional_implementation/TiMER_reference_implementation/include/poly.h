#ifndef TiMER_POLY_H
#define TiMER_POLY_H

#include "parameters.h"
#include "toomcook.h"

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

typedef struct {
    int16_t coeffs[LWE_N];
} poly;

typedef struct {
    poly vec[MODULE_RANK];
} polyvec;

#define vec_vec_mult_add TiMER_NAMESPACE(vec_vec_mult_add)
void vec_vec_mult_add(poly *r, const polyvec *a, uint16_t sw[MODULE_RANK][7][64],
                      const uint8_t mod);
#define matrix_vec_mult_add TiMER_NAMESPACE(matrix_vec_mult_add)
void matrix_vec_mult_add(polyvec *r, const polyvec a[MODULE_RANK],
                         uint16_t sw[MODULE_RANK][7][64]);
#define matrix_vec_mult_sub TiMER_NAMESPACE(matrix_vec_mult_sub)
void matrix_vec_mult_sub(polyvec *r, const polyvec a[MODULE_RANK],
                         uint16_t sw[MODULE_RANK][7][64]);
                         
uint16_t flipabs(uint16_t x);
void poly_frommsg(poly *r, const uint8_t *msg);
void poly_tomsg(unsigned char *msg, const poly *x);


#endif // TiMER_POLY_H
