#ifndef TiMER_TOOMCOOK_H
#define TiMER_TOOMCOOK_H

#include "parameters.h"
#include <stdint.h>

/*#define poly_mul_acc TiMER_NAMESPACE(poly_mul_acc)
void poly_mul_acc(const int16_t a[LWE_N], const int16_t b[LWE_N],
                  int16_t res[LWE_N]);*/

#define evaluation_single TiMER_NAMESPACE(evaluation_single)
void evaluation_single(const uint16_t *b, uint16_t sw[7][64]);

#define toom_cook_4way_evaluate TiMER_NAMESPACE(toom_cook_4way_evaluate)
void toom_cook_4way_evaluate(const uint16_t *a1, const uint16_t sw[7][64], uint16_t *res);

#define toom_cook_4way_interpol TiMER_NAMESPACE(toom_cook_4way_interpol)
void toom_cook_4way_interpol(const uint16_t *a_i, uint16_t *res_i);

#endif // TiMER_TOOMCOOK_H
