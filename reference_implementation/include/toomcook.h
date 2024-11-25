#ifndef SMAUG_TOOMCOOK_H
#define SMAUG_TOOMCOOK_H

#include "parameters.h"
#include <stdint.h>

/*#define poly_mul_acc SMAUG_NAMESPACE(poly_mul_acc)
void poly_mul_acc(const int16_t a[LWE_N], const int16_t b[LWE_N],
                  int16_t res[LWE_N]);*/

#define evaluation_single SMAUG_NAMESPACE(evaluation_single)
void evaluation_single(const uint16_t *b, uint16_t sw[7][64]);

#define toom_cook_4way_evaluate SMAUG_NAMESPACE(toom_cook_4way_evaluate)
void toom_cook_4way_evaluate(const uint16_t *a1, const uint16_t sw[7][64], uint16_t *res);

#define toom_cook_4way_interpol SMAUG_NAMESPACE(toom_cook_4way_interpol)
void toom_cook_4way_interpol(const uint16_t *a_i, uint16_t *res_i);

#endif // SMAUG_TOOMCOOK_H
