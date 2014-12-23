/*=============================================================================

    This file is part ofELLMUL

   ELLMULis free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

   ELLMULis distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along withELLMUL if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Jean-Pierre Flori

******************************************************************************/

#ifdef T

#include "templates.h"

/*
  Weierstrass XZ ladder
  P1=P3-P2, P2, P3 |-> P4=[2]P2, P5=P2+P3
  ladd-2002-it-3
*/
void
TEMPLATE(T, weierstrass_xz_ladd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5, TEMPLATE(T, t) x4, TEMPLATE(T, t) z4,
                                 const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                                 const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                                 const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                                 const TEMPLATE(T, weierstrass_xz_t) E,
                                 const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) xx, zz, e, azz, a, b, c, d;
    TEMPLATE(T, init)(xx, K);
    TEMPLATE(T, init)(zz, K);
    TEMPLATE(T, init)(azz, K);
    TEMPLATE(T, init)(e, K);
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);

    /* XX = X2² */
    TEMPLATE(T, sqr)(xx, x2, K);
    /* ZZ = Z2² */
    TEMPLATE(T, sqr)(zz, z2, K);
    /* aZZ = a*ZZ */
    TEMPLATE(T, mul)(azz, E->a, zz, K);
    /* E = (X2+Z2)²-XX-ZZ */
    TEMPLATE(T, add)(a, x2, z2, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, add)(a, xx, zz, K);
    TEMPLATE(T, sub)(e, b, a, K);
    /* X4 = (XX-aZZ)²-b4*E*ZZ */
    TEMPLATE(T, sub)(a, xx, azz, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(a, e, zz, K);
    TEMPLATE(T, mul)(c, a, E->b4, K);
    TEMPLATE(T, sub)(x4, b, c, K);
    /* Z4 = 2*E*(XX+aZZ)+b4*ZZ² */
    TEMPLATE(T, add)(a, xx, azz, K);
    TEMPLATE(T, mul)(b, a, e, K);
    TEMPLATE(T, mul_ui)(c, b, 2, K);
    TEMPLATE(T, sqr)(a, zz, K);
    TEMPLATE(T, mul)(b, E->b4, a, K);
    TEMPLATE(T, add)(z4, c, b, K);
    /* A = X2*X3 */
    TEMPLATE(T, mul)(a, x2, x3, K);
    /* B = Z2*Z3 */
    TEMPLATE(T, mul)(b, z2, z3, K);
    /* C = X2*Z3 */
    TEMPLATE(T, mul)(c, x2, z3, K);
    /* D = X3*Z2 */
    TEMPLATE(T, mul)(d, x3, z2, K);
    /* X5 = Z1*((A-a*B)²-b4*B*(C+D)) */
    TEMPLATE(T, mul)(azz, E->a, b, K);
    TEMPLATE(T, sub)(x5, a, azz, K);
    TEMPLATE(T, sqr)(a, x5, K);
    TEMPLATE(T, add)(x5, c, d, K);
    TEMPLATE(T, mul)(e, b, x5, K);
    TEMPLATE(T, mul)(b, E->b4, e, K);
    TEMPLATE(T, sub)(e, a, b, K);
    TEMPLATE(T, mul)(x5, z1, e, K);
    /* Z5 = X1*(C-D)² */
    TEMPLATE(T, sub)(a, c, d, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(z5, x1, b, K);

    TEMPLATE(T, clear)(xx, K);
    TEMPLATE(T, clear)(zz, K);
    TEMPLATE(T, clear)(azz, K);
    TEMPLATE(T, clear)(e, K);
    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
}

#endif
