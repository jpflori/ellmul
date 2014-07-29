/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Jean-Pierre Flori

******************************************************************************/

#ifdef T

#include "templates.h"

/*
  Weierstrass XZ differential addition
  P1=P3-P2, P2, P3 |-> P5=P2+P3
  dadd-2002-it-3
*/
void
TEMPLATE(T, weierstrass_xz_dadd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5,
                                 const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                                 const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                                 const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                                 const TEMPLATE(T, weierstrass_xz_t) E,
                                 const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, c, d, e, f;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);
    TEMPLATE(T, init)(e, K);
    TEMPLATE(T, init)(f, K);

    /* A = X2*X3 */
    TEMPLATE(T, mul)(a, x2, x3, K);
    /* B = Z2*Z3 */
    TEMPLATE(T, mul)(b, z2, z3, K);
    /* C = X2*Z3 */
    TEMPLATE(T, mul)(c, x2, z3, K);
    /* D = X3*Z2 */
    TEMPLATE(T, mul)(d, x3, z2, K);
    /* E = (A-a*B)² */
    TEMPLATE(T, mul)(e, E->a, b, K);
    TEMPLATE(T, sub)(f, a, e, K);
    TEMPLATE(T, sqr)(e, f, K);
    /* F = 4*b*B*(C+D) */
    TEMPLATE(T, add)(f, c, d, K);
    TEMPLATE(T, mul)(a, b, f, K);
    TEMPLATE(T, mul)(b, E->b, a, K);
    TEMPLATE(T, mul_ui)(f, b, 4, K);
    /* X5 = Z1*(E-F) */
    TEMPLATE(T, sub)(a, e, f, K);
    TEMPLATE(T, mul)(x5, z1, a, K);
    /* Z5 = X1*(C-D)² */
    TEMPLATE(T, sub)(a, c, d, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(z5, x1, b, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
    TEMPLATE(T, clear)(e, K);
    TEMPLATE(T, clear)(f, K);
}

#endif
