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
  Montgomery XZ differential addition
  P1=P3-P2, P2, P3 |-> P5=P2+P3
  dadd-1987-m-3
*/
void
TEMPLATE(T, montgomery_xz_dadd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5,
                                const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                                const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                                const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                                const TEMPLATE(T, montgomery_xz_t) E,
                                const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, c, d, da, cb;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);
    TEMPLATE(T, init)(da, K);
    TEMPLATE(T, init)(cb, K);

    /* A = X2+Z2 */
    TEMPLATE(T, add)(a, x2, z2, K);
    /* B = X2-Z2 */
    TEMPLATE(T, sub)(b, x2, z2, K);
    /* C = X3+Z3 */
    TEMPLATE(T, add)(c, x3, z3, K);
    /* D = X3-Z3 */
    TEMPLATE(T, sub)(d, x3, z3, K);
    /* DA = D*A */
    TEMPLATE(T, mul)(da, d, a, K);
    /* CB = C*B */
    TEMPLATE(T, mul)(cb, c, b, K);
    /* X5 = Z1*(DA+CB)² */
    TEMPLATE(T, add)(x5, da, cb, K);
    TEMPLATE(T, sqr)(a, x5, K);
    TEMPLATE(T, mul)(x5, z1, a, K);
    /* Z5 = X1*(DA-CB)² */
    TEMPLATE(T, sub)(z5, da, cb, K);
    TEMPLATE(T, sqr)(a, z5, K);
    TEMPLATE(T, mul)(z5, x1, a, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
    TEMPLATE(T, clear)(da, K);
    TEMPLATE(T, clear)(cb, K);
}

#endif
