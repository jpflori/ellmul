/*=============================================================================

    This file is part of ELLMUL

    ELLMUL is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    ELLMUL is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with ELLMUL if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2014 Jean-Pierre Flori

******************************************************************************/

#ifdef T

#include "templates.h"

/*
  Weierstrass XZ doubling
  P1 |-> P3=[2]P1
  dbl-2002-bj-3
*/
void
TEMPLATE(T, weierstrass_xz_dbl)(TEMPLATE(T, t) x3, TEMPLATE(T, t) z3,
                                const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                                const TEMPLATE(T, weierstrass_xz_t) E,
                                const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) xx, zz, a, azz;
    TEMPLATE(T, init)(xx, K);
    TEMPLATE(T, init)(zz, K);
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(azz, K);

    /* XX = X1² */
    TEMPLATE(T, sqr)(xx, x1, K);
    /* ZZ = Z1² */
    TEMPLATE(T, sqr)(zz, z1, K);
    /* A = 2*((X1+Z1)²-XX-ZZ) */
    TEMPLATE(T, add)(azz, x1, z1, K);
    TEMPLATE(T, sqr)(a, azz, K);
    TEMPLATE(T, sub)(a, a, xx, K);
    TEMPLATE(T, sub)(a, a, zz, K);
    TEMPLATE(T, mul_ui)(a, a, 2, K);
    /* aZZ = a*ZZ */
    TEMPLATE(T, mul)(azz, E->a, zz, K);
    /* X3 = (XX-aZZ)²-b2*A*ZZ */
    TEMPLATE(T, sub)(x3, xx, azz, K);
    TEMPLATE(T, sqr)(x3, x3, K);
    TEMPLATE(T, mul)(z3, a, zz, K);
    TEMPLATE(T, mul)(z3, E->b2, z3, K);
    TEMPLATE(T, sub)(x3, x3, z3, K);
    /* Z3 = A*(XX+aZZ)+b4*ZZ² */
    TEMPLATE(T, add)(xx, xx, azz, K);
    TEMPLATE(T, mul)(azz, a, xx, K);
    TEMPLATE(T, sqr)(a, zz, K);
    TEMPLATE(T, mul)(zz, E->b4, a, K);
    TEMPLATE(T, add)(z3, azz, zz, K);

    TEMPLATE(T, clear)(xx, K);
    TEMPLATE(T, clear)(zz, K);
    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(azz, K);
}

#endif
