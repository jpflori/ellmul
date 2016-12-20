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
  Montgomery XZ doubling
  P1 |-> P3=[2]P1
  dbl-1987-m-3
*/
void
TEMPLATE(T, montgomery_xz_dbl)(TEMPLATE(T, t) x3, TEMPLATE(T, t) z3,
                               const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                               const TEMPLATE(T, montgomery_xz_t) E,
                               const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, aa, bb;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(aa, K);
    TEMPLATE(T, init)(bb, K);

    /* A = X1+Z1 */
    TEMPLATE(T, add)(a, x1, z1, K);
    /* AA = A² */
    TEMPLATE(T, sqr)(aa, a, K);
    /* B = X1-Z1 */
    TEMPLATE(T, sub)(b, x1, z1, K);
    /* BB = B² */
    TEMPLATE(T, sqr)(bb, b, K);
    /* C = AA-BB */
    TEMPLATE(T, sub)(x3, aa, bb, K);
    /* Z3 = C*(BB+a24*C) */
    TEMPLATE(T, mul)(a, E->a24, x3, K);
    TEMPLATE(T, add)(b, bb, a, K);
    TEMPLATE(T, mul)(z3, x3, b, K);
    /* X3 = AA*BB */
    TEMPLATE(T, mul)(x3, aa, bb, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(aa, K);
    TEMPLATE(T, clear)(bb, K);
}

#endif
