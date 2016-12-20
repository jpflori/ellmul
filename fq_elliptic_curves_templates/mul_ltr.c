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
#ifdef EC

#include "templates.h"

/*
  Elliptic curve point multiplication
  Left-to-right binary exponentiation
*/
void
TEMPLATE3(T, EC, mul_ltr)(TEMPLATE(T, t) x, TEMPLATE(T, t) z,
                          const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                          const fmpz_t m,
                          const TEMPLATE3(T, EC, t) E,
                          const TEMPLATE(T, ctx_t) K)
{
    if (fmpz_is_zero(m))
    {
        TEMPLATE(T, one)(x, K);
        TEMPLATE(T, zero)(z, K);
    }
    else if (fmpz_is_one(m))
    {
        TEMPLATE(T, set)(x, x1, K);
        TEMPLATE(T, set)(z, z1, K);
    }
    else
    {
        ulong bit;
        TEMPLATE(T, t) x2, z2, x3, z3, x4, z4, x5, z5;
        TEMPLATE(T, init)(x2, K);
        TEMPLATE(T, init)(z2, K);
        TEMPLATE(T, init)(x3, K);
        TEMPLATE(T, init)(z3, K);
        TEMPLATE(T, init)(x4, K);
        TEMPLATE(T, init)(z4, K);
        TEMPLATE(T, init)(x5, K);
        TEMPLATE(T, init)(z5, K);

        /* P */
        TEMPLATE(T, set)(x2, x1, K);
        TEMPLATE(T, set)(z2, z1, K);
        /* 2*P */
        TEMPLATE3(T, EC, dbl)(x3, z3, x1, z1, E, K);

        for (bit = fmpz_bits(m) - 2; bit > 0; bit--)
        {
            if (fmpz_tstbit(m, bit))
            {
                /* P, [k+1]P, [k]P -> [2k+2]P, [2k+1]P */
                TEMPLATE3(T, EC, ladd)(x5, z5, x4, z4, x2, z2, x3, z3, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x5, K);
                TEMPLATE(T, swap)(z2, z5, K);
                TEMPLATE(T, swap)(x3, x4, K);
                TEMPLATE(T, swap)(z3, z4, K);
            }
            else
            {
                /* P, [k]P, [k+1]P -> [2k]P, [2k+1]P */
                TEMPLATE3(T, EC, ladd)(x5, z5, x4, z4, x3, z3, x2, z2, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x4, K);
                TEMPLATE(T, swap)(z2, z4, K);
                TEMPLATE(T, swap)(x3, x5, K);
                TEMPLATE(T, swap)(z3, z5, K);
            }
        }

        /* Last iteration, no need for a ladder */
        {
            if (fmpz_tstbit(m, 0))
            {
                TEMPLATE3(T, EC, dadd)(x, z, x3, z3, x2, z2, x1, z1, E, K);
            }
            else
            {
                TEMPLATE3(T, EC, dbl)(x, z, x2, z2, E, K);
            }
        }

        TEMPLATE(T, clear)(x2, K);
        TEMPLATE(T, clear)(z2, K);
        TEMPLATE(T, clear)(x3, K);
        TEMPLATE(T, clear)(z3, K);
        TEMPLATE(T, clear)(x4, K);
        TEMPLATE(T, clear)(z4, K);
        TEMPLATE(T, clear)(x5, K);
        TEMPLATE(T, clear)(z5, K);
    }
}

#endif
#endif
