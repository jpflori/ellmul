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

void
TEMPLATE3(T, EC, init)(TEMPLATE3(T, EC, t) E, TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, EC, clear)(TEMPLATE3(T, EC, t) E, TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, EC, dbl)(TEMPLATE(T, t) x3, TEMPLATE(T, t) z3,
                      const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                      const TEMPLATE3(T, EC, t) E,
                      const TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, EC, dadd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5,
                       const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                       const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const TEMPLATE3(T, EC, t) E,
                       const TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, EC, ladd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5, TEMPLATE(T, t) x4, TEMPLATE(T, t) z4,
                       const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                       const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const TEMPLATE3(T, EC, t) E,
                       const TEMPLATE(T, ctx_t) K);
void
TEMPLATE3(T, EC, mul_ltr)(TEMPLATE(T, t) x, TEMPLATE(T, t) z,
                          const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                          const fmpz_t m,
                          const TEMPLATE3(T, EC, t) E,
                          const TEMPLATE(T, ctx_t) K);

#endif
#endif
