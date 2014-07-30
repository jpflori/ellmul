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

void
TEMPLATE3(T, weierstrass_xz_set, T)(TEMPLATE(T, weierstrass_xz_t) E,
                                   const TEMPLATE(T, t) a, const TEMPLATE(T, t) b,
                                   const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, set)(E->a, a, K);
    TEMPLATE(T, set)(E->b, b, K);
    TEMPLATE(T, mul_ui)(E->b2, E->b, 2, K);
    TEMPLATE(T, mul_ui)(E->b4, E->b, 4, K);
}

#endif
