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

/* requires the inverse of 4 in K */
void
TEMPLATE3(T, montgomery_xz_set, T)(TEMPLATE(T, montgomery_xz_t) E,
                                   const TEMPLATE(T, t) a, const TEMPLATE(T, t) b,
                                   const TEMPLATE(T, t) inv4,
                                   const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, set)(E->a, a, K);
    TEMPLATE(T, set_ui)(E->a24, 2, K);
    TEMPLATE(T, add)(E->a24, E->a, E->a24, K);
    TEMPLATE(T, mul)(E->a24, E->a24, inv4, K);
    TEMPLATE(T, set)(E->b, b, K);
}

#endif
