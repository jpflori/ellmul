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

/*
  Montgomery XZ type
*/
typedef struct {
    TEMPLATE(T, t) a;
    TEMPLATE(T, t) a24;
    TEMPLATE(T, t) b;
} TEMPLATE(T, montgomery_xz_struct);

typedef TEMPLATE(T, montgomery_xz_struct) TEMPLATE(T, montgomery_xz_t)[1];

void
TEMPLATE(T, montgomery_xz_set_ui)(TEMPLATE(T, montgomery_xz_t) E,
                                  ulong a, ulong b,
                                  const TEMPLATE(T, t) inv4,
                                  const TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, montgomery_xz_set, T)(TEMPLATE(T, montgomery_xz_t) E,
                                   const TEMPLATE(T, t) a, const TEMPLATE(T, t) b,
                                   const TEMPLATE(T, t) inv4,
                                   const TEMPLATE(T, ctx_t) K);

#ifdef EC
#undef EC
#endif

#define EC montgomery_xz
#include "fq_elliptic_curves_templates.h"
#undef EC

#endif
