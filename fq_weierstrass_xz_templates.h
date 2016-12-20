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
  Weierstrass XZ type
*/
typedef struct {
    TEMPLATE(T, t) a;
    TEMPLATE(T, t) b;
    TEMPLATE(T, t) b2;
    TEMPLATE(T, t) b4;
} TEMPLATE(T, weierstrass_xz_struct);

typedef TEMPLATE(T, weierstrass_xz_struct) TEMPLATE(T, weierstrass_xz_t)[1];

void
TEMPLATE(T, weierstrass_xz_set_ui)(TEMPLATE(T, weierstrass_xz_t) E,
                                   ulong a, ulong b,
                                   const TEMPLATE(T, ctx_t) K);

void
TEMPLATE3(T, weierstrass_xz_set, T)(TEMPLATE(T, weierstrass_xz_t) E,
                                    const TEMPLATE(T, t) a, const TEMPLATE(T, t) b,
                                    const TEMPLATE(T, ctx_t) K);

#ifdef EC
#undef EC
#endif

#define EC weierstrass_xz
#include "fq_elliptic_curves_templates.h"
#undef EC

#endif
