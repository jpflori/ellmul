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

#include "fp_weierstrass_xz.h"

#ifdef T
#undef T
#endif

#define T fp
#ifdef USE_REDC
#define mul CONCAT(mul_redc, USE_REDC)
#define sqr CONCAT(sqr_redc, USE_REDC)
#endif
#include "fq_weierstrass_xz_templates/dadd.c"
#ifdef USE_REDC
#undef mul
#undef sqr
#endif
#undef T
