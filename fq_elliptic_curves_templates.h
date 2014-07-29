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
#ifdef EC

void
TEMPLATE3(T, EC, mul_ltr)(TEMPLATE(T, t) x, TEMPLATE(T, t) z,
                          const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                          const fmpz_t m,
                          const TEMPLATE3(T, EC, t) E,
                          const TEMPLATE(T, ctx_t) K);

#endif
#endif
