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

    Copyright (C) 2017 Jean-Pierre Flori

******************************************************************************/

#ifndef FP_H
#define FP_H

#include "fmpz.h"

typedef fmpz_t fp_t;
typedef fmpz fp;

typedef struct
{
    fmpz_t p;

    /* For sparse p; currently unused */
    int sparse_modulus; 
    slong *j;
    slong *len;

    /* For faster division */
    fmpz_preinvn_t pinv;

    /* For REDC */
    mp_bitcnt_t b_exp;
    fmpz_t b_square;
    fmpz_t mpinv;
}
fp_ctx_struct;

typedef fp_ctx_struct fp_ctx_t[1];

void fp_ctx_init(fp_ctx_t ctx, const fmpz_t p);
void fp_ctx_clear(fp_ctx_t ctx);

static __inline__
void fp_init(fp_t rop, const fp_ctx_t p)
{
    fmpz_init(rop);
}

static __inline__
void fp_clear(fp_t rop, const fp_ctx_t p)
{
    fmpz_clear(rop);
}

static __inline__
void fp_set(fp_t rop, const fp_t op, const fp_ctx_t p)
{
    fmpz_set(rop, op);
}

static __inline__
void fp_set_fmpz(fp_t rop, const fmpz_t op, const fp_ctx_t p)
{
    fmpz_set(rop, op); /* FIXME */
}

static __inline__
void fp_set_ui(fp_t rop, slong op, const fp_ctx_t p)
{
    fmpz_set_ui(rop, op); /* FIXME */
}

static __inline__
void fp_set_str(fp_t rop, char * str, slong base, const fp_ctx_t p)
{
    fmpz_set_str(rop, str, base); /* FIXME */
}

static __inline__
void fp_zero(fp_t rop, const fp_ctx_t p)
{
    fmpz_zero(rop);
}

static __inline__
void fp_one(fp_t rop, const fp_ctx_t p)
{
    fmpz_one(rop);
}

static __inline__
void fp_swap(fp_t rop, fp_t op, const fp_ctx_t p)
{
    fmpz_swap(rop, op);
}

static __inline__
void fp_add(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    fmpz_add(rop, op1, op2);
    if (fmpz_cmp(rop, ctx->p) >= 0)
        fmpz_sub(rop, rop, ((ctx)->p));
}

static __inline__
void fp_sub(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    fmpz_sub(rop, op1, op2);
    if (fmpz_sgn(rop) < 0)
        fmpz_add(rop, rop, ((ctx)->p));
}

void fp_mul(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx);
void fp_mul_ui(fp_t rop, const fp_t op1, ulong op2, const fp_ctx_t ctx);

static __inline__
void fp_sqr(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fp_mul(rop, op, op, ctx);
}

void fp_redc(fp_t rop, const fp_t op, const fp_ctx_t ctx);
void fp_mul_redc(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx);
static __inline__
void fp_sqr_redc(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fp_mul_redc(rop, op, op, ctx);
}
void fp_redcify(fp_t rop, const fp_t op, const fp_ctx_t ctx);
void fp_unredcify(fp_t rop, const fp_t op, const fp_ctx_t ctx);

void fp_mul_redc_mpir(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx);
void fp_sqr_redc_mpir(fp_t rop, const fp_t op, const fp_ctx_t ctx);

void fp_mul_redc_cios(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx);
static __inline__
void fp_sqr_redc_cios(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fp_mul_redc_cios(rop, op, op, ctx);
}

static __inline__
void fp_inv(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fmpz_invmod(rop, op, ((ctx)->p));
}

#endif
