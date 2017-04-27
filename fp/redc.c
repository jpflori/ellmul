#include "gmp.h"
/* From gmp-impl.h
 * void mpn_mul_basecase (mp_ptr, mp_srcptr, mp_size_t, mp_srcptr, mp_size_t);
 */
/* Not exposed by GMP but by MPIR (with different signatures...)
 * mp_limb_t/void mpn_redc_1 (mp_ptr, mp_ptr, mp_srcptr, mp_size_t, mp_limb_t);
 */
#include "fmpz.h"
#include "fp.h"

void fp_redc(fp_t rop, const fmpz_t op, const fp_ctx_t ctx)
{
    fmpz_t t, m;
    fmpz_init(t);
    fmpz_init(m);
    fmpz_fdiv_r_2exp(t, op, ctx->b_exp);
    fmpz_mul(m, t, ctx->mpinv);
    fmpz_fdiv_r_2exp(m, m, ctx->b_exp);
    fmpz_mul(t, m, ctx->p);
    fmpz_add(rop, t, op);
    fmpz_fdiv_q_2exp(rop, rop, ctx->b_exp);
    if (fmpz_cmp(rop, ctx->p) >= 0)
        fmpz_sub(rop, rop, ctx->p);
    fmpz_clear(t);
    fmpz_clear(m);
}

void fp_mul_redc(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    fmpz_mul(rop, op1, op2);
    fp_redc(rop, rop, ctx);
}

void fp_redcify(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fp_mul_redc(rop, op, ctx->b_square, ctx);
}

void fp_unredcify(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    fp_redc(rop, op, ctx);
}

void fp_redc_mpir(fp_t rop, const fmpz_t op, const fp_ctx_t ctx)
{
    ulong len_in, len_p, len_out;

    __mpz_struct *zrop, *zp, *zmpinv;
    mp_ptr in, p;
    mp_ptr out;
    mp_limb_t mpinv;

    zp = COEFF_TO_PTR(*(ctx->p));
    len_p = zp->_mp_size;
    zmpinv = COEFF_TO_PTR(*(ctx->mpinv));
    in = flint_malloc(2*len_p*sizeof(mp_limb_t));
    p = zp->_mp_d;
    mpinv = zmpinv->_mp_d[0];

    if (!COEFF_IS_MPZ(*op))
    {
        len_in = 1;
        in[0] = *op;
    } else {
        __mpz_struct *zop = COEFF_TO_PTR(*op);
        len_in = zop->_mp_size;
        mpn_copyi(in, zop->_mp_d, len_in);
    }
    mpn_zero(in + len_in, 2*len_p - len_in);

    /* We don't need input anymore so we don't care about aliasing */
    zrop = _fmpz_promote(rop);
    if (zrop->_mp_alloc < len_p)
        _mpz_realloc(zrop, len_p);
    out = zrop->_mp_d;
    mpn_redc_1(out, in, p, len_p, mpinv);
    if (mpn_cmp(out, p, len_p) >= 0)
        mpn_sub_n(out, out, p, len_p);
    len_out = len_p;
    while(!out[--len_out])
        ;
    zrop->_mp_size = ++len_out;
    if (len_out <= 1)
        _fmpz_demote_val(rop);

    flint_free(in);
}

void fp_mul_redc_mpir(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    ulong len_a, len_b, len_p, len_out;

    __mpz_struct *zop1, *zop2, *zrop, *zp, *zmpinv;
    mp_srcptr a, b, p;
    mp_ptr out, prod;
    mp_limb_t mpinv;

    zp = COEFF_TO_PTR(*(ctx->p));
    len_p = zp->_mp_size;
    zmpinv = COEFF_TO_PTR(*(ctx->mpinv));
    prod = flint_malloc(2*len_p*sizeof(mp_limb_t));
    p = zp->_mp_d;
    mpinv = zmpinv->_mp_d[0];

    if (!COEFF_IS_MPZ(*op1))
    {
        if (!COEFF_IS_MPZ(*op2))
        {
            len_a = 1;
            len_b = 1;
            umul_ppmm(prod[1], prod[0], *op1, *op2);
        } else {
            len_a = 1;
            zop2 = COEFF_TO_PTR(*op2);
            b = zop2->_mp_d;
            len_b = zop2->_mp_size;
            prod[len_b] = mpn_mul_1(prod, b, len_b, *op1);
        }
    } else {
        if (!COEFF_IS_MPZ(*op2))
        {
            zop1 = COEFF_TO_PTR(*op1);
            a = zop1->_mp_d;
            len_a = zop1->_mp_size;
            len_b = 1;
            prod[len_a] = mpn_mul_1(prod, a, len_a, *op2);
        } else {
            zop1 = COEFF_TO_PTR(*op1);
            a = zop1->_mp_d;
            len_a = zop1->_mp_size;
            zop2 = COEFF_TO_PTR(*op2);
            b = zop2->_mp_d;
            len_b = zop2->_mp_size;

            /* mpn_mul_basecase(prod, a, len_a, b, len_b); */
            if (len_a >= len_b)
                mpn_mul(prod, a, len_a, b, len_b);
            else
                mpn_mul(prod, b, len_b, a, len_a);
        }
    }
    mpn_zero(prod + len_a + len_b, 2*len_p - len_a - len_b);

    /* We don't need inputs anymore so we don't care about aliasing */
    zrop = _fmpz_promote(rop);
    if (zrop->_mp_alloc < len_p)
        _mpz_realloc(zrop, len_p);
    out = zrop->_mp_d;
    mpn_redc_1(out, prod, p, len_p, mpinv);
    if (mpn_cmp(out, p, len_p) >= 0)
        mpn_sub_n(out, out, p, len_p);
    len_out = len_p;
    while(!out[--len_out])
        ;
    zrop->_mp_size = ++len_out;
    if (len_out <= 1)
        _fmpz_demote_val(rop);

    flint_free(prod);
}

void fp_sqr_redc_mpir(fp_t rop, const fp_t op, const fp_ctx_t ctx)
{
    ulong len_a, len_p, len_out;

    __mpz_struct *zop, *zrop, *zp, *zmpinv;
    mp_srcptr a, p;
    mp_ptr out, prod;
    mp_limb_t mpinv;

    zp = COEFF_TO_PTR(*(ctx->p));
    len_p = zp->_mp_size;
    zmpinv = COEFF_TO_PTR(*(ctx->mpinv));
    prod = flint_malloc(2*len_p*sizeof(mp_limb_t));
    p = zp->_mp_d;
    mpinv = zmpinv->_mp_d[0];

    if (!COEFF_IS_MPZ(*op))
    {
        len_a = 1;
        umul_ppmm(prod[1], prod[0], *op, *op);
    } else {
        zop = COEFF_TO_PTR(*op);
        a = zop->_mp_d;
        len_a = zop->_mp_size;
        mpn_sqr(prod, a, len_a);
    }
    mpn_zero(prod + 2*len_a, 2*len_p - 2*len_a);

    /* We don't need inputs anymore so we don't care about aliasing */
    zrop = _fmpz_promote(rop);
    if (zrop->_mp_alloc < len_p)
        _mpz_realloc(zrop, len_p);
    out = zrop->_mp_d;
    mpn_redc_1(out, prod, p, len_p, mpinv);
    if (mpn_cmp(out, p, len_p) >= 0)
        mpn_sub_n(out, out, p, len_p);
    len_out = len_p;
    while(!out[--len_out])
        ;
    zrop->_mp_size = ++len_out;

    flint_free(prod);
}

void fp_mul_redc_cios(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    ulong i, len_a, len_b, len_p, len_out;

    __mpz_struct *zop1, *zop2, *zrop, *zp, *zmpinv;
    mp_srcptr a, b, p;
    mp_ptr out, prod;
    mp_limb_t mpinv, m, carry;

    zp = COEFF_TO_PTR(*(ctx->p));
    len_p = zp->_mp_size;
    zmpinv = COEFF_TO_PTR(*(ctx->mpinv));
    prod = flint_malloc(2*len_p*sizeof(mp_limb_t) + 1);
    p = zp->_mp_d;
    mpinv = zmpinv->_mp_d[0];

    /* i = 0; */
    if (!COEFF_IS_MPZ(*op1))
    {
        if (!COEFF_IS_MPZ(*op2))
        {
            len_a = 1;
            len_b = 1;
            umul_ppmm(prod[1], prod[0], *op1, *op2);
        } else {
            len_a = 1;
            zop2 = COEFF_TO_PTR(*op2);
            b = zop2->_mp_d;
            len_b = zop2->_mp_size;
            prod[len_b] = mpn_mul_1(prod, b, len_b, *op1);
        }
    } else {
        if (!COEFF_IS_MPZ(*op2))
        {
            len_a = 1;
            zop1 = COEFF_TO_PTR(*op1);
            b = zop1->_mp_d;
            len_b = zop1->_mp_size;
            prod[len_b] = mpn_mul_1(prod, b, len_b, *op2);
        } else {
            zop1 = COEFF_TO_PTR(*op1);
            a = zop1->_mp_d;
            len_a = zop1->_mp_size;
            zop2 = COEFF_TO_PTR(*op2);
            b = zop2->_mp_d;
            len_b = zop2->_mp_size;
            prod[len_b] = mpn_mul_1(prod, b, len_b, a[0]);
        }
    }
    {
        if (len_b < len_p)
            mpn_zero(prod + len_b + 1, len_p - len_b);

        /* REDC */
        m = prod[0] * mpinv;
        carry = mpn_addmul_1(prod, p, len_p, m);
        prod[len_p + 1] = mpn_add_1(prod + len_p, prod + len_p, 1, carry);

        /* Shift */
        prod++;
    }

    i = 1;
    for (; i < len_a; i++) {
        /* Product */
        carry = mpn_addmul_1(prod, b, len_b, a[i]);
        prod[len_p + 1] = mpn_add_1(prod + len_b, prod + len_b, len_p + 1 - len_b, carry);

        /* REDC */
        m = prod[0] * mpinv;
        carry = mpn_addmul_1(prod, p, len_p, m);
        prod[len_p + 1] += mpn_add_1(prod + len_p, prod + len_p, 1, carry);

        /* Shift */
        prod++;
    }

    /* i = len_a; */
    for (; i < len_p; i++) {
        /* REDC */
        m = prod[0] * mpinv;
        carry = mpn_addmul_1(prod, p, len_p, m);
        prod[len_p + 1] = mpn_add_1(prod + len_p, prod + len_p, 1, carry);

        /* Shift */
        prod++;
    }

    /* We don't need inputs anymore so we don't care about aliasing */
    zrop = _fmpz_promote(rop);
    if (zrop->_mp_alloc < len_p)
        _mpz_realloc(zrop, len_p);
    out = zrop->_mp_d;
    if (prod[len_p] || mpn_cmp(prod, p, len_p) >= 0)
    {
        mpn_sub_n(out, prod, p, len_p);
    } else {
        mpn_copyi(out, prod, len_p);
    }
    len_out = len_p;
    while(!out[--len_out])
        ;
    zrop->_mp_size = ++len_out;
    if (len_out <= 1)
        _fmpz_demote_val(rop);

    flint_free(prod-len_p);
}
