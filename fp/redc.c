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
