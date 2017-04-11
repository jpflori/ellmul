#include "fmpz.h"
#include "fp.h"

void fp_ctx_init(fp_ctx_t ctx, const fmpz_t p)
{
    fmpz_t q, b;

    fmpz_init(ctx->p);
    fmpz_set(ctx->p, p);

    ctx->sparse_modulus = 0;
    ctx->j = NULL;
    ctx->len = NULL;

    fmpz_preinvn_init(ctx->pinv, p);

    ctx->b_exp = fmpz_bits(p);

    fmpz_init(ctx->b_square);
    fmpz_one(ctx->b_square);
    fmpz_mul_2exp(ctx->b_square, ctx->b_square, 2*ctx->b_exp);
    fmpz_init(q);
    fmpz_fdiv_qr_preinvn(q, ctx->b_square, ctx->b_square, p, ctx->pinv);
    fmpz_clear(q);

    fmpz_init(ctx->mpinv);
    fmpz_init(b);
    fmpz_one(b);
    fmpz_mul_2exp(b, b, ctx->b_exp);
    fmpz_invmod(ctx->mpinv, p, b);
    fmpz_sub(ctx->mpinv, b, ctx->mpinv);
}

