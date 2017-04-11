#include "fmpz.h"
#include "fp.h"

void fp_mul(fp_t rop, const fp_t op1, const fp_t op2, const fp_ctx_t ctx)
{
    fmpz_t q;
    fmpz_init(q);

    fmpz_mul(rop, op1, op2);
    fmpz_fdiv_qr_preinvn(q, rop, rop, ((ctx)->p), ((ctx)->pinv));

    fmpz_clear(q);
}

