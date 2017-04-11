#include "fmpz.h"
#include "fp.h"

void fp_ctx_clear(fp_ctx_t ctx)
{
    fmpz_clear(((ctx)->p));

    fmpz_preinvn_clear(((ctx)->pinv));

    fmpz_clear(((ctx)->b_square));
    fmpz_clear(((ctx)->mpinv));
}
