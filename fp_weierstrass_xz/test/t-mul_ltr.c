#include <stdio.h>
#include <time.h>
#include "fp_weierstrass_xz.h"

#define ITER 1000

int main(void)
{
    clock_t start;
    int i;
    fmpz_t p;
    fp_ctx_t K;
    fp_t a, b, x, z, xx, zz;
    fmpz_t m;
    fp_weierstrass_xz_t E;

    fmpz_init(p);
    fmpz_set_str(p, "4000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000f", 16);
    fp_ctx_init(K, p);

    fp_init(a, K);
    fp_init(b, K);
#ifndef USE_REDC
    fp_set_str(a, "9e601ddfcbe3574cc0670fd0348941ce6e902b051acd739d1785b968f1e2a8d6e3a4d3e7c9fe4ec54262122e5599d916b7356b160c838ce3081b9ae58d5bc3d", 16, K);
    fp_set_str(b, "1daada6902372a47141b1b9c63e05d15ba94c315c8661fc6c8130a9993cadec32f8940ea3ffaafd6f3dbad98f4a8a354481dad96846b54ed4884a5345da6286e", 16, K);
#else
    fp_set_str(a, "2e178ff8c36b78a02e7d84b33afd4939e1635eacdb7d8e72e7ca88b674ee06da2a55e55aca86589c47103bd23eff11ead0f7ae6d5112afaca1987b232de7e24a", 16, K);
    fp_set_str(b, "bf4cf637b12175749a58758976a2ee8452246e508108d691b8984015c73ca40dbd4c919013ec99ed8835026a879b83f190b50b8f6d8186300e947ba0d0e87dc", 16, K);
#endif
    fp_init(x, K);
    fp_set_str(x, "3705594e2087e6d6f157f683944a68a384e776e11c81ec24e9b8a47ac1b42e7dce41acb2cc3de1fe22400d9d386ea2069873662fdcd64acf1176d68dccf018f6", 16, K);
    fp_init(z, K);
    fp_set_str(z, "1", 16, K);
    fp_init(xx, K);
    fp_init(zz, K);

    fp_weierstrass_xz_init(E, K);
    fp_weierstrass_xz_set_fp(E, a, b, K);

    fmpz_init(m);
    fmpz_set_str(m, "3236dc8b8ad8959a9f5a83fb72bbb8b673d30893dfd94cece36a4c45960e890d7d7c6de58b25d420e89a6b9ab41844f658387bfe15aadf7896e1c58084a4c348", 16);

    start = clock();
    for (i = 0; i < ITER; i++)
        fp_weierstrass_xz_mul_ltr(xx, zz, x, z, m, E, K);
    printf("Took %f seconds.\n", (double) (clock() - start)/(ITER*CLOCKS_PER_SEC));

    printf("xx: ");
    fmpz_print(xx);
    printf("\n");
    printf("zz: ");
    fmpz_print(zz);
    printf("\n");

    fmpz_clear(m);
    fp_weierstrass_xz_clear(E, K);
    fp_clear(a, K);
    fp_clear(b, K);
    fp_clear(x, K);
    fp_clear(z, K);
    fp_clear(xx, K);
    fp_clear(zz, K);
    fp_ctx_clear(K);
    fmpz_clear(p);

    return 0;
}
