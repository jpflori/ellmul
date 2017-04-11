#include <stdio.h>
#include <time.h>
#include "fp_weierstrass_xz.h"

#define ITER 100

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
    fp_set_str(a, "170bc7fc61b5bc50173ec2599d7ea49cf0b1af566dbec73973e5445b3a77036d152af2ad65432c4e23881de91f7f88f5687bd736a88957d650cc3d9196f3f125", 16, K);
    fp_set_str(b, "5fa67b1bd890baba4d2c3ac4bb5177422912372840846b48dc4c200ae39e5206dea648c809f64cf6c41a813543cdc1f8c85a85c7b6c0c318074a3dd068743ee", 16, K);
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
