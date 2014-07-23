#include "flint/fq.h"

#define T fq

#include "ellmul.c"

// Test routine
int main(int argc, char* argv[])
{
    unsigned long n;
    n = 300;

    ulong plong = 1031;
    fmpz_t p;
    fmpz_init(p);
    fmpz_set_ui(p, plong);

    fmpz_t q;
    fmpz_init(q);
    fmpz_pow_ui(q, p, n);
    //fmpz_set_ui(q, 4);
    fmpz_print(q);
    printf("\n");

    fmpz_mod_poly_t f;
    fmpz_mod_poly_init(f, p);
    fmpz_mod_poly_set_coeff_ui(f, 300, 1);
    fmpz_mod_poly_set_coeff_ui(f, 1, 1);
    fmpz_mod_poly_set_coeff_ui(f, 0, 321);

    TEMPLATE(T, ctx_t) K;
    TEMPLATE(T, ctx_init_modulus)(K, f, "a");

    TEMPLATE(T, t) inv4;
    TEMPLATE(T, init)(inv4, K);
    TEMPLATE(T, set_ui)(inv4, 4, K);
    TEMPLATE(T, inv)(inv4, inv4, K);

    TEMPLATE(T, t) x1, z1;
    TEMPLATE(T, init)(x1, K);
    TEMPLATE(T, init)(z1, K);
    TEMPLATE(T, one)(z1, K);
    TEMPLATE(T, gen)(x1, K);
    TEMPLATE(T, add)(x1, x1, z1, K);

    TEMPLATE(T, t) x, z;
    TEMPLATE(T, init)(x, K);
    TEMPLATE(T, init)(z, K);

    clock_t t;

    // Weierstrass XZ
    TEMPLATE(T, weierstrass_xz_t) wxz;
    TEMPLATE(T, weierstrass_xz_init)(wxz, K);
    TEMPLATE(T, weierstrass_xz_set_ui)(wxz, 207, 693, K);

    t = clock();
    TEMPLATE(T, weierstrass_xz_mul_ltr)(x, z, x1, z1, q, wxz, K);
    printf("Weierstrass XZ: %lf\n", ((double) clock() - t) / CLOCKS_PER_SEC);
    TEMPLATE(T, inv)(z1, z, K);
    TEMPLATE(T, mul)(x1, x, z1, K);
    TEMPLATE(T, print_pretty)(x1, K);
    printf("\n");

    TEMPLATE(T, weierstrass_xz_clear)(wxz, K);

    // Montgomery XZ
    TEMPLATE(T, montgomery_xz_t) mxz;
    TEMPLATE(T, montgomery_xz_init)(mxz, K);
    TEMPLATE(T, montgomery_xz_set_ui)(mxz, 207, 693, K, inv4);

    t = clock();
    TEMPLATE(T, montgomery_xz_mul_ltr)(x, z, x1, z1, q, mxz, K);
    printf("Montgomery XZ: %lf\n", ((double) clock() - t) / CLOCKS_PER_SEC);
    TEMPLATE(T, inv)(z1, z, K);
    TEMPLATE(T, mul)(x1, x, z1, K);
    TEMPLATE(T, print_pretty)(x1, K);
    printf("\n");

    TEMPLATE(T, montgomery_xz_clear)(mxz, K);

    TEMPLATE(T, clear)(x, K);
    TEMPLATE(T, clear)(z, K);
    TEMPLATE(T, clear)(x1, K);
    TEMPLATE(T, clear)(z1, K);
    TEMPLATE(T, ctx_clear)(K);
    fmpz_mod_poly_clear(f);
    fmpz_clear(q);
    fmpz_clear(p);

    return 0;
}
