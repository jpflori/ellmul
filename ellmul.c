// copyright blablabla GPL v3 or later at your option

#ifdef T

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpz_mod_poly.h"
#include "flint/fq.h"

// support nails?
// make the curve types ctx like? and add a point type?
// store finit field ctx within curve ctx?
// make the multiplication routine generic? using function pointers?
// pass tmp space to functions?
// use windows?
// return z==0 for bdl, dadd, ladd, mul functions?
// add z=1 formulas?
// add conversion function, extract y from x?

// Montgomery XZ type
typedef struct {
    TEMPLATE(T, t) a;
    TEMPLATE(T, t) a24;
    TEMPLATE(T, t) b;
//    TEMPLATE(T, ctx_t) K;
} TEMPLATE(T, montgomery_xz_struct);

typedef TEMPLATE(T, montgomery_xz_struct) TEMPLATE(T, montgomery_xz_t)[1];

void __inline__
TEMPLATE(T, montgomery_xz_init)(TEMPLATE(T, montgomery_xz_t) E, TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, init)(E->a, K);
    TEMPLATE(T, init)(E->a24, K);
    TEMPLATE(T, init)(E->b, K);
}

// requires the inverse of 4 in K
void TEMPLATE(T, montgomery_xz_set_ui)(TEMPLATE(T, montgomery_xz_t) E, ulong a, ulong b, TEMPLATE(T, ctx_t) K, TEMPLATE(T, t) inv4)
{
    TEMPLATE(T, set_ui)(E->a, a, K);
    TEMPLATE(T, set_ui)(E->a24, 2, K);
    TEMPLATE(T, add)(E->a24, E->a, E->a24, K);
    TEMPLATE(T, mul)(E->a24, E->a24, inv4, K);
    TEMPLATE(T, set_ui)(E->b, b, K);
}

void TEMPLATE(T, montgomery_xz_clear)(TEMPLATE(T, montgomery_xz_t) E, TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, clear)(E->a, K);
    TEMPLATE(T, clear)(E->a24, K);
    TEMPLATE(T, clear)(E->b, K);
}

// Montgomery XZ doubling
// P1 |-> P3=[2]P1
// dbl-1987-m-3
void TEMPLATE(T, montgomery_xz_dbl)(TEMPLATE(T, t) x3, TEMPLATE(T, t) z3,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const TEMPLATE(T, montgomery_xz_t) E,
                       const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, aa, bb;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(aa, K);
    TEMPLATE(T, init)(bb, K);

    // A = X1+Z1
    TEMPLATE(T, add)(a, x1, z1, K);
    // AA = A²
    TEMPLATE(T, sqr)(aa, a, K);
    // B = X1-Z1
    TEMPLATE(T, sub)(b, x1, z1, K);
    // BB = B²
    TEMPLATE(T, sqr)(bb, b, K);
    // C = AA-BB
    TEMPLATE(T, sub)(x3, aa, bb, K);
    // Z3 = C*(BB+a24*C)
    TEMPLATE(T, mul)(a, E->a24, x3, K);
    TEMPLATE(T, add)(b, bb, a, K);
    TEMPLATE(T, mul)(z3, x3, b, K);
    // X3 = AA*BB
    TEMPLATE(T, mul)(x3, aa, bb, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(aa, K);
    TEMPLATE(T, clear)(bb, K);
}

// Montgomery XZ differential addition
// P1=P3-P2, P2, P3 |-> P5=P2+P3
// dadd-1987-m-3
void TEMPLATE(T, montgomery_xz_dadd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5,
                        const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                        const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                        const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                        const TEMPLATE(T, montgomery_xz_t) E,
                        const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, c, d, da, cb;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);
    TEMPLATE(T, init)(da, K);
    TEMPLATE(T, init)(cb, K);

    // A = X2+Z2
    TEMPLATE(T, add)(a, x2, z2, K);
    // B = X2-Z2
    TEMPLATE(T, sub)(b, x2, z2, K);
    // C = X3+Z3
    TEMPLATE(T, add)(c, x3, z3, K);
    // D = X3-Z3
    TEMPLATE(T, sub)(d, x3, z3, K);
    // DA = D*A
    TEMPLATE(T, mul)(da, d, a, K);
    // CB = C*B
    TEMPLATE(T, mul)(cb, c, b, K);
    // X5 = Z1*(DA+CB)²
    TEMPLATE(T, add)(x5, da, cb, K);
    TEMPLATE(T, sqr)(a, x5, K);
    TEMPLATE(T, mul)(x5, z1, a, K);
    // Z5 = X1*(DA-CB)²
    TEMPLATE(T, sub)(z5, da, cb, K);
    TEMPLATE(T, sqr)(a, z5, K);
    TEMPLATE(T, mul)(z5, x1, a, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
    TEMPLATE(T, clear)(da, K);
    TEMPLATE(T, clear)(cb, K);
}

// Montgomery XZ ladder
// P1=P3-P2, P2, P3 |-> P4=[2]P2, P5=P2+P3
// ladd-1987-m-3
void TEMPLATE(T, montgomery_xz_ladd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5, TEMPLATE(T, t) x4, TEMPLATE(T, t) z4,
                        const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                        const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                        const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                        const TEMPLATE(T, montgomery_xz_t) E,
                        const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, c, d, e, aa, bb;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);
    TEMPLATE(T, init)(e, K);
    TEMPLATE(T, init)(aa, K);
    TEMPLATE(T, init)(bb, K);

    // A = X2+Z2
    TEMPLATE(T, add)(a, x2, z2, K);
    // AA = A²
    TEMPLATE(T, sqr)(aa, a, K);
    // B = X2-Z2
    TEMPLATE(T, sub)(b, x2, z2, K);
    // BB = B²
    TEMPLATE(T, sqr)(bb, b, K);
    // E = AA-BB
    TEMPLATE(T, sub)(e, aa, bb, K);
    // X4 = AA*BB
    TEMPLATE(T, mul)(x4, aa, bb, K);
    // Z4 = E*(BB+a24*E)
    TEMPLATE(T, mul)(z4, E->a24, e, K);
    TEMPLATE(T, add)(aa, bb, z4, K);
    TEMPLATE(T, mul)(z4, e, aa, K);
    // C = X3+Z3
    TEMPLATE(T, add)(c, x3, z3, K);
    // D = X3-Z3
    TEMPLATE(T, sub)(d, x3, z3, K);
    // DA = D*A
    TEMPLATE(T, mul)(aa, d, a, K);
    // CB =  C*B
    TEMPLATE(T, mul)(bb, c, b, K);
    // X5 = Z1*(DA+CB)²
    TEMPLATE(T, add)(x5, aa, bb, K);
    TEMPLATE(T, sqr)(e, x5, K);
    TEMPLATE(T, mul)(x5, z1, e, K);
    // Z5 = X1*(DA-CB)²
    TEMPLATE(T, sub)(z5, aa, bb, K);
    TEMPLATE(T, sqr)(e, z5, K);
    TEMPLATE(T, mul)(z5, x1, e, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
    TEMPLATE(T, clear)(e, K);
    TEMPLATE(T, clear)(aa, K);
    TEMPLATE(T, clear)(bb, K);
}

// Montgomery XZ point multiplication
// Left-to-right binary exponentiation
TEMPLATE(T, montgomery_xz_mul_ltr)(TEMPLATE(T, t) x, TEMPLATE(T, t) z,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const fmpz_t m,
                       const TEMPLATE(T, montgomery_xz_t) E,
                       const TEMPLATE(T, ctx_t) K)
{
    if (fmpz_is_zero(m))
    {
        TEMPLATE(T, one)(x, K);
        TEMPLATE(T, zero)(z, K);
    }
    else if (fmpz_is_one(m))
    {
        TEMPLATE(T, set)(x, x1, K);
        TEMPLATE(T, set)(z, z1, K);
    }
    else
    {
        ulong bit;
        TEMPLATE(T, t) x2, z2, x3, z3, x4, z4, x5, z5;
        TEMPLATE(T, init)(x2, K);
        TEMPLATE(T, init)(z2, K);
        TEMPLATE(T, init)(x3, K);
        TEMPLATE(T, init)(z3, K);
        TEMPLATE(T, init)(x4, K);
        TEMPLATE(T, init)(z4, K);
        TEMPLATE(T, init)(x5, K);
        TEMPLATE(T, init)(z5, K);

        // P
        TEMPLATE(T, set)(x2, x1, K);
        TEMPLATE(T, set)(z2, z1, K);
        // 2*P
        TEMPLATE(T, montgomery_xz_dbl)(x3, z3, x1, z1, E, K);

        for (bit = fmpz_bits(m) - 2; bit > 0; bit--)
        {
            if (fmpz_tstbit(m, bit))
            {
                // P, [k+1]P, [k]P -> [2k+2]P, [2k+1]P
                TEMPLATE(T, montgomery_xz_ladd)(x5, z5, x4, z4, x2, z2, x3, z3, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x5, K);
                TEMPLATE(T, swap)(z2, z5, K);
                TEMPLATE(T, swap)(x3, x4, K);
                TEMPLATE(T, swap)(z3, z4, K);
            }
            else
            {
                // P, [k]P, [k+1]P -> [2k]P, [2k+1]P
                TEMPLATE(T, montgomery_xz_ladd)(x5, z5, x4, z4, x3, z3, x2, z2, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x4, K);
                TEMPLATE(T, swap)(z2, z4, K);
                TEMPLATE(T, swap)(x3, x5, K);
                TEMPLATE(T, swap)(z3, z5, K);
            }
        }

        // Last iteration, no need for a ladder
        {
            if (fmpz_tstbit(m, 0))
            {
                TEMPLATE(T, montgomery_xz_dadd)(x, z, x3, z3, x2, z2, x1, z1, E, K);
            }
            else
            {
                TEMPLATE(T, montgomery_xz_dbl)(x, z, x2, z2, E, K);
            }
        }

        TEMPLATE(T, clear)(x2, K);
        TEMPLATE(T, clear)(z2, K);
        TEMPLATE(T, clear)(x3, K);
        TEMPLATE(T, clear)(z3, K);
        TEMPLATE(T, clear)(x4, K);
        TEMPLATE(T, clear)(z4, K);
        TEMPLATE(T, clear)(x5, K);
        TEMPLATE(T, clear)(z5, K);
    }
}

// Weierstrass XZ type
typedef struct {
    TEMPLATE(T, t) a;
    TEMPLATE(T, t) b;
    TEMPLATE(T, t) b2;
    TEMPLATE(T, t) b4;
//    TEMPLATE(T, ctx_t) K;
} TEMPLATE(T, weierstrass_xz_struct);

typedef TEMPLATE(T, weierstrass_xz_struct) TEMPLATE(T, weierstrass_xz_t)[1];

void TEMPLATE(T, weierstrass_xz_init)(TEMPLATE(T, weierstrass_xz_t) E, TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, init)(E->a, K);
    TEMPLATE(T, init)(E->b, K);
    TEMPLATE(T, init)(E->b2, K);
    TEMPLATE(T, init)(E->b4, K);
}

void TEMPLATE(T, weierstrass_xz_set_ui)(TEMPLATE(T, weierstrass_xz_t) E, ulong a, ulong b, TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, set_ui)(E->a, a, K);
    TEMPLATE(T, set_ui)(E->b, b, K);
    TEMPLATE(T, mul_ui)(E->b2, E->b, 2, K);
    TEMPLATE(T, mul_ui)(E->b4, E->b, 4, K);
}

void TEMPLATE(T, weierstrass_xz_clear)(TEMPLATE(T, weierstrass_xz_t) E, TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, clear)(E->a, K);
    TEMPLATE(T, clear)(E->b, K);
    TEMPLATE(T, clear)(E->b2, K);
    TEMPLATE(T, clear)(E->b4, K);
}

// Weierstrass XZ doubling
// P1 |-> P3=[2]P1
// dbl-2002-bj-3
void TEMPLATE(T, weierstrass_xz_dbl)(TEMPLATE(T, t) x3, TEMPLATE(T, t) z3,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const TEMPLATE(T, weierstrass_xz_t) E,
                       const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) xx, zz, a, azz;
    TEMPLATE(T, init)(xx, K);
    TEMPLATE(T, init)(zz, K);
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(azz, K);

    // XX = X1²
    TEMPLATE(T, sqr)(xx, x1, K);
    // ZZ = Z1²
    TEMPLATE(T, sqr)(zz, z1, K);
    // A = 2*((X1+Z1)²-XX-ZZ)
    TEMPLATE(T, add)(azz, x1, z1, K);
    TEMPLATE(T, sqr)(a, azz, K);
    TEMPLATE(T, sub)(a, a, xx, K);
    TEMPLATE(T, sub)(a, a, zz, K);
    TEMPLATE(T, mul_ui)(a, a, 2, K);
    // aZZ = a*ZZ
    TEMPLATE(T, mul)(azz, E->a, zz, K);
    // X3 = (XX-aZZ)²-b2*A*ZZ
    TEMPLATE(T, sub)(x3, xx, azz, K);
    TEMPLATE(T, sqr)(x3, x3, K);
    TEMPLATE(T, mul)(z3, a, zz, K);
    TEMPLATE(T, mul)(z3, E->b2, z3, K);
    TEMPLATE(T, sub)(x3, x3, z3, K);
    // Z3 = A*(XX+aZZ)+b4*ZZ²
    TEMPLATE(T, add)(xx, xx, azz, K);
    TEMPLATE(T, mul)(azz, a, xx, K);
    TEMPLATE(T, sqr)(a, zz, K);
    TEMPLATE(T, mul)(zz, E->b4, a, K);
    TEMPLATE(T, add)(z3, azz, zz, K);

    TEMPLATE(T, clear)(xx, K);
    TEMPLATE(T, clear)(zz, K);
    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(azz, K);
}

// Weierstrass XZ differential addition
// P1=P3-P2, P2, P3 |-> P5=P2+P3
// dadd-2002-it-3
void TEMPLATE(T, weierstrass_xz_dadd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5,
                         const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                         const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                         const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                         const TEMPLATE(T, weierstrass_xz_t) E,
                         const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) a, b, c, d, e, f;
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);
    TEMPLATE(T, init)(e, K);
    TEMPLATE(T, init)(f, K);

    // A = X2*X3
    TEMPLATE(T, mul)(a, x2, x3, K);
    // B = Z2*Z3
    TEMPLATE(T, mul)(b, z2, z3, K);
    // C = X2*Z3
    TEMPLATE(T, mul)(c, x2, z3, K);
    // D = X3*Z2
    TEMPLATE(T, mul)(d, x3, z2, K);
    // E = (A-a*B)²
    TEMPLATE(T, mul)(e, E->a, b, K);
    TEMPLATE(T, sub)(f, a, e, K);
    TEMPLATE(T, sqr)(e, f, K);
    // F = 4*b*B*(C+D)
    TEMPLATE(T, add)(f, c, d, K);
    TEMPLATE(T, mul)(a, b, f, K);
    TEMPLATE(T, mul)(b, E->b, a, K);
    TEMPLATE(T, mul_ui)(f, b, 4, K);
    // X5 = Z1*(E-F)
    TEMPLATE(T, sub)(a, e, f, K);
    TEMPLATE(T, mul)(x5, z1, a, K);
    // Z5 = X1*(C-D)²
    TEMPLATE(T, sub)(a, c, d, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(z5, x1, b, K);

    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
    TEMPLATE(T, clear)(e, K);
    TEMPLATE(T, clear)(f, K);
}

// Weierstrass XZ ladder
// P1=P3-P2, P2, P3 |-> P4=[2]P2, P5=P2+P3
// ladd-2002-it-3
void TEMPLATE(T, weierstrass_xz_ladd)(TEMPLATE(T, t) x5, TEMPLATE(T, t) z5, TEMPLATE(T, t) x4, TEMPLATE(T, t) z4,
                         const TEMPLATE(T, t) x3, const TEMPLATE(T, t) z3,
                         const TEMPLATE(T, t) x2, const TEMPLATE(T, t) z2,
                         const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                         const TEMPLATE(T, weierstrass_xz_t) E,
                         const TEMPLATE(T, ctx_t) K)
{
    TEMPLATE(T, t) xx, zz, e, azz, a, b, c, d;
    TEMPLATE(T, init)(xx, K);
    TEMPLATE(T, init)(zz, K);
    TEMPLATE(T, init)(azz, K);
    TEMPLATE(T, init)(e, K);
    TEMPLATE(T, init)(a, K);
    TEMPLATE(T, init)(b, K);
    TEMPLATE(T, init)(c, K);
    TEMPLATE(T, init)(d, K);

    // XX = X2²
    TEMPLATE(T, sqr)(xx, x2, K);
    // ZZ = Z2²
    TEMPLATE(T, sqr)(zz, z2, K);
    // aZZ = a*ZZ
    TEMPLATE(T, mul)(azz, E->a, zz, K);
    // E = (X2+Z2)²-XX-ZZ
    TEMPLATE(T, add)(a, x2, z2, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, add)(a, xx, zz, K);
    TEMPLATE(T, sub)(e, b, a, K);
    // X4 = (XX-aZZ)²-b4*E*ZZ
    TEMPLATE(T, sub)(a, xx, azz, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(a, e, zz, K);
    TEMPLATE(T, mul)(c, a, E->b4, K);
    TEMPLATE(T, sub)(x4, b, c, K);
    // Z4 = 2*E*(XX+aZZ)+b4*ZZ²
    TEMPLATE(T, add)(a, xx, azz, K);
    TEMPLATE(T, mul)(b, a, e, K);
    TEMPLATE(T, mul_ui)(c, b, 2, K);
    TEMPLATE(T, sqr)(a, zz, K);
    TEMPLATE(T, mul)(b, E->b4, a, K);
    TEMPLATE(T, add)(z4, c, b, K);
    // A = X2*X3
    TEMPLATE(T, mul)(a, x2, x3, K);
    // B = Z2*Z3
    TEMPLATE(T, mul)(b, z2, z3, K);
    // C = X2*Z3
    TEMPLATE(T, mul)(c, x2, z3, K);
    // D = X3*Z2
    TEMPLATE(T, mul)(d, x3, z2, K);
    // X5 = Z1*((A-a*B)²-b4*B*(C+D))
    TEMPLATE(T, mul)(azz, E->a, b, K);
    TEMPLATE(T, sub)(x5, a, azz, K);
    TEMPLATE(T, sqr)(a, x5, K);
    TEMPLATE(T, add)(x5, c, d, K);
    TEMPLATE(T, mul)(e, b, x5, K);
    TEMPLATE(T, mul)(b, E->b4, e, K);
    TEMPLATE(T, sub)(e, a, b, K);
    TEMPLATE(T, mul)(x5, z1, e, K);
    // Z5 = X1*(C-D)²
    TEMPLATE(T, sub)(a, c, d, K);
    TEMPLATE(T, sqr)(b, a, K);
    TEMPLATE(T, mul)(z5, x1, b, K);

    TEMPLATE(T, clear)(xx, K);
    TEMPLATE(T, clear)(zz, K);
    TEMPLATE(T, clear)(azz, K);
    TEMPLATE(T, clear)(e, K);
    TEMPLATE(T, clear)(a, K);
    TEMPLATE(T, clear)(b, K);
    TEMPLATE(T, clear)(c, K);
    TEMPLATE(T, clear)(d, K);
}

// Make the following generic? using hackish pointers?
// Weierstrass XZ point multiplication
// Left-to-right binary exponentiation
TEMPLATE(T, weierstrass_xz_mul_ltr)(TEMPLATE(T, t) x, TEMPLATE(T, t) z,
                       const TEMPLATE(T, t) x1, const TEMPLATE(T, t) z1,
                       const fmpz_t m,
                       const TEMPLATE(T, weierstrass_xz_t) E,
                       const TEMPLATE(T, ctx_t) K)
{
    if (fmpz_is_zero(m))
    {
        TEMPLATE(T, one)(x, K);
        TEMPLATE(T, zero)(z, K);
    }
    else if (fmpz_is_one(m))
    {
        TEMPLATE(T, set)(x, x1, K);
        TEMPLATE(T, set)(z, z1, K);
    }
    else
    {
        ulong bit;
        TEMPLATE(T, t) x2, z2, x3, z3, x4, z4, x5, z5;
        TEMPLATE(T, init)(x2, K);
        TEMPLATE(T, init)(z2, K);
        TEMPLATE(T, init)(x3, K);
        TEMPLATE(T, init)(z3, K);
        TEMPLATE(T, init)(x4, K);
        TEMPLATE(T, init)(z4, K);
        TEMPLATE(T, init)(x5, K);
        TEMPLATE(T, init)(z5, K);

        // P
        TEMPLATE(T, set)(x2, x1, K);
        TEMPLATE(T, set)(z2, z1, K);
        // 2*P
        TEMPLATE(T, weierstrass_xz_dbl)(x3, z3, x1, z1, E, K);

        for (bit = fmpz_bits(m) - 2; bit > 0; bit--)
        {
            if (fmpz_tstbit(m, bit))
            {
                // P, [k+1]P, [k]P -> [2k+2]P, [2k+1]P
                TEMPLATE(T, weierstrass_xz_ladd)(x5, z5, x4, z4, x2, z2, x3, z3, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x5, K);
                TEMPLATE(T, swap)(z2, z5, K);
                TEMPLATE(T, swap)(x3, x4, K);
                TEMPLATE(T, swap)(z3, z4, K);
            }
            else
            {
                // P, [k]P, [k+1]P -> [2k]P, [2k+1]P
                TEMPLATE(T, weierstrass_xz_ladd)(x5, z5, x4, z4, x3, z3, x2, z2, x1, z1, E, K);
                TEMPLATE(T, swap)(x2, x4, K);
                TEMPLATE(T, swap)(z2, z4, K);
                TEMPLATE(T, swap)(x3, x5, K);
                TEMPLATE(T, swap)(z3, z5, K);
            }
        }

        // Last iteration, no need for a ladder
        {
            if (fmpz_tstbit(m, 0))
            {
                TEMPLATE(T, weierstrass_xz_dadd)(x, z, x3, z3, x2, z2, x1, z1, E, K);
            }
            else
            {
                TEMPLATE(T, weierstrass_xz_dbl)(x, z, x2, z2, E, K);
            }
        }

        TEMPLATE(T, clear)(x2, K);
        TEMPLATE(T, clear)(z2, K);
        TEMPLATE(T, clear)(x3, K);
        TEMPLATE(T, clear)(z3, K);
        TEMPLATE(T, clear)(x4, K);
        TEMPLATE(T, clear)(z4, K);
        TEMPLATE(T, clear)(x5, K);
        TEMPLATE(T, clear)(z5, K);
    }
}

#endif
