/**
** @file mruby/bigint.c - Multi-precision Integer
**
** See Copyright Notice in mruby.h
*/

#include <mruby.h>
#include <mruby/object.h>
#include <mruby/numeric.h>
#include <mruby/array.h>
#include <mruby/string.h>
#include <mruby/error.h>
#include <mruby/internal.h>
#include <string.h>
#include "bigint.h"

#define DIG_SIZE (MPZ_DIG_SIZE)
#define DIG_BASE (1ULL << DIG_SIZE)
#define DIG_MASK (DIG_BASE - 1)
#define HIGH(x) ((x) >> DIG_SIZE)
#define LOW(x)  ((x) & DIG_MASK)

#define iabs(x) (((x)>0)?(x):(-x))
#define imax(x,y) (((x)>(y))?(x):(y))
#define imin(x,y) (((x)<(y))?(x):(y))
#define dg(x,i) (((size_t)i < (x)->sz)?(x)->p[i]:0)

typedef mrb_state mpz_ctx_t;
#define MPZ_MRB(ctx) (ctx)
#define MPZ_CTX_INIT(mrb_ptr, ctx, pool_ptr) mrb_state *ctx = (mrb_ptr);

/* Zero n limbs at p */
static inline void
mpn_zero(mp_limb *p, size_t n)
{
  memset(p, 0, n * sizeof(mp_limb));
}

static void
mpz_init(mpz_ctx_t *ctx, mpz_t *s)
{
  s->p = NULL;
  s->sn = 0;
  s->sz = 0;
}

/* Heap-preferred allocation */
static void
mpz_init_heap(mpz_ctx_t *ctx, mpz_t *s, size_t hint)
{
  s->sn = 0;
  if (hint > 0) {
    /* Check for overflow in size calculation (same check as mpz_realloc) */
    if (hint > SIZE_MAX / sizeof(mp_limb)) {
      mrb_state *mrb = MPZ_MRB(ctx);
      mrb_raise(mrb, E_RUNTIME_ERROR, "bigint size too large");
    }
    s->p = (mp_limb*)mrb_malloc(MPZ_MRB(ctx), hint * sizeof(mp_limb));
    mpn_zero(s->p, hint);
    s->sz = hint;
  }
  else {
    s->p = NULL;  /* Lazy allocation via mpz_realloc later */
    s->sz = 0;
  }
}

#define mpz_init_temp(ctx, s, hint) mpz_init_heap(ctx, s, hint)

static void
mpz_realloc(mpz_ctx_t *ctx, mpz_t *x, size_t size)
{
  if (x->sz < size) {
    /* Check for overflow in size calculation */
    if (size > SIZE_MAX / sizeof(mp_limb)) {
      mrb_state *mrb = MPZ_MRB(ctx);
      mrb_raise(mrb, E_RUNTIME_ERROR, "bigint size too large");
    }

    size_t old_sz = x->sz;
    x->p = (mp_limb*)mrb_realloc(MPZ_MRB(ctx), x->p, size * sizeof(mp_limb));
    /* Zero-initialize new limbs */
    mpn_zero(x->p + old_sz, size - old_sz);
    x->sz = size;
  }
}

static void
mpz_set(mpz_ctx_t *ctx, mpz_t *y, mpz_t *x)
{
  size_t i, k = x->sz;

  mpz_realloc(ctx, y, k);
  for (i=0;i < k; i++)
    y->p[i] = x->p[i];

  y->sz = k;
  y->sn = x->sn;
}

static void
mpz_init_set(mpz_ctx_t *ctx, mpz_t *s, mpz_t *t)
{
  mpz_init(ctx, s);
  mpz_set(ctx, s, t);
}

static void
mpz_set_int(mpz_ctx_t *ctx, mpz_t *y, mrb_int v)
{
  mrb_uint u;

  if (v == 0) {
    y->sn = 0;
    y->sz = 0;
    return;
  }
  else if (v > 0) {
    y->sn = 1;
    u = v;
  }
  else /* if (v < 0) */ {
    y->sn = -1;
    if (v == MRB_INT_MIN) u = v;
    else u = -v;
  }
#if MRB_INT_BIT > DIG_SIZE
  if ((u & ~DIG_MASK) != 0) {
    mpz_realloc(ctx, y, 2);
    y->p[1] = (mp_limb)HIGH(u);
    y->p[0] = (mp_limb)LOW(u);
    y->sz = 2;
    return;
  }
#endif
  mpz_realloc(ctx, y, 1);
  y->p[0] = (mp_limb)u;
  y->sz = 1;
}


static void
mpz_set_uint64(mpz_ctx_t *ctx, mpz_t *y, uint64_t u)
{
  size_t len = 0;

  for (uint64_t u0=u; u0; u0>>=DIG_SIZE,len++)
    ;
  y->sn = (u != 0);
  mpz_realloc(ctx, y, len);
  y->sz = len;
  for (size_t i=0; i<len; i++) {
    y->p[i] = (mp_limb)LOW(u);
    u >>= DIG_SIZE;
  }
}

#ifdef MRB_INT32
static void
mpz_set_int64(mpz_ctx_t *ctx, mpz_t *y, int64_t v)
{
  uint64_t u;

  if (v < 0) {
    if (v == INT64_MIN) u = v;
    else u = -v;
  }
  else {
    u = v;
  }
  mpz_set_uint64(ctx, y, u);
  if (v < 0) {
    y->sn = -1;
  }
}
#endif

static void
mpz_init_set_int(mpz_ctx_t *ctx, mpz_t *y, mrb_int v)
{
  mpz_init(ctx, y);
  mpz_set_int(ctx, y, v);
}

static void
mpz_clear(mpz_ctx_t *ctx, mpz_t *s)
{
  if (s->p) {
    mrb_free(MPZ_MRB(ctx), s->p);
    s->p = NULL;
  }
  s->sn = 0;
  s->sz = 0;
}

static void
mpz_move(mpz_ctx_t *ctx, mpz_t *y, mpz_t *x)
{
  /* Normal move: transfer ownership */
  mpz_clear(ctx, y);
  y->sn = x->sn;
  y->sz = x->sz;
  y->p = x->p;
  x->p = NULL;
  x->sn = 0;
  x->sz = 0;
}

static size_t
digits(mpz_t *x)
{
  size_t i;

  if (x->sz == 0) return 0;
  for (i = x->sz - 1; x->p[i] == 0 && i > 0; i--)
    ;
  return i+1;
}

static void
trim(mpz_t *x)
{
  while (x->sz && x->p[x->sz-1] == 0) {
    x->sz--;
  }
}

/* z = x + y, without regard for sign */
/* Core addition algorithm for unsigned operands */
/* Note: mpn_add/mpn_sub are defined later in the file, forward declare here */
static mp_limb mpn_add(mp_limb*, const mp_limb*, size_t, const mp_limb*, size_t);
static mp_limb mpn_sub(mp_limb*, const mp_limb*, size_t, const mp_limb*, size_t);

static void
uadd(mpz_t *z, mpz_t *x, mpz_t *y)
{
  mp_limb carry;
  size_t max_sz = (x->sz > y->sz) ? x->sz : y->sz;

  /* Ensure larger array is first argument to mpn_add */
  if (x->sz >= y->sz) {
    carry = mpn_add(z->p, x->p, x->sz, y->p, y->sz);
  }
  else {
    carry = mpn_add(z->p, y->p, y->sz, x->p, x->sz);
  }

  /* Store final carry */
  z->p[max_sz] = carry;
}

/* z = y - x, ignoring sign */
/* precondition: abs(y) >= abs(x) */
/* Core subtraction algorithm for unsigned operands */
static void
usub(mpz_t *z, mpz_t *y, mpz_t *x)
{
  /* y->sz >= x->sz is guaranteed by precondition */
  mpn_sub(z->p, y->p, y->sz, x->p, x->sz);

  /* Normalize result size */
  z->sz = digits(z);
}

/* Compare two same-length limb arrays: returns <0, 0, or >0 */
static inline int
mpn_cmp(const mp_limb *ap, const mp_limb *bp, size_t n)
{
  while (n-- > 0) {
    if (ap[n] != bp[n]) {
      return (ap[n] > bp[n]) ? 1 : -1;
    }
  }
  return 0;
}

/* compare abs(x) and abs(y) */
static int
ucmp(mpz_t *y, mpz_t *x)
{
  if (y->sz < x->sz) return -1;
  if (y->sz > x->sz) return 1;
  if (x->sz == 0) return 0;
  return mpn_cmp(y->p, x->p, x->sz);
}

#define zero_p(x) ((x)->sn == 0)

/* check if all digits are zero */
static int
uzero_p(mpz_t *x)
{
  if (x->sz == 0) return 1;
  for (size_t i=x->sz-1;; i--) {
    if (x->p[i] != 0)
      return 0;
    if (i == 0) break;
  }
  return 1;
}

static void
zero(mpz_t *x)
{
  x->sn=0;
  if (x->p) {
    x->sz=1;
    x->p[0]=0;
  }
  else {
    x->sz=0;
  }
}

/* z = x + y */
static void
mpz_add(mpz_ctx_t *ctx, mpz_t *zz, mpz_t *x, mpz_t *y)
{
  if (zero_p(x)) {
    mpz_set(ctx, zz, y);
    trim(zz);
    return;
  }
  if (zero_p(y)) {
    mpz_set(ctx, zz, x);
    trim(zz);
    return;
  }

  /* Fast path: single-limb + multi-limb */
  if (y->sz == 1 && x->sz > 1) {
    mp_limb y_limb = y->p[0];
    mpz_t z;
    mpz_init_heap(ctx, &z, x->sz + 1);

    if ((x->sn > 0 && y->sn > 0) || (x->sn < 0 && y->sn < 0)) {
      /* Same signs: addition */
      mp_dbl_limb carry = y_limb;
      carry += x->p[0];
      z.p[0] = (mp_limb)carry;
      carry >>= DIG_SIZE;

      /* Propagate carry through remaining limbs */
      for (size_t i = 1; i < x->sz; i++) {
        carry += x->p[i];
        z.p[i] = (mp_limb)carry;
        carry >>= DIG_SIZE;
      }
      z.p[x->sz] = (mp_limb)carry;
      z.sn = x->sn;
    }
    else {
      /* Different signs: subtraction */
      if (x->sz == 1 && y_limb == x->p[0]) {
        /* Equal magnitude: result is zero */
        zero(&z);
      }
      else if (x->sz == 1 && x->p[0] > y_limb) {
        /* |x| > |y|: result has sign of x */
        z.p[0] = x->p[0] - y_limb;
        z.p[1] = 0;
        z.sn = x->sn;
      }
      else {
        /* |x| > |y|: subtract y from x */
        mp_dbl_limb borrow = y_limb;
        if (x->p[0] >= borrow) {
          z.p[0] = x->p[0] - (mp_limb)borrow;
          borrow = 0;
        }
        else {
          z.p[0] = (mp_limb)(((mp_dbl_limb)1 << DIG_SIZE) + x->p[0] - (mp_limb)borrow);
          borrow = 1;
        }

        /* Propagate borrow through remaining limbs */
        for (size_t i = 1; i < x->sz; i++) {
          if (x->p[i] >= borrow) {
            z.p[i] = x->p[i] - (mp_limb)borrow;
            borrow = 0;
          }
          else {
            z.p[i] = (mp_limb)(((mp_dbl_limb)1 << DIG_SIZE) + x->p[i] - (mp_limb)borrow);
            borrow = 1;
          }
        }
        z.sn = x->sn;
      }
    }
    trim(&z);
    mpz_move(ctx, zz, &z);
    return;
  }

  if (x->sz == 1 && y->sz > 1) {
    /* Swap and use the same fast path */
    mpz_add(ctx, zz, y, x);
    return;
  }

  mpz_t z;
  size_t estimated_size = ((x->sz > y->sz) ? x->sz : y->sz) + 1;
  mpz_init_heap(ctx, &z, estimated_size);

  if (x->sn > 0 && y->sn > 0) {
    uadd(&z, x, y);
    z.sn = 1;
  }
  else if (x->sn < 0 && y->sn < 0) {
    uadd(&z, x, y);
    z.sn = -1;
  }
  else {
    int mg;

    /* signs differ */
    if ((mg = ucmp(x,y)) == 0) {
      zero(&z);
    }
    else if (mg > 0) {  /* abs(y) < abs(x) */
      usub(&z, x, y);
      z.sn = (x->sn > 0 && y->sn < 0) ? 1 : (-1);
    }
    else { /* abs(y) > abs(x) */
      usub(&z, y, x);
      z.sn = (x->sn < 0 && y->sn > 0) ? 1 : (-1);
    }
  }
  trim(&z);
  mpz_move(ctx, zz, &z);
}

/* x += n                                              */
/*   ignores sign of x                                 */
/*   assumes n is positive and small (fits in mp_limb) */
static void
mpz_add_int(mpz_ctx_t *ctx, mpz_t *x, mrb_int n)
{
  // If n is zero, no operation is needed
  if (n == 0) return;

  // Assume x is positive and n is a small positive integer
  mp_dbl_limb carry = n; // Initialize carry with n
  for (size_t i = 0; i < x->sz && carry; i++) {
    carry += (mp_dbl_limb)x->p[i]; // Add current limb and carry
    x->p[i] = LOW(carry);          // Store lower 32 bits in current limb
    carry = HIGH(carry);           // Update carry with higher bits
  }

  if (carry != 0) {
    mpz_realloc(ctx, x, x->sz + 1);
    x->p[x->sz-1] = (mp_limb)carry;
    x->sn = 1;
  }
  trim(x);
}

/* z = x - y */
static void
mpz_sub(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mpz_t *y)
{
  /* In-place optimization: z == x, both positive, x >= y */
  if (z == x && x->sn > 0 && y->sn > 0 && ucmp(x, y) >= 0) {
    mpn_sub(x->p, x->p, x->sz, y->p, y->sz);
    x->sz = digits(x);
    if (x->sz == 0) x->sn = 0;
    return;
  }

  /* General case: create view of y with negated sign and use mpz_add */
  mpz_t u;
  u.p = y->p;
  u.sz = y->sz;
  u.sn = -(y->sn);
  mpz_add(ctx, z, x, &u);
}

/* x -= n                                              */
/*   ignores sign of x                                 */
/*   assumes n is positive and small (fits in mp_limb) */
static void
mpz_sub_int(mpz_ctx_t *ctx, mpz_t *x, mrb_int n)
{
  // If n is zero, no operation is needed
  if (n == 0) return;

  // If x is zero, set x to n
  if (zero_p(x) || x->sz == 0) {
    mpz_set_int(ctx, x, n);
    return;
  }

  // Initialize borrow and start decrement
  mp_dbl_limb_signed borrow = (mp_limb)n;
  size_t i = 0;

  // Subtract 1 from the least significant limb and propagate if necessary
  borrow = (mp_dbl_limb_signed)x->p[i] - borrow;
  x->p[i] = LOW(borrow);
  borrow = (borrow < 0) ? 1 : 0;

  // Continue through limbs while there is a borrow
  for (i = 1; i < x->sz && borrow; i++) {
    borrow = (mp_dbl_limb_signed)x->p[i] - borrow;
    x->p[i] = LOW(borrow);
    borrow = (borrow < 0) ? 1 : 0;
  }

  // Trim any unnecessary leading zeros
  trim(x);
}

/* Multiply-and-add: rp[0..n-1] += s1p[0..n-1] * limb; return carry (high limb) */
static inline mp_limb
mpn_addmul_1(mp_limb *rp, const mp_limb *s1p, size_t n, mp_limb limb)
{
  mp_dbl_limb acc = 0;
  for (size_t i = 0; i < n; i++) {
    acc += (mp_dbl_limb)rp[i] + (mp_dbl_limb)s1p[i] * (mp_dbl_limb)limb;
    rp[i] = LOW(acc);
    acc = HIGH(acc);
  }
  return (mp_limb)acc;
}

/* w = u * v (optimized schoolbook using mpn_addmul_1) */
static void
mpz_mul(mpz_ctx_t *ctx, mpz_t *ww, mpz_t *u, mpz_t *v)
{
  if (zero_p(u) || zero_p(v)) {
    zero(ww);
    return;
  }

  /* Ensure outer loop iterates over the shorter operand for better cache use */
  mpz_t *a, *b;
  if (v->sz > u->sz) {
    a = v; b = u;
  }
  else {
    a = u; b = v;
  }

  /* Fast path: single-limb * multi-limb */
  if (b->sz == 1) {
    mp_limb scalar = b->p[0];
    mpz_t w;
    mpz_init_heap(ctx, &w, a->sz + 1);
    mpn_zero(w.p, a->sz + 1);

    mp_limb carry = mpn_addmul_1(w.p, a->p, a->sz, scalar);
    w.p[a->sz] = carry;

    w.sn = a->sn * b->sn;
    trim(&w);
    mpz_move(ctx, ww, &w);
    return;
  }

  mpz_t w;
  mpz_init_heap(ctx, &w, a->sz + b->sz);
  mpn_zero(w.p, a->sz + b->sz);

  for (size_t j = 0; j < a->sz; j++) {
    mp_limb a_limb = a->p[j];
    if (a_limb == 0) continue;

    mp_limb carry = mpn_addmul_1(w.p + j, b->p, b->sz, a_limb);

    /* Properly handle carry propagation to avoid overflow */
    size_t k = j + b->sz;
    while (carry && k < a->sz + b->sz) {
      mp_dbl_limb sum = (mp_dbl_limb)w.p[k] + (mp_dbl_limb)carry;
      w.p[k] = LOW(sum);
      carry = HIGH(sum);
      k++;
    }
  }

  w.sn = a->sn * b->sn;
  trim(&w);
  mpz_move(ctx, ww, &w);
}

/* Count set bits in a limb */
static inline int
limb_popcount(mp_limb x)
{
#if defined(__GNUC__) || __has_builtin(__builtin_popcount)
  if (sizeof(mp_limb) == sizeof(unsigned long long))
    return __builtin_popcountll(x);
  else
    return __builtin_popcount(x);
#else
  int count = 0;
  while (x) {
    count++;
    x &= x - 1;  /* Clear lowest set bit */
  }
  return count;
#endif
}

/* number of leading zero bits in digit */
static int
lzb(mp_limb x)
{
  if (x == 0) return 0;
#if (defined(__GNUC__) || __has_builtin(__builtin_clz))
  if (sizeof(mp_limb) == sizeof(int64_t))
    return __builtin_clzll(x);
  else if (sizeof(mp_limb) == sizeof(int32_t))
    return __builtin_clz(x);
#endif

  int j=0;

  for (mp_limb i = ((mp_limb)1 << (DIG_SIZE-1)); i && !(x&i); j++,i>>=1)
    ;
  return j;
}

/*
 * mpn-style low-level shift functions.
 * These operate directly on limb arrays without mpz_t overhead.
 */

/*
 * Right shift limb array by cnt bits (0 < cnt < DIG_SIZE).
 * rp[0..n-1] = ap[0..n-1] >> cnt
 * Returns the bits shifted out from the low end.
 * Supports in-place operation (rp == ap).
 *
 * Processing order: LOW to HIGH
 * - rp[i] depends on ap[i] and ap[i+1]
 * - Writing rp[i] doesn't affect ap[i+1], so in-place is safe
 */
static mp_limb
mpn_rshift(mp_limb *rp, const mp_limb *ap, size_t n, unsigned int cnt)
{
  mp_limb shifted_out;
  size_t i;

  mrb_assert(cnt > 0 && cnt < DIG_SIZE);

  shifted_out = (ap[0] << (DIG_SIZE - cnt)) & DIG_MASK;

  for (i = 0; i < n - 1; i++) {
    rp[i] = ((ap[i] >> cnt) | (ap[i + 1] << (DIG_SIZE - cnt))) & DIG_MASK;
  }
  rp[n - 1] = (ap[n - 1] >> cnt) & DIG_MASK;

  return shifted_out;
}

/*
 * Left shift limb array by cnt bits (0 < cnt < DIG_SIZE).
 * rp[0..n-1] = ap[0..n-1] << cnt (lower n limbs)
 * Returns the bits shifted out from the high end (carry).
 * Supports in-place operation (rp == ap).
 *
 * Processing order: HIGH to LOW
 * - rp[i] depends on ap[i] and ap[i-1]
 * - Writing rp[i] doesn't affect ap[i-1], so in-place is safe
 */
static mp_limb
mpn_lshift(mp_limb *rp, const mp_limb *ap, size_t n, unsigned int cnt)
{
  mp_limb carry;
  size_t i;

  mrb_assert(cnt > 0 && cnt < DIG_SIZE);

  carry = (ap[n - 1] >> (DIG_SIZE - cnt)) & DIG_MASK;

  for (i = n - 1; i > 0; i--) {
    rp[i] = ((ap[i] << cnt) | (ap[i - 1] >> (DIG_SIZE - cnt))) & DIG_MASK;
  }
  rp[0] = (ap[0] << cnt) & DIG_MASK;

  return carry;
}

/*
 * Add two limb arrays of the same size.
 * rp[0..n-1] = ap[0..n-1] + bp[0..n-1]
 * Returns carry (0 or 1).
 * Supports in-place operation (rp == ap or rp == bp).
 */
static inline mp_limb
mpn_add_n(mp_limb *rp, const mp_limb *ap, const mp_limb *bp, size_t n)
{
  mp_dbl_limb carry = 0;

  for (size_t i = 0; i < n; i++) {
    carry += (mp_dbl_limb)ap[i] + (mp_dbl_limb)bp[i];
    rp[i] = LOW(carry);
    carry = HIGH(carry);
  }
  return (mp_limb)carry;
}

/*
 * Add two limb arrays of different sizes.
 * rp[0..an-1] = ap[0..an-1] + bp[0..bn-1]
 * Precondition: an >= bn
 * Returns carry (0 or 1).
 * Supports in-place operation (rp == ap).
 */
static mp_limb
mpn_add(mp_limb *rp, const mp_limb *ap, size_t an, const mp_limb *bp, size_t bn)
{
  mp_dbl_limb carry;

  mrb_assert(an >= bn);

  /* Add overlapping part */
  carry = mpn_add_n(rp, ap, bp, bn);

  /* Propagate carry through remaining limbs of ap */
  for (size_t i = bn; i < an; i++) {
    carry += (mp_dbl_limb)ap[i];
    rp[i] = LOW(carry);
    carry = HIGH(carry);
  }
  return (mp_limb)carry;
}

/*
 * Subtract two limb arrays of the same size.
 * rp[0..n-1] = ap[0..n-1] - bp[0..n-1]
 * Returns borrow (0 or 1).
 * Supports in-place operation (rp == ap).
 */
static inline mp_limb
mpn_sub_n(mp_limb *rp, const mp_limb *ap, const mp_limb *bp, size_t n)
{
  mp_dbl_limb_signed borrow = 0;

  for (size_t i = 0; i < n; i++) {
    borrow += (mp_dbl_limb_signed)ap[i] - (mp_dbl_limb_signed)bp[i];
    rp[i] = LOW(borrow);
    borrow = HIGH(borrow);
  }
  return (mp_limb)(-borrow);
}

/*
 * Subtract two limb arrays of different sizes.
 * rp[0..an-1] = ap[0..an-1] - bp[0..bn-1]
 * Precondition: an >= bn and ap >= bp (result is non-negative)
 * Returns borrow (0 or 1, should be 0 if precondition met).
 * Supports in-place operation (rp == ap).
 */
static mp_limb
mpn_sub(mp_limb *rp, const mp_limb *ap, size_t an, const mp_limb *bp, size_t bn)
{
  mp_dbl_limb_signed borrow;

  mrb_assert(an >= bn);

  /* Subtract overlapping part */
  borrow = -(mp_dbl_limb_signed)mpn_sub_n(rp, ap, bp, bn);

  /* Propagate borrow through remaining limbs of ap */
  for (size_t i = bn; i < an; i++) {
    borrow += (mp_dbl_limb_signed)ap[i];
    rp[i] = LOW(borrow);
    borrow = HIGH(borrow);
  }
  return (mp_limb)(-borrow);
}

/* c1 = a>>n */
/* n must be < DIG_SIZE */
static void
urshift(mpz_ctx_t *ctx, mpz_t *c1, mpz_t *a, size_t n)
{
  mrb_assert(n < DIG_SIZE);

  if (n == 0) {
    mpz_set(ctx, c1, a);
    trim(c1);
  }
  else if (uzero_p(a)) {
    zero(c1);
  }
  else {
    mpz_realloc(ctx, c1, a->sz);
    mpn_rshift(c1->p, a->p, a->sz, (unsigned int)n);
    c1->sz = a->sz;
    trim(c1);
  }
}

/* c1 = a<<n */
/* n must be < DIG_SIZE */
static void
ulshift(mpz_ctx_t *ctx, mpz_t *c1, mpz_t *a, size_t n)
{
  mrb_assert(n < DIG_SIZE);
  if (n == 0) {
    mpz_set(ctx, c1, a);
    trim(c1);
  }
  else if (uzero_p(a)) {
    zero(c1);
  }
  else if (c1 == a) {
    /* In-place optimization: mpn_lshift works from high to low, safe for aliasing */
    mp_limb carry;
    size_t old_sz = a->sz;

    mpz_realloc(ctx, c1, old_sz + 1);
    carry = mpn_lshift(c1->p, c1->p, old_sz, (unsigned int)n);
    c1->p[old_sz] = carry;
    c1->sz = old_sz + 1;
    trim(c1);
  }
  else {
    mpz_t c;
    mp_limb carry;

    mpz_init_heap(ctx, &c, a->sz + 1);
    carry = mpn_lshift(c.p, a->p, a->sz, (unsigned int)n);
    c.p[a->sz] = carry;
    c.sz = a->sz + 1;
    trim(&c);
    mpz_move(ctx, c1, &c);
  }
}

/* Fast division by single limb */
static void
div_limb(mpz_ctx_t *ctx, mpz_t *q, mpz_t *r, mpz_t *x, mp_limb d)
{
  mrb_state *mrb = MPZ_MRB(ctx);
  mpz_t temp_q, temp_r;
  size_t n;
  mp_dbl_limb remainder;

  if (zero_p(x)) {
    zero(q);
    zero(r);
    return;
  }

  if (d == 0) {
    mrb_raise(mrb, E_ZERODIV_ERROR, "divided by 0");
  }

  /* Power-of-2 divisor optimization */
  if ((d & (d - 1)) == 0) {
    /* d is power of 2, use bit operations */
    int shift = 0;
    mp_limb temp = d;
    while (temp > 1) {
      temp >>= 1;
      shift++;
    }

    /* Quotient = x >> shift */
    if (shift == 0) {
      mpz_init(ctx, &temp_q);
      mpz_init(ctx, &temp_r);
      mpz_set(ctx, &temp_q, x);
    }
    else {
      /* Manual right shift implementation */
      size_t limb_shift = shift / DIG_SIZE;
      size_t bit_shift = shift % DIG_SIZE;

      if (limb_shift >= x->sz) {
        mpz_init(ctx, &temp_q);
        mpz_init(ctx, &temp_r);
        zero(&temp_q);
      }
      else {
        size_t new_size = x->sz - limb_shift;
        mpz_init_temp(ctx, &temp_q, new_size);
        mpz_init(ctx, &temp_r);

        if (bit_shift == 0) {
          /* Simple limb copy */
          for (size_t i = 0; i < new_size; i++) {
            temp_q.p[i] = x->p[i + limb_shift];
          }
        }
        else {
          /* Bit shift within limbs */
          mp_limb carry = 0;
          for (size_t i = new_size; i > 0; i--) {
            mp_limb current = x->p[i - 1 + limb_shift];
            temp_q.p[i - 1] = (current >> bit_shift) | carry;
            carry = (current << (DIG_SIZE - bit_shift)) & DIG_MASK;
          }
        }
        temp_q.sz = new_size;
        trim(&temp_q);
        temp_q.sn = (temp_q.sz == 0) ? 0 : 1;
      }
    }

    /* Remainder = x & (d - 1) */
    /* temp_r is already initialized in all code paths above */
    mpz_realloc(ctx, &temp_r, 1);
    temp_r.p[0] = x->p[0] & (d - 1);
    temp_r.sz = (temp_r.p[0] == 0) ? 0 : 1;
    temp_r.sn = (temp_r.sz == 0) ? 0 : 1;
    mpz_move(ctx, q, &temp_q);
    mpz_move(ctx, r, &temp_r);
    return;
  }

  /* General single-limb division */
  if (x->sz == 1) {
    /* Both dividend and divisor are single limb */
    mpz_init_temp(ctx, &temp_q, 1);
    mpz_init_temp(ctx, &temp_r, 1);

    temp_q.p[0] = x->p[0] / d;
    temp_r.p[0] = x->p[0] % d;

    temp_q.sz = (temp_q.p[0] == 0) ? 0 : 1;
    temp_q.sn = (temp_q.sz == 0) ? 0 : 1;

    temp_r.sz = (temp_r.p[0] == 0) ? 0 : 1;
    temp_r.sn = (temp_r.sz == 0) ? 0 : 1;
    mpz_move(ctx, q, &temp_q);
    mpz_move(ctx, r, &temp_r);
    return;
  }

  /* Multi-limb dividend, single-limb divisor */
  n = x->sz;
  mpz_init_temp(ctx, &temp_q, n);
  mpz_init_temp(ctx, &temp_r, 1);

  remainder = 0;

  /* Process from most significant limb to least significant */
  for (size_t i = n; i > 0; i--) {
    remainder = (remainder << DIG_SIZE) + x->p[i-1];
    temp_q.p[i-1] = (mp_limb)(remainder / d);
    remainder = remainder % d;
  }

  /* Set remainder */
  temp_r.p[0] = (mp_limb)remainder;
  temp_r.sz = (remainder == 0) ? 0 : 1;
  temp_r.sn = (temp_r.sz == 0) ? 0 : 1;

  /* Trim leading zeros from quotient */
  trim(&temp_q);
  temp_q.sn = (temp_q.sz == 0) ? 0 : 1;

  /* Copy results to avoid pool/heap mixing */
  mpz_move(ctx, q, &temp_q);
  mpz_move(ctx, r, &temp_r);
}

static void
udiv(mpz_ctx_t *ctx, mpz_t *qq, mpz_t *rr, mpz_t *xx, mpz_t *yy)
{
  /* Handle simple cases */
  int cmp = ucmp(xx, yy);
  if (cmp == 0) {
    mpz_set_int(ctx, qq, 1);
    zero(rr);
    return;
  }
  else if (cmp < 0) {
    zero(qq);
    mpz_set(ctx, rr, xx);
    return;
  }

  /* Fast path for single-limb divisor */
  if (yy->sz == 1) {
    div_limb(ctx, qq, rr, xx, yy->p[0]);
    return;
  }

  mrb_assert(yy->sn != 0);      /* divided by zero */
  mrb_assert(yy->sz > 0);       /* divided by zero */
  mrb_assert(!uzero_p(yy));     /* divided by zero */

  /* Pre-check size constraints to avoid memory leaks on exception.
   * Check both dividend and divisor sizes as ulshift uses size+1 for both.
   * Fail early before any allocation to prevent leaks. */
  size_t max_mpz_size = SIZE_MAX / sizeof(mp_limb);
  if (xx->sz + 1 > max_mpz_size || yy->sz + 1 > max_mpz_size) {
    mrb_state *mrb = MPZ_MRB(ctx);
    mrb_raise(mrb, E_RUNTIME_ERROR, "bigint size too large");
  }

  mpz_t q, x, y;
  mpz_init_temp(ctx, &q, xx->sz - yy->sz + 1);  /* Quotient size estimate */
  mpz_init_temp(ctx, &x, xx->sz + 1);           /* Dividend with potential carry */
  mpz_init_temp(ctx, &y, yy->sz);               /* Divisor copy */
  mpz_realloc(ctx, &x, xx->sz+1);
  size_t yd = digits(yy);
  size_t ns = lzb(yy->p[yd-1]);
  ulshift(ctx, &x, xx, ns);
  ulshift(ctx, &y, yy, ns);
  trim(&y);  /* Trim after shift to remove any zero limbs */
  size_t xd = digits(&x);
  yd = digits(&y);

  /* Handle edge case: divisor became zero after normalization */
  if (yd == 0 || y.p[yd-1] == 0) {
    /* This should not happen with valid inputs, but handle gracefully */
    zero(qq);
    zero(rr);
    mpz_clear(ctx, &q);
    mpz_clear(ctx, &x);
    mpz_clear(ctx, &y);
    return;
  }

  mpz_realloc(ctx, &q, xd-yd+1);  // Quotient has xd-yd+1 digits maximum

  /* Core Knuth Algorithm D division loop */
  mp_dbl_limb z = y.p[yd-1];
  mrb_assert(z != 0);  /* Divisor high limb must be non-zero after normalization */

  if (xd >= yd) {
    for (size_t j = xd - yd;; j--) {
      mp_dbl_limb qhat;
      mp_dbl_limb rhat;

      if (j + yd == xd) {
        /* Only one high limb available */
        mp_dbl_limb dividend_val = (((mp_dbl_limb)0 << DIG_SIZE) + x.p[j+yd-1]);
        qhat = dividend_val / z;
        rhat = dividend_val % z;
      }
      else {
        /* Two limbs available - standard Knuth estimation */
        mp_dbl_limb dividend_val = ((mp_dbl_limb)x.p[j+yd] << DIG_SIZE) + x.p[j+yd-1];
        qhat = dividend_val / z;
        rhat = dividend_val % z;
      }

      /* Standard Knuth Algorithm D qhat refinement (2-limb check) */
      if (yd >= 2) {
        mp_dbl_limb y_second = y.p[yd-2];
        mp_dbl_limb x_third = (j+yd-2 < x.sz) ? x.p[j+yd-2] : 0;
        mp_dbl_limb left_side = qhat * y_second;
        mp_dbl_limb right_side = (rhat << DIG_SIZE) + x_third;

        while (qhat >= ((mp_dbl_limb)1 << DIG_SIZE) || (left_side > right_side)) {
          qhat--;
          rhat += z;
          if (rhat >= ((mp_dbl_limb)1 << DIG_SIZE)) break;
          left_side -= y_second;
          right_side = (rhat << DIG_SIZE) + x_third;
        }
      }

      if (qhat > 0) {
        /* Subtract qhat * divisor from dividend */
        /* Use pointers to avoid repeated i+j index calculation */
        mp_dbl_limb_signed borrow = 0;
        mp_limb *xp = x.p + j;
        mp_limb *xp_end = x.p + x.sz;
        const mp_limb *yp = y.p;

        for (size_t i = 0; i < yd; i++) {
          mp_dbl_limb product = qhat * *yp++;
          mp_dbl_limb_signed diff = (mp_dbl_limb_signed)*xp - (mp_dbl_limb_signed)LOW(product) + borrow;
          *xp++ = LOW(diff);
          borrow = HIGH(diff) - (mp_dbl_limb_signed)HIGH(product);
        }

        /* Handle final borrow propagation */
        if (xp < xp_end) {
          borrow += (mp_dbl_limb_signed)*xp;
          *xp = LOW(borrow);
          borrow = HIGH(borrow);
        }

        /* Correction: if borrow is negative, qhat was too large, add back */
        if (borrow < 0) {
          qhat--;
          mp_dbl_limb carry = 0;
          xp = x.p + j;
          yp = y.p;
          for (size_t i = 0; i < yd; i++) {
            carry += (mp_dbl_limb)*xp + (mp_dbl_limb)*yp++;
            *xp++ = LOW(carry);
            carry = HIGH(carry);
          }
          if (xp < xp_end && carry > 0) {
            *xp += (mp_limb)carry;
          }
        }
      }

      q.p[j] = (mp_limb)qhat;
      if (j == 0) break;
    }
  }
  x.sz = yd;
  urshift(ctx, rr, &x, ns);
  trim(&q);
  mpz_move(ctx, qq, &q);
  mpz_clear(ctx, &q);
  mpz_clear(ctx, &x);
  mpz_clear(ctx, &y);
}

static void
mpz_mdiv(mpz_ctx_t *ctx, mpz_t *q, mpz_t *x, mpz_t *y)
{
  mpz_t r;
  short sn1 = x->sn, sn2 = y->sn, qsign;

  if (zero_p(x)) {
    mpz_init_set_int(ctx, q, 0);
    return;
  }
  mpz_init(ctx, &r);
  udiv(ctx, q, &r, x, y);
  qsign = q->sn = sn1 * sn2;
  if (uzero_p(q))
    q->sn = 0;
  /* now if r != 0 and q < 0 we need to round q towards -inf */
  if (!uzero_p(&r) && qsign < 0) {
    /* add 1 to magnitude */
    mpz_add_int(ctx, q, 1);
    /* force negative sign in case the value of q was zero before rounding */
    q->sn = -1;
  }
  mpz_clear(ctx, &r);
}

static void
mpz_mmod(mpz_ctx_t *ctx, mpz_t *r, mpz_t *x, mpz_t *y)
{
  mpz_t q;
  short sn1 = x->sn, sn2 = y->sn, sn3;

  mpz_init(ctx, &q);
  if (sn1 == 0) {
    zero(r);
    return;
  }
  udiv(ctx, &q, r, x, y);
  mpz_clear(ctx, &q);
  if (uzero_p(r)) {
    r->sn = 0;
    return;
  }
  sn3 = sn1 * sn2;
  if (sn3 > 0)
    r->sn = sn1;
  else if (sn1 < 0 && sn2 > 0) {
    r->sn = 1;
    mpz_sub(ctx, r, y, r);
  }
  else {
    r->sn = 1;
    mpz_add(ctx, r, y, r);
  }
}

static void
mpz_mdivmod(mpz_ctx_t *ctx, mpz_t *q, mpz_t *r, mpz_t *x, mpz_t *y)
{
  short sn1 = x->sn, sn2 = y->sn, qsign;

  if (sn1 == 0) {
    zero(q);
    zero(r);
    return;
  }
  udiv(ctx, q, r, x, y);
  qsign = q->sn = sn1 * sn2;
  if (uzero_p(r)) {
    /* q != 0, since q=r=0 would mean x=0, which was tested above */
    r->sn = 0;
    return;
  }
  if (q->sn > 0)
    r->sn = sn1;
  else if (sn1 < 0 && sn2 > 0) {
    r->sn = 1;
    mpz_sub(ctx, r, y, r);
  }
  else {
    r->sn = 1;
    mpz_add(ctx, r, y, r);
  }
  if (uzero_p(q))
    q->sn = 0;
  /* now if r != 0 and q < 0 we need to round q towards -inf */
  if (!uzero_p(r) && qsign < 0) {
    /* add 1 to magnitude */
    mpz_add_int(ctx, q, 1);
    /* force negative sign in case the value of q was zero before rounding */
    q->sn = -1;
  }
}

/* Fast modular reduction for single-limb modulus */
static void
mpz_mod_limb(mpz_ctx_t *ctx, mpz_t *r, mpz_t *x, mp_limb m)
{
  if (zero_p(x)) {
    zero(r);
    return;
  }

  if (x->sz == 1) {
    /* Single limb case - simple modulo */
    mp_limb result = x->p[0] % m;
    mpz_set_int(ctx, r, result);
    if (result == 0)
      r->sn = 0;
    else
      r->sn = x->sn;
    return;
  }

  /* Multi-limb case - use repeated division */
  mp_dbl_limb remainder = 0;
  for (size_t i = x->sz; i > 0; i--) {
    remainder = (remainder << DIG_SIZE) | x->p[i-1];
    remainder %= m;
  }

  mpz_set_int(ctx, r, (mp_limb)remainder);
  r->sn = x->sn;
  if (remainder == 0)
    r->sn = 0;
}

static void
mpz_mod(mpz_ctx_t *ctx, mpz_t *r, mpz_t *x, mpz_t *y)
{
  short sn = x->sn;

  if (zero_p(x)) {
    zero(r);
    return;
  }

  /* Fast path for single-limb modulus */
  if (y->sz == 1) {
    mpz_mod_limb(ctx, r, x, y->p[0]);
    if (y->sn < 0) r->sn = -r->sn;
    return;
  }

  /* General division fallback */
  mpz_t q;
  mpz_init_temp(ctx, &q, x->sz);
  mpz_realloc(ctx, r, y->sz);
  udiv(ctx, &q, r, x, y);
  r->sn = sn;
  if (uzero_p(r))
    r->sn = 0;
  mpz_clear(ctx, &q);
}

static mrb_int
mpz_cmp(mpz_ctx_t *ctx, mpz_t *x, mpz_t *y)
{
  if (x->sn < 0 && y->sn > 0)
    return (-1);
  if (x->sn > 0 && y->sn < 0)
    return 1;
  int abscmp=ucmp(x, y);
  if (x->sn >=0 && y->sn >=0)
    return abscmp;
  return (-abscmp);          // if (x->sn <=0 && y->sn <=0)
}

/* 2<=base<=36 - this overestimates the optimal value, which is OK */
static size_t
mpz_sizeinbase(mpz_t *x, mrb_int base)
{
  size_t i, j;

  size_t bits = digits(x) * DIG_SIZE;
  mrb_assert(2 <= base && base <= 36);

  if (zero_p(x) || x->sz == 0) return 0;
  for (j=0,i=1; i<=(size_t)base; i*=2,j++)
    ;
  return bits/(j-1)+1;
}

/* x = y * n (only called from mpz_init_set_str) */
/*   assumes x and n are positive or zero        */
/*   assumes n is small (fits in mp_limb)        */
static void
mpz_mul_int(mpz_ctx_t *ctx, mpz_t *x, mrb_int n)
{
  if (n == 0 || zero_p(x)) {
    zero(x);
    return;
  }

  size_t x_sz = x->sz;
  size_t new_sz = x_sz + 1; // Maximum possible size after multiplication

  // Reallocate x if necessary
  mpz_realloc(ctx, x, new_sz);

  mp_dbl_limb cc = 0;
  mp_limb n_limb = (mp_limb)n;

  for (size_t i = 0; i < x_sz; i++) {
    // Multiply each limb and add carry
    cc += (mp_dbl_limb)x->p[i] * n_limb;
    x->p[i] = LOW(cc);
    cc = HIGH(cc);
  }

  if (cc) {
    // If there is a remaining carry, store it and update size
    x->p[x_sz] = (mp_limb)cc;
    x->sz = x_sz + 1;
  }
  else {
    x->sz = x_sz;
  }

  x->sn = 1;
  trim(x);
}

static int
mpz_init_set_str(mpz_ctx_t *ctx, mpz_t *x, const char *s, mrb_int len, mrb_int base)
{
  int retval = 0;
  short sn;
  uint8_t k;

  mpz_init(ctx, x);
  if (*s == '-') {
    sn = -1; s++; len--;
  }
  else if (*s == '+') {
    sn = 1; s++; len--;
  }
  else
    sn = 1;

  for (mrb_int i = 0; i < len; i++) {
    if (s[i]=='_') continue;
    if (s[i] >= '0' && s[i] <= '9')
      k = (uint8_t)s[i] - (uint8_t)'0';
    else if (s[i] >= 'A' && s[i] <= 'Z')
      k = (uint8_t)s[i] - (uint8_t)'A'+10;
    else if (s[i] >= 'a' && s[i] <= 'z')
      k = (uint8_t)s[i] - (uint8_t)'a'+10;
    else {
      retval = (-1);
      break;
    }
    if (k >= base) {
      retval = (-1);
      break;
    }
    mpz_mul_int(ctx, x, base);
    mpz_add_int(ctx, x, k);
  }
  x->sn = x->sz == 0 ? 0 : sn;

  return retval;
}

/* power of base no bigger than DIG_BASE */
/* power of 2 is handled differently */
static const mp_limb base_limit[34*2] = {
#ifdef MRB_NO_MPZ64BIT
  59049,        // 3^10
  0,            // 4^8 (skip)
  15625,        // 5^6
  46656,        // 6^6
  16807,        // 7^5
  0,            // 8^5 (skip)
  59049,        // 9^5
  10000,        // 10^4
  14641,        // 11^4
  20736,        // 12^4
  28561,        // 13^4
  38416,        // 14^4
  50625,        // 15^4
  0,            // 16^4 (skip)
  4913,         // 17^3
  5832,         // 18^3
  6859,         // 19^3
  8000,         // 20^3
  9261,         // 21^3
  10648,        // 22^3
  12167,        // 23^3
  13824,        // 24^3
  15625,        // 25^3
  17576,        // 26^3
  19683,        // 27^3
  21952,        // 28^3
  24389,        // 29^3
  27000,        // 30^3
  29791,        // 31^3
  0,            // 32^3 (skip)
  35937,        // 33^3
  39304,        // 34^3
  42875,        // 35^3
  46656,        // 36^3
#else
  3486784401UL, // 3^20
  0,            // 4^16 (skip)
  1220703125UL, // 5^13
  2176782336UL, // 6^12
  1977326743UL, // 7^11
  0,            // 8^10 (skip)
  3486784401UL, // 9^10
  1000000000UL, // 10^9
  2357947691UL, // 11^9
  429981696UL,  // 12^8
  815730721UL,  // 13^8
  1475789056UL, // 14^8
  2562890625UL, // 15^8
  0,            // 16^8 (skip)
  410338673UL,  // 17^7
  612220032UL,  // 18^7
  893871739UL,  // 19^7
  1280000000UL, // 20^7
  1801088541UL, // 21^7
  2494357888UL, // 22^7
  3404825447UL, // 23^7
  191102976UL,  // 24^6
  244140625UL,  // 25^6
  308915776UL,  // 26^6
  387420489UL,  // 27^6
  481890304UL,  // 28^6
  594823321UL,  // 29^6
  729000000UL,  // 30^6
  887503681UL,  // 31^6
  0,            // 32^6 (skip)
  1291467969UL, // 33^6
  1544804416UL, // 34^6
  1838265625UL, // 35^6
  2176782336UL, // 36^6
#endif
};

static char*
mpz_get_str(mpz_ctx_t *ctx, char *s, mrb_int sz, mrb_int base, mpz_t *x)
{
  mrb_state *mrb = MPZ_MRB(ctx);

  mrb_assert(2 <= base && base <= 36);
  if (zero_p(x)) {
    *s='0';
    *(s+1)='\0';
    return s;
  }

  char *ps = s;
  char *se = s+sz;
  int xlen = (int)digits(x);

  if ((base & (base - 1)) == 0) {  // base is a power of 2
    int shift = 0;
    while (((uint64_t)1 << shift) < (uint64_t)base) shift++;
    mp_limb mask = (mp_limb)base - 1;
    mp_dbl_limb value = 0;
    int bits = 0;

    /* Process all limbs */
    for (int i = 0; i < xlen; i++) {
      value |= (mp_dbl_limb)x->p[i] << bits;
      bits += DIG_SIZE;
      while (bits >= shift) {
        mp_limb digit = value & mask;
        value >>= shift;
        bits -= shift;

        if (digit < 10) *s++ = '0' + digit;
        else *s++ = 'a' + digit - 10;
      }
    }

    /* Handle any remaining bits */
    while (bits > 0) {
      mp_limb digit = value & mask;
      value >>= shift;
      bits -= shift;

      if (digit < 10) *s++ = '0' + digit;
      else *s++ = 'a' + digit - 10;
    }
  }
  else {
    /* Check for overflow in size calculation */
    if ((size_t)xlen > SIZE_MAX / sizeof(mp_limb)) {
      mrb_raise(mrb, E_RUNTIME_ERROR, "bigint size too large for string conversion");
    }
    mp_limb *t = (mp_limb*)mrb_malloc(mrb, xlen * sizeof(mp_limb));

    mp_limb *tend = t + xlen;
    memcpy(t, x->p, xlen * sizeof(mp_limb));
    mp_limb b2 = base_limit[base-3];

    for (;;) {
      mp_limb *d = tend;
      mp_dbl_limb a = 0;
      while (--d >= t) {
        mp_limb d0 = *d;
        a = (a<<DIG_SIZE) | d0;
        *d = (mp_limb)(a / b2);
        a %= b2;
      }

      // convert to character
      for (mp_limb b=b2; b>=base; b/=(mp_limb)base) {
        char a0 = (char)(a % base);
        if (a0 < 10) a0 += '0';
        else a0 += 'a' - 10;
        if (s == se) break;
        *s++ = a0;
        a /= base;
      }

      // check if number is zero
      for (d = t; d < tend; d++) {
        if (*d != 0) break;
      }
      if (d == tend) break;
    }
    mrb_free(mrb, t);
  }

  while (ps<s && s[-1]=='0') s--;
  if (x->sn < 0) {
    *s++ = '-';
  }

  /* reverse string */
  for (char *u = ps,*v=s-1; u < v; u++,v--) {
    char temp = *u;
    *u = *v;
    *v = temp;
  }
  *s = '\0'; /* null termination */
  return ps;
}

static int
mpz_get_int(mpz_t *y, mrb_int *v)
{
  if (zero_p(y)) {
    *v = 0;
    return TRUE;
  }

#ifdef MRB_NO_MPZ64BIT
  /* When using 16-bit limbs, we need to handle larger accumulation */
  mrb_uint i = 0;
  mp_limb *d = y->p + y->sz;

  while (d-- > y->p) {
    /* Check for overflow before shifting */
    if (i > (mrb_uint)(MRB_INT_MAX >> DIG_SIZE)) {
      return FALSE;
    }
    i = (i << DIG_SIZE) | *d;
  }

  if (i > (mrb_uint)MRB_INT_MAX) {
    return FALSE;
  }
#else
  /* Original logic for 32-bit limbs */
  mp_dbl_limb i = 0;
  mp_limb *d = y->p + y->sz;

  while (d-- > y->p) {
    if (HIGH(i) != 0) {
      /* will overflow */
      return FALSE;
    }
    i = (i << DIG_SIZE) | *d;
  }
  if (i > MRB_INT_MAX) {
    /* overflow */
    return FALSE;
  }
#endif

  if (y->sn < 0) {
    *v = -(mrb_int)i;
  }
  else {
    *v = (mrb_int)i;
  }
  return TRUE;
}

/* Maximum bits for bigint operations (128MB of limb data) */
#define MRB_BIGINT_BIT_LIMIT (128 * 1024 * 1024 * (size_t)8)

static void
mpz_mul_2exp(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mrb_int e)
{
  if (e==0) {
    mpz_set(ctx, z, x);
    trim(z);
  }
  else {
    short sn = x->sn;
    size_t digs = e / DIG_SIZE;
    size_t bs = e % DIG_SIZE;
    mpz_t y;

    /* Check for result size overflow before allocation */
    size_t result_bits = (size_t)x->sz * DIG_SIZE + (size_t)e;
    if ((size_t)e > MRB_BIGINT_BIT_LIMIT || result_bits > MRB_BIGINT_BIT_LIMIT) {
      mrb_state *mrb = MPZ_MRB(ctx);
      mrb_raise(mrb, E_RANGE_ERROR, "shift width too large");
    }
    mpz_init_heap(ctx, &y, x->sz+digs);
    for (size_t i=0;i<x->sz;i++)
      y.p[i+digs] = x->p[i];
    if (bs) {
      ulshift(ctx, z, &y, bs);
      mpz_clear(ctx, &y);
    }
    else {
      mpz_move(ctx, z, &y);
      trim(z);
    }
    if (uzero_p(z))
      z->sn = 0;
    else
      z->sn = sn;
  }
}

static void
mpz_div_2exp(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mrb_int e)
{
  short sn = x->sn;
  if (e == 0) {
    if (z != x) {
      mpz_clear(ctx, z);
      mpz_init_heap(ctx, z, x->sz);
      mpz_set(ctx, z, x);
      trim(z);
    }
    /* else: z == x, nothing to do */
  }
  else {
    size_t digs = e / DIG_SIZE;
    size_t bs = e % DIG_SIZE;

    /* If shifting by more limbs than we have, result is zero */
    if (digs >= x->sz) {
      zero(z);
      return;
    }

    size_t new_size = x->sz - digs;

    /* In-place optimization: when z == x, no allocation needed */
    if (z == x) {
      /* Shift by whole limbs: memmove in place */
      if (digs > 0) {
        memmove(z->p, z->p + digs, new_size * sizeof(mp_limb));
      }
      z->sz = new_size;
      /* Shift by remaining bits: mpn_rshift supports in-place */
      if (bs) {
        mpn_rshift(z->p, z->p, new_size, (unsigned int)bs);
      }
      trim(z);
      if (uzero_p(z))
        z->sn = 0;
      else
        z->sn = sn;
      return;
    }

    /* General case: z != x, use temporary */
    mpz_t y;
    mpz_init_temp(ctx, &y, new_size);
    mpz_realloc(ctx, &y, new_size);
    for (size_t i = 0; i < new_size; i++)
      y.p[i] = x->p[i + digs];
    if (bs) {
      mpz_clear(ctx, z);
      mpz_init_heap(ctx, z, new_size);
      urshift(ctx, z, &y, bs);
      mpz_clear(ctx, &y);
    }
    else {
      mpz_move(ctx, z, &y);
      trim(z);
    }
    if (uzero_p(z))
      z->sn = 0;
    else {
      z->sn = sn;
    }
  }
}

static void
mpz_neg(mpz_ctx_t *ctx, mpz_t *x, mpz_t *y)
{
  /* In-place optimization: just flip the sign */
  if (x == y) {
    x->sn = -(y->sn);
    return;
  }
  mpz_init_heap(ctx, x, y->sz);
  mpz_set(ctx, x, y);
  trim(x);
  x->sn = -(y->sn);
}

#define make_2comp(v,c) do { v=~(v)+(c); c=((v)==0 && (c));} while (0)

static void
mpz_and(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mpz_t *y)
{
  if (zero_p(x) || zero_p(y)) {
    mpz_init(ctx, z);
    zero(z);
    return;
  }
  mrb_assert(x->sz > 0 || y->sz > 0);

  /* Fast path: both positive - just AND the limbs */
  if (x->sn > 0 && y->sn > 0) {
    size_t min_sz = (x->sz < y->sz) ? x->sz : y->sz;
    mpz_init_heap(ctx, z, min_sz);
    for (size_t i = 0; i < min_sz; i++) {
      z->p[i] = x->p[i] & y->p[i];
    }
    z->sn = 1;
    trim(z);
    return;
  }

  /* Slow path: at least one negative operand */
  size_t max_sz = (x->sz > y->sz) ? x->sz : y->sz;
  mpz_init_heap(ctx, z, max_sz);
  z->sn = (x->sn < 0 && y->sn < 0) ? -1 : 1;

  char c1 = 1, c2 = 1, c3 = 1;
  for (size_t i = 0; i < max_sz; i++) {
    mp_limb xv = (i < x->sz) ? x->p[i] : 0;
    mp_limb yv = (i < y->sz) ? y->p[i] : 0;

    if (x->sn < 0) make_2comp(xv, c1);
    if (y->sn < 0) make_2comp(yv, c2);
    mp_limb zv = xv & yv;
    if (z->sn < 0) make_2comp(zv, c3);
    z->p[i] = zv;
  }
  trim(z);
}

static void
mpz_or(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mpz_t *y)
{
  if (zero_p(x)) {
    mpz_init_heap(ctx, z, y->sz);
    mpz_set(ctx, z, y);
    trim(z);
    return;
  }
  if (zero_p(y)) {
    mpz_init_heap(ctx, z, x->sz);
    mpz_set(ctx, z, x);
    trim(z);
    return;
  }
  mrb_assert(x->sz > 0 || y->sz > 0);

  size_t max_sz = (x->sz > y->sz) ? x->sz : y->sz;

  /* Fast path: both positive - just OR the limbs */
  if (x->sn > 0 && y->sn > 0) {
    mpz_init_heap(ctx, z, max_sz);
    size_t min_sz = (x->sz < y->sz) ? x->sz : y->sz;
    for (size_t i = 0; i < min_sz; i++) {
      z->p[i] = x->p[i] | y->p[i];
    }
    /* Copy remaining limbs from the longer operand */
    if (x->sz > y->sz) {
      for (size_t i = min_sz; i < max_sz; i++) z->p[i] = x->p[i];
    }
    else {
      for (size_t i = min_sz; i < max_sz; i++) z->p[i] = y->p[i];
    }
    z->sn = 1;
    trim(z);
    return;
  }

  /* Slow path: at least one negative operand */
  mpz_init_heap(ctx, z, max_sz);
  z->sn = (x->sn == y->sn) ? x->sn : -1;

  char c1 = 1, c2 = 1, c3 = 1;
  for (size_t i = 0; i < max_sz; i++) {
    mp_limb xv = (i < x->sz) ? x->p[i] : 0;
    mp_limb yv = (i < y->sz) ? y->p[i] : 0;

    if (x->sn < 0) make_2comp(xv, c1);
    if (y->sn < 0) make_2comp(yv, c2);
    mp_limb zv = xv | yv;
    if (z->sn < 0) make_2comp(zv, c3);
    z->p[i] = zv;
  }
  trim(z);
}

static void
mpz_xor(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x, mpz_t *y)
{
  if (zero_p(x)) {
    mpz_init_heap(ctx, z, y->sz);
    mpz_set(ctx, z, y);
    trim(z);
    return;
  }
  if (zero_p(y)) {
    mpz_init_heap(ctx, z, x->sz);
    mpz_set(ctx, z, x);
    trim(z);
    return;
  }
  mrb_assert(x->sz > 0 || y->sz > 0);

  size_t max_sz = (x->sz > y->sz) ? x->sz : y->sz;

  /* Fast path: both positive - just XOR the limbs */
  if (x->sn > 0 && y->sn > 0) {
    mpz_init_heap(ctx, z, max_sz);
    size_t min_sz = (x->sz < y->sz) ? x->sz : y->sz;
    for (size_t i = 0; i < min_sz; i++) {
      z->p[i] = x->p[i] ^ y->p[i];
    }
    /* Copy remaining limbs from the longer operand (XOR with 0) */
    if (x->sz > y->sz) {
      for (size_t i = min_sz; i < max_sz; i++) z->p[i] = x->p[i];
    }
    else {
      for (size_t i = min_sz; i < max_sz; i++) z->p[i] = y->p[i];
    }
    z->sn = 1;
    trim(z);
    return;
  }

  /* Slow path: at least one negative operand */
  mpz_init_heap(ctx, z, max_sz);
  z->sn = (x->sn == y->sn) ? 1 : -1;

  char c1 = 1, c2 = 1, c3 = 1;
  for (size_t i = 0; i < max_sz; i++) {
    mp_limb xv = (i < x->sz) ? x->p[i] : 0;
    mp_limb yv = (i < y->sz) ? y->p[i] : 0;

    if (x->sn < 0) make_2comp(xv, c1);
    if (y->sn < 0) make_2comp(yv, c2);
    mp_limb zv = xv ^ yv;
    if (z->sn < 0) make_2comp(zv, c3);
    z->p[i] = zv;
  }
  trim(z);
}

static void
mpz_pow(mpz_ctx_t *ctx, mpz_t *zz, mpz_t *x, mrb_int e)
{
  if (e == 0) {
    mpz_init_set_int(ctx, zz, 1L);
    return;
  }

  mrb_uint mask = 1ULL << (sizeof(mrb_int) * 8 - 1);
  while (mask != 0 && !(mask & e)) {
    mask >>= 1;
  }

  /* Set initial value to x for exponentiation */
  mpz_init_set(ctx, zz, x);

  if (mask == 0) { /* e is 0 or 1 */
    if (e == 0) mpz_set_int(ctx, zz, 1L);
    return;
  }

  mask >>= 1;

  /* Pre-allocate a single temporary variable */
  mpz_t temp;
  mpz_init(ctx, &temp);

  for (; mask != 0; mask >>= 1) {
    /* squaring: temp = zz * zz */
    mpz_mul(ctx, &temp, zz, zz);

    if (e & mask) {
      /* multiplication: zz = temp * x */
      mpz_mul(ctx, zz, &temp, x);
    }
    else {
      /* move result: zz = temp */
      mpz_move(ctx, zz, &temp);
    }
  }
  mpz_clear(ctx, &temp);
}

static void
mpz_powm(mpz_ctx_t *ctx, mpz_t *zz, mpz_t *x, mpz_t *ex, mpz_t *n)
{
  /* Handle special cases */
  if (zero_p(ex) || uzero_p(ex)) {
    mpz_set_int(ctx, zz, 1);
    return;
  }

  if (ex->sn < 0) {
    return;
  }

  /* Check modulus size before allocating large temporaries. */
  {
    size_t mod_bits = (size_t)n->sz * DIG_SIZE;
    if (mod_bits > MRB_BIGINT_BIT_LIMIT / 2) {
      mrb_state *mrb = MPZ_MRB(ctx);
      mrb_raise(mrb, E_RANGE_ERROR, "modulus too large");
    }
  }

  mpz_t t, b;
  mpz_init_set_int(ctx, &t, 1);
  mpz_init_set(ctx, &b, x);

  mpz_t temp;
  mpz_init_temp(ctx, &temp, n->sz * 2);  /* For intermediate calculations */

  size_t len = digits(ex);
  for (size_t i = 0; i < len; i++) {
    mp_limb e = ex->p[i];
    for (size_t j = 0; j < sizeof(mp_limb) * 8; j++) {
      if ((e & 1) == 1) {
        mpz_mul(ctx, &temp, &t, &b);
        mpz_mod(ctx, &t, &temp, n);
      }
      e >>= 1;
      mpz_mul(ctx, &temp, &b, &b);
      mpz_mod(ctx, &b, &temp, n);
    }
  }

  mpz_move(ctx, zz, &t);
  mpz_clear(ctx, &t);
  mpz_clear(ctx, &b);
  mpz_clear(ctx, &temp);
}

static void
mpz_powm_i(mpz_ctx_t *ctx, mpz_t *zz, mpz_t *x, mrb_int ex, mpz_t *n)
{
  if (ex == 0) {
    mpz_set_int(ctx, zz, 1);
    return;
  }

  if (ex < 0) {
    return;
  }

  mpz_t t, b;
  mpz_init_set_int(ctx, &t, 1);
  mpz_init_set(ctx, &b, x);
  mpz_t temp;

  mpz_init_temp(ctx, &temp, n->sz * 2);  /* For intermediate calculations */
  while (ex > 0) {
    if ((ex & 1) == 1) {
      mpz_mul(ctx, &temp, &t, &b);
      mpz_mod(ctx, &t, &temp, n);
    }
    ex >>= 1;
    if (ex > 0) {  /* Skip final squaring when ex becomes 0 */
      mpz_mul(ctx, &temp, &b, &b);
      mpz_mod(ctx, &b, &temp, n);
    }
  }

  mpz_move(ctx, zz, &t);
  mpz_clear(ctx, &t);
  mpz_clear(ctx, &b);
  mpz_clear(ctx, &temp);
}

/* Helper functions for pool-based GCD operations */
static int
mpz_abs_copy(mpz_ctx_t *ctx, mpz_t *result, mpz_t *operand) {
  if (!operand || operand->sz == 0) {
    result->sz = 0;
    result->sn = 0;
    return 1;
  }

  /* Copy limbs */
  for (size_t i = 0; i < operand->sz && i < result->sz; i++) {
    result->p[i] = operand->p[i];
  }
  result->sz = (operand->sz < result->sz) ? operand->sz : result->sz;
  result->sn = (operand->sn < 0) ? -operand->sn : operand->sn; /* Always positive */

  return 1;
}

static void
mpz_abs(mpz_ctx_t *ctx, mpz_t *x, mpz_t *y)
{
  /* In-place optimization: just make sign positive */
  if (x == y) {
    if (x->sn < 0) x->sn = -x->sn;
    return;
  }
  mpz_init_heap(ctx, x, y->sz);
  mpz_abs_copy(ctx, x, y);
}

/* Fast GCD for single limbs using binary algorithm */
static mp_limb
limb_gcd(mp_limb a, mp_limb b)
{
  if (a == 0) return b;
  if (b == 0) return a;

  /* Find power of 2 dividing both a and b */
  int shift = 0;
  while (((a | b) & 1) == 0) {
    a >>= 1;
    b >>= 1;
    shift++;
  }

  /* Make a odd */
  while ((a & 1) == 0) {
    a >>= 1;
  }

  /* From here on, a is always odd */
  do {
    /* Make b odd */
    while ((b & 1) == 0) {
      b >>= 1;
    }

    /* Now both a and b are odd. Ensure a >= b */
    if (a < b) {
      mp_limb temp = a;
      a = b;
      b = temp;
    }

    /* Replace b with (b - a) */
    b = b - a;

  } while (b != 0);

  /* Restore common factors of 2 */
  return a << shift;
}

/* Count trailing zero bits in a multi-precision integer */
static size_t
mpz_trailing_zeros(mpz_t *x)
{
  if (zero_p(x) || x->sz == 0) return 0;

  size_t zeros = 0;

  /* Count complete zero limbs */
  size_t i = 0;
  while (i < x->sz && x->p[i] == 0) {
    zeros += DIG_SIZE;
    i++;
  }

  /* Count trailing zeros in first non-zero limb */
  if (i < x->sz) {
    mp_limb limb = x->p[i];
#if (defined(__GNUC__) || __has_builtin(__builtin_ctzll))
    if (sizeof(mp_limb) == sizeof(unsigned long long)) {
      zeros += __builtin_ctzll(limb);
    }
    else if (sizeof(mp_limb) == sizeof(unsigned long)) {
      zeros += __builtin_ctzl(limb);
    }
    else {
      zeros += __builtin_ctz(limb);
    }
#else
    /* Fallback bit counting */
    while ((limb & 1) == 0) {
      limb >>= 1;
      zeros++;
    }
#endif
  }

  return zeros;
}

/* Check if a number is a power of 2 */
static int
mpz_power_of_2_p(mpz_t *x)
{
  if (zero_p(x) || x->sz == 0) return 0;

  /* Count non-zero limbs */
  size_t non_zero_limbs = 0;
  size_t non_zero_index = 0;

  for (size_t i = 0; i < x->sz; i++) {
    if (x->p[i] != 0) {
      non_zero_limbs++;
      non_zero_index = i;
      if (non_zero_limbs > 1) return 0; /* More than one non-zero limb */
    }
  }

  if (non_zero_limbs == 0) return 0; /* All zero */
  if (non_zero_limbs > 1) return 0;  /* Multiple non-zero limbs */

  /* Check if the single non-zero limb is a power of 2 */
  mp_limb limb = x->p[non_zero_index];
  return (limb != 0) && ((limb & (limb - 1)) == 0);
}

/* Binary GCD algorithm (Stein's algorithm) - faster than Euclidean GCD */
static void
mpz_gcd(mpz_ctx_t *ctx, mpz_t *gg, mpz_t *aa, mpz_t *bb)
{
  mpz_t a, b;
  size_t shift;
  size_t a_zeros;
  size_t b_zeros;

  /* Handle special cases */
  if (zero_p(aa)) {
    mpz_abs(ctx, gg, bb);
    return;
  }
  if (zero_p(bb)) {
    mpz_abs(ctx, gg, aa);
    return;
  }

  /* Fast path for single-limb numbers */
  if (aa->sz <= 1 && bb->sz <= 1) {
    mp_limb a_limb = (aa->sz == 0) ? 0 : aa->p[0];
    mp_limb b_limb = (bb->sz == 0) ? 0 : bb->p[0];
    mp_limb result = limb_gcd(a_limb, b_limb);

    mpz_init(ctx, gg);
    if (result == 0) {
      gg->sn = 0;
      gg->sz = 0;
    }
    else {
      mpz_realloc(ctx, gg, 1);
      gg->p[0] = result;
      gg->sn = 1;
    }
    return;
  }

  /* Fast path for powers of 2 */
  if (mpz_power_of_2_p(aa)) {
    a_zeros = mpz_trailing_zeros(aa);
    b_zeros = mpz_trailing_zeros(bb);
    size_t min_zeros = (a_zeros < b_zeros) ? a_zeros : b_zeros;

    mpz_init_set_int(ctx, gg, 1);
    mpz_mul_2exp(ctx, gg, gg, min_zeros);
    return;
  }
  if (mpz_power_of_2_p(bb)) {
    a_zeros = mpz_trailing_zeros(aa);
    b_zeros = mpz_trailing_zeros(bb);
    size_t min_zeros = (a_zeros < b_zeros) ? a_zeros : b_zeros;

    mpz_init_set_int(ctx, gg, 1);
    mpz_mul_2exp(ctx, gg, gg, min_zeros);
    return;
  }

  mpz_init(ctx, &a);
  mpz_abs(ctx, &a, aa);
  mpz_init(ctx, &b);
  mpz_abs(ctx, &b, bb);

  shift = 0;
  a_zeros = mpz_trailing_zeros(&a);
  b_zeros = mpz_trailing_zeros(&b);
  shift = (a_zeros < b_zeros) ? a_zeros : b_zeros;

  mpz_div_2exp(ctx, &a, &a, a_zeros);
  mpz_div_2exp(ctx, &b, &b, b_zeros);

  /* Euclidean algorithm for multi-limb numbers */
  while (!zero_p(&b)) {
    mpz_t temp;
    mpz_init_temp(ctx, &temp, a.sz);
    mpz_mod(ctx, &temp, &a, &b);
    mpz_move(ctx, &a, &b);
    mpz_move(ctx, &b, &temp);
    mpz_clear(ctx, &temp);
  }
  mpz_mul_2exp(ctx, gg, &a, shift);
  mpz_clear(ctx, &a);
  mpz_clear(ctx, &b);
}


static size_t
mpz_bits(const mpz_t *x)
{
  if (x->sz == 0 || x->sn == 0) return 0;

  size_t limb_bits = sizeof(mp_limb) * 8;

  // Get the most significant limb
  size_t i = x->sz - 1;
  mp_limb high = x->p[i];

  // Number of bits = total full limbs + significant bits in top limb
  return i * limb_bits + (limb_bits - lzb(high));
}

static void
mpz_sqrt(mpz_ctx_t *ctx, mpz_t *z, mpz_t *x)
{
  mrb_assert(x->sn >= 0);

  if (x->sz == 0) {
    // sqrt(0) = 0
    mpz_init(ctx, z);
    z->sn = 0;
    z->sz = 0;
    return;
  }

  // Use heap-only implementation for now
  size_t xbits = mpz_bits(x);
  size_t sbit = (xbits + 1) / 2;
  mpz_t s, t;
  mpz_init_set_int(ctx, &s, 1);
  mpz_mul_2exp(ctx, &s, &s, sbit);

  mpz_init_temp(ctx, &t, x->sz + 1);

  // Iteratively refine s using Newton-Raphson method:
  // s = (s + x / s) / 2
  for (;;) {
    mpz_mdiv(ctx, &t, x, &s);     // t = x / s
    mpz_add(ctx, &t, &t, &s);     // t = s + x/s
    mpz_div_2exp(ctx, &t, &t, 1); // t = (s + x/s) / 2

    if (mpz_cmp(ctx, &t, &s) >= 0) {
      // Converged: t >= s
      break;
    }

    mpz_set(ctx, &s, &t);
  }

  mpz_move(ctx, z, &s);
  mpz_clear(ctx, &t);
}

/* --- mruby functions --- */
/* initialize mpz_t from RBigint (not need to clear) */
static void
bint_as_mpz(struct RBigint *b, mpz_t *x)
{
  x->p = RBIGINT_ARY(b);
  x->sz = RBIGINT_SIZE(b);
  x->sn = RBIGINT_SIGN(b);
}

/* Transfer mpz_t data to RBigint structure */
static void
bint_set(mpz_ctx_t *ctx, struct RBigint *b, mpz_t *x)
{
  if (x->sz <= RBIGINT_EMBED_SIZE_MAX) {
    RBIGINT_SET_EMBED_SIZE(b, x->sz);
    RBIGINT_SET_EMBED_SIGN(b, x->sn);
    if (x->p) {
      memcpy(RBIGINT_EMBED_ARY(b), x->p, x->sz*sizeof(mp_limb));
    }
    else {
      /* Initialize embedded array to zero when x->p is NULL */
      memset(RBIGINT_EMBED_ARY(b), 0, x->sz*sizeof(mp_limb));
    }
    mpz_clear(ctx, x);
  }
  else {
    RBIGINT_SET_HEAP(b);
    mpz_init(ctx, &b->as.heap);  /* Initialize before mpz_move */
    mpz_move(ctx, &b->as.heap, x);
  }
}

static struct RBigint*
bint_new(mpz_ctx_t *ctx, mpz_t *x)
{
  struct RBigint *b = MRB_OBJ_ALLOC(MPZ_MRB(ctx), MRB_TT_BIGINT, MPZ_MRB(ctx)->integer_class);
  bint_set(ctx, b, x);
  return b;
}

static struct RBigint*
bint_new_int(mpz_ctx_t *ctx, mrb_int n)
{
  mpz_t x;
  mpz_init_set_int(ctx, &x, n);
  return bint_new(ctx, &x);
}

mrb_value
mrb_bint_new_int(mrb_state *mrb, mrb_int x)
{
  MPZ_CTX_INIT(mrb, ctx, pool);
  struct RBigint *b = bint_new_int(ctx, x);
  return mrb_obj_value(b);
}

#ifdef MRB_INT32
mrb_value
mrb_bint_new_int64(mrb_state *mrb, int64_t n)
{
  mpz_t x;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_set_int64(ctx, &x, n);
  struct RBigint *b = bint_new(ctx, &x);
  return mrb_obj_value(b);
}
#endif

mrb_value
mrb_bint_new_uint64(mrb_state *mrb, uint64_t x)
{
  mpz_t z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &z);
  mpz_set_uint64(ctx, &z, x);
  struct RBigint *b = bint_new(ctx, &z);
  return mrb_obj_value(b);
}

static mrb_value
bint_norm(mrb_state *mrb, struct RBigint *b)
{
  mrb_int i;
  mpz_t a;

  bint_as_mpz(b, &a);
  if (mpz_get_int(&a, &i)) {
    return mrb_int_value(mrb, i);
  }
  return mrb_obj_value(b);
}

mrb_value
mrb_bint_new_str(mrb_state *mrb, const char *x, mrb_int len, mrb_int base)
{
  mpz_t z;
  int sn = 1;

  if (base < 0) {
    base = -base;
    sn = -1;
  }
  mrb_assert(2 <= base && base <= 36);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init_set_str(ctx, &z, x, len, base);
  if (sn < 0) {
    z.sn = sn;
  }
  return bint_norm(mrb, bint_new(ctx, &z));
}

void
mrb_gc_free_bint(mrb_state *mrb, struct RBasic *x)
{
  struct RBigint *b = (struct RBigint*)x;
  MPZ_CTX_INIT(mrb, ctx, pool);

  if (!RBIGINT_EMBED_P(b)) {
    mpz_clear(ctx, &b->as.heap);
  }
}

#ifndef MRB_NO_FLOAT
mrb_value
mrb_bint_new_float(mrb_state *mrb, mrb_float x)
{
  /* x should not be NaN nor Infinity */
  mrb_assert(x == x && x != x * 0.5);

  if (FIXABLE_FLOAT(x)) {
    return mrb_int_value(mrb, (mrb_int)x);
  }

  int sn;
  if (x < 0.0) {
    x = -x;
    sn = -1;
  }
  else {
    sn = 1;
  }
  if (x < 1.0) {
    return mrb_fixnum_value(0);
  }

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_t r;
  mpz_init(ctx, &r);
  r.sn = sn;

  mrb_float b = (double)DIG_BASE;
  mrb_float bi = 1.0 / b;
  size_t rn;

  for (rn = 1; x >= b; rn++)
    x *= bi;

  mpz_realloc(ctx, &r, rn);
  mp_limb *rp = r.p;
  for (size_t i=rn-1;;i--) {
    mp_limb f = LOW((mp_limb)x);
    x -= f;
    mrb_assert(x < 1.0);
    rp[i] = f;
    if (i == 0) break;
  }
  return bint_norm(mrb, bint_new(ctx, &r));
}

mrb_float
mrb_bint_as_float(mrb_state *mrb, mrb_value self)
{
  mpz_t m;
  bint_as_mpz(RBIGINT(self), &m);

  mp_limb *d = m.p + m.sz;
  mrb_float val = 0;

  while (d-- > m.p) {
    val = val * DIG_BASE + *d;
  }

  if (m.sn < 0) {
    val = -val;
  }
  return val;
}
#endif

mrb_value
mrb_as_bint(mrb_state *mrb, mrb_value x)
{
  if (mrb_bigint_p(x)) return x;
  return mrb_bint_new_int(mrb, mrb_as_int(mrb, x));
}

mrb_int
mrb_bint_as_int(mrb_state *mrb, mrb_value x)
{
  mpz_t m;
  mrb_int i;

  bint_as_mpz(RBIGINT(x), &m);
  if (!mpz_get_int(&m, &i)) {
    mrb_raise(mrb, E_RANGE_ERROR, "integer out of range");
  }
  return i;
}

#ifdef MRB_INT32
int64_t
mrb_bint_as_int64(mrb_state *mrb, mrb_value x)
{
  mpz_t m;
  bint_as_mpz(RBIGINT(x), &m);

  uint64_t u = 0;
  size_t len = digits(&m);

  if (len*sizeof(mp_limb) > sizeof(uint64_t)) {
  out_of_range:
    mrb_raise(mrb, E_RANGE_ERROR, "integer out of range");
  }
  for (size_t i=len-1; ; i--) {
    u <<= DIG_SIZE;
    u |= m.p[i];
    if (i==0) break;
  }
  if (u > INT64_MAX) goto out_of_range;
  if (m.sn < 0) return -(int64_t)u;
  return (int64_t)u;
}
#endif

uint64_t
mrb_bint_as_uint64(mrb_state *mrb, mrb_value x)
{
  mpz_t m;
  bint_as_mpz(RBIGINT(x), &m);

  uint64_t u = 0;
  size_t len = digits(&m);

  if (m.sn < 0 || len*sizeof(mp_limb) > sizeof(uint64_t)) {
    mrb_raise(mrb, E_RANGE_ERROR, "integer out of range");
  }
  for (size_t i=len-1; ; i--) {
    u <<= DIG_SIZE;
    u |= m.p[i];
    if (i==0) break;
  }
  return u;
}

static mrb_bool
int_fit_limb_p(mrb_int i)
{
#if DIG_SIZE == 32
# ifdef MRB_INT64
  // if mp_limb is int32_t
  return (i > INT32_MIN && i <= INT32_MAX);
# else
  // if mp_limb is also int32_t, it always fits
  return TRUE;
# endif
#else /* if DIG_SIZE == 16 */
  // if mp_limb is int16_t
  return (i > INT16_MIN && i <= INT16_MAX);
#endif
}

/* unnormalize version of mrb_bint_add */
mrb_value
mrb_bint_add_n(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b, z;

  bint_as_mpz(RBIGINT(x), &a);

  MPZ_CTX_INIT(mrb, ctx, pool);

  if (mrb_integer_p(y)) {
    mrb_int n = mrb_integer(y);
    if (int_fit_limb_p(n)) {
      mpz_init_set(ctx, &z, &a);
      if ((n > 0) ^ (z.sn > 0)) {
        mpz_sub_int(ctx, &z, n<0 ? -n : n);
      }
      else {
        mpz_add_int(ctx, &z, n<0 ? -n : n);
      }
      struct RBigint *v = bint_new(ctx, &z);
      return mrb_obj_value(v);
    }
  }
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  mpz_init(ctx, &z);
  mpz_add(ctx, &z, &a, &b);
  struct RBigint *v = bint_new(ctx, &z);
  return mrb_obj_value(v);
}

mrb_value
mrb_bint_add(mrb_state *mrb, mrb_value x, mrb_value y)
{
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    return mrb_float_value(mrb,v1+v2);
  }
#endif
  x = mrb_bint_add_n(mrb, x, y);
  return bint_norm(mrb, RBIGINT(x));
}

/* unnormalize version of mrb_bint_sub */
mrb_value
mrb_bint_sub_n(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  if (mrb_integer_p(y)) {
    mrb_int n = mrb_integer(y);
    if (int_fit_limb_p(n)) {
      mpz_init_set(ctx, &z, &a);
      if ((n > 0) ^ (z.sn > 0)) {
        mpz_add_int(ctx, &z, n<0 ? -n : n);
      }
      else {
        mpz_sub_int(ctx, &z, n<0 ? -n : n);
      }
      struct RBigint *v = bint_new(ctx, &z);
      return mrb_obj_value(v);
    }
  }
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  mpz_init(ctx, &z);
  mpz_sub(ctx, &z, &a, &b);
  struct RBigint *v = bint_new(ctx, &z);
  return mrb_obj_value(v);
}

mrb_value
mrb_bint_sub(mrb_state *mrb, mrb_value x, mrb_value y)
{
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    return mrb_float_value(mrb,v1-v2);
  }
#endif
  x = mrb_bint_sub_n(mrb, x, y);
  return bint_norm(mrb, RBIGINT(x));
}

struct bint_mul_data {
  mpz_ctx_t *ctx;
  mpz_t *a;
  mpz_t *b;
  mpz_t z;  /* cleanup target */
};

static mrb_value
bint_mul_body(mrb_state *mrb, void *userdata)
{
  struct bint_mul_data *d = (struct bint_mul_data *)userdata;
  mpz_init(d->ctx, &d->z);
  mpz_mul(d->ctx, &d->z, d->a, d->b);
  return mrb_nil_value();
}

static struct RBigint*
bint_mul(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b;

  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(x), &a);
  bint_as_mpz(RBIGINT(y), &b);

  MPZ_CTX_INIT(mrb, ctx, pool);
  struct bint_mul_data d = {ctx, &a, &b, {0,0,0}};
  mrb_value exc;
  MRB_ENSURE(mrb, exc, bint_mul_body, &d) {
    /* On exception, cleanup z (mpz_clear is safe on zero-initialized mpz_t) */
    if (mrb->exc) {
      mpz_clear(ctx, &d.z);
    }
  }
  return bint_new(ctx, &d.z);
}

mrb_value
mrb_bint_mul(mrb_state *mrb, mrb_value x, mrb_value y)
{
  if (mrb_integer_p(y)) {
    if (mrb_integer(y) == 0) return mrb_fixnum_value(0);
    if (mrb_integer(y) == 1) return bint_norm(mrb, RBIGINT(x));
  }
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    return mrb_float_value(mrb,v1*v2);
  }
#endif
  return bint_norm(mrb, bint_mul(mrb, x, y));
}

mrb_value
mrb_bint_mul_n(mrb_state *mrb, mrb_value x, mrb_value y)
{
  struct RBigint *b = bint_mul(mrb, x, y);
  return mrb_obj_value(b);
}

mrb_value
mrb_bint_div(mrb_state *mrb, mrb_value x, mrb_value y)
{
  if (mrb_integer_p(y)) {
    if (mrb_integer(y) == 0) mrb_int_zerodiv(mrb);
    if (mrb_integer(y) == 1) return bint_norm(mrb, RBIGINT(x));
  }
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    return mrb_float_value(mrb,v1*v2);
  }
#endif
  mpz_t a, b, z;

  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&b) || uzero_p(&b)) {
    mrb_int_zerodiv(mrb);
  }
  bint_as_mpz(RBIGINT(x), &a);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init(ctx, &z);
  mpz_mdiv(ctx, &z, &a, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_add_ii(mrb_state *mrb, mrb_int x, mrb_int y)
{
  mpz_t a, b, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &z);
  mpz_init_set_int(ctx, &a, x);
  mpz_init_set_int(ctx, &b, y);
  mpz_add(ctx, &z, &a, &b);
  mpz_clear(ctx, &a);
  mpz_clear(ctx, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_sub_ii(mrb_state *mrb, mrb_int x, mrb_int y)
{
  mpz_t a, b, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &z);
  mpz_init_set_int(ctx, &a, x);
  mpz_init_set_int(ctx, &b, y);
  mpz_sub(ctx, &z, &a, &b);
  mpz_clear(ctx, &a);
  mpz_clear(ctx, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_mul_ii(mrb_state *mrb, mrb_int x, mrb_int y)
{
  mpz_t a, b, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &z);
  mpz_init_set_int(ctx, &a, x);
  mpz_init_set_int(ctx, &b, y);
  mpz_mul(ctx, &z, &a, &b);
  mpz_clear(ctx, &a);
  mpz_clear(ctx, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_mod(mrb_state *mrb, mrb_value x, mrb_value y)
{
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    return mrb_float_value(mrb, fmod(v1, v2));
  }
#endif
  if (mrb_integer_p(y) && mrb_integer(y) == 0) {
    mrb_int_zerodiv(mrb);
  }
  mpz_t a, b, z;
  x = mrb_as_bint(mrb, x);
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&b) || uzero_p(&b)) {
    mrb_int_zerodiv(mrb);
  }
  bint_as_mpz(RBIGINT(x), &a);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init(ctx, &z);
  mpz_mmod(ctx, &z, &a, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_rem(mrb_state *mrb, mrb_value x, mrb_value y)
{
  /* called from mrbgems/mruby-numeric-ext/src/numeric_ext.c */
  /* y should not be float */
  if (mrb_integer_p(y) && mrb_integer(y) == 0) {
    mrb_int_zerodiv(mrb);
  }
  mpz_t a, b, z;
  x = mrb_as_bint(mrb, x);
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&b) || uzero_p(&b)) {
    mrb_int_zerodiv(mrb);
  }
  bint_as_mpz(RBIGINT(x), &a);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init(ctx, &z);
  mpz_mod(ctx, &z, &a, &b);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_divmod(mrb_state *mrb, mrb_value x, mrb_value y)
{
  /* called from src/numeric.c */
  /* y should not be float */
  if (mrb_integer_p(y) && mrb_integer(y) == 0) {
    mrb_int_zerodiv(mrb);
  }
  y = mrb_as_bint(mrb, y);
  mpz_t a, b, c, d;
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&b) || uzero_p(&b)) {
    mrb_int_zerodiv(mrb);
  }
  bint_as_mpz(RBIGINT(x), &a);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init(ctx, &c);
  mpz_init(ctx, &d);
  mpz_mdivmod(ctx, &c, &d, &a, &b);
  return mrb_assoc_new(mrb, bint_norm(mrb, bint_new(ctx, &c)), bint_norm(mrb, bint_new(ctx, &d)));
}

mrb_int
mrb_bint_cmp(mrb_state *mrb, mrb_value x, mrb_value y)
{
#ifndef MRB_NO_FLOAT
  if (mrb_float_p(y)) {
    mrb_float v1 = mrb_bint_as_float(mrb, x);
    mrb_float v2 = mrb_float(y);
    if (v1 == v2) return 0;
    if (v1 > v2)  return 1;
    return -1;
  }
#endif
  mpz_t a;

  bint_as_mpz(RBIGINT(x), &a);
  if (!mrb_bigint_p(y)) {
    if (!mrb_integer_p(y)) return -2; /* type mismatch */

    mrb_int i1, i2 = mrb_integer(y);
    if (mpz_get_int(&a, &i1)) {
      if (i1 == i2) return 0;
      if (i1 > i2) return 1;
      return -1;
    }
    if (a.sn > 0) return 1;
    return -1;
  }
  mpz_t b;
  bint_as_mpz(RBIGINT(y), &b);
  MPZ_CTX_INIT(mrb, ctx, pool);
  return mpz_cmp(ctx, &a, &b);
}

/* Maximum bits for power result to prevent resource exhaustion */
/* 1 million bits = ~125KB per bigint, ~300,000 decimal digits */
#define MRB_BIGINT_POW_MAX_BITS 1000000

mrb_value
mrb_bint_pow(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a;

  bint_as_mpz(RBIGINT(x), &a);
  switch (mrb_type(y)) {
  case MRB_TT_INTEGER:
    break;
  case MRB_TT_BIGINT:
    mrb_raise(mrb, E_TYPE_ERROR, "too big power");
  default:
    mrb_raisef(mrb, E_TYPE_ERROR, "%Y cannot be convert to integer", y);
  }

  mrb_int exp = mrb_integer(y);
  if (exp < 0) {
    mrb_raise(mrb, E_ARGUMENT_ERROR, "negative exponent");
  }

  /* Check if result would be too large */
  /* result_bits ≈ base_bits * exp */
  size_t base_bits = mpz_bits(&a);
  if (base_bits == 0) base_bits = 1;  /* handle 0 and 1 */
  if (exp > 0 && (size_t)exp > MRB_BIGINT_POW_MAX_BITS / base_bits) {
    mrb_raise(mrb, E_RANGE_ERROR, "exponent too large");
  }

  mpz_t z;
  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_pow(ctx, &z, &a, exp);

  struct RBigint *b = bint_new(ctx, &z);
  return mrb_obj_value(b);
}

mrb_value
mrb_bint_powm(mrb_state *mrb, mrb_value x, mrb_value exp, mrb_value mod)
{
  mpz_t a, b, c, z;
  mrb_bool neg_mod = FALSE;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  if (mrb_integer_p(mod)) {
    mrb_int m = mrb_integer(mod);
    if (m == 0) mrb_int_zerodiv(mrb);
    if (m < 0) {
      neg_mod = TRUE;
      m = -m;
    }
    mpz_init_set_int(ctx, &c, m);
  }
  else {
    mod = mrb_as_bint(mrb, mod);
    bint_as_mpz(RBIGINT(mod), &c);
    if (zero_p(&c) || uzero_p(&c)) {
      mrb_int_zerodiv(mrb);
    }
    if (c.sn < 0) {
      neg_mod = TRUE;
      c.sn = 1;  /* use absolute value */
    }
  }

  /* Check for zero base case: 0^n = 0 for n > 0 */
  if (zero_p(&a) || uzero_p(&a)) {
    mrb_bool exp_positive;
    if (mrb_bigint_p(exp)) {
      bint_as_mpz(RBIGINT(exp), &b);
      exp_positive = (b.sn > 0) && !uzero_p(&b);
    }
    else {
      exp_positive = mrb_integer(exp) > 0;
    }
    if (exp_positive) {
      /* 0^n mod m = 0 for n > 0 */
      if (mrb_integer_p(mod)) mpz_clear(ctx, &c);
      return mrb_fixnum_value(0);
    }
  }

  mpz_init(ctx, &z);
  if (mrb_bigint_p(exp)) {
    bint_as_mpz(RBIGINT(exp), &b);
    if (b.sn < 0) goto raise;
    mpz_powm(ctx, &z, &a, &b, &c);
  }
  else {
    mrb_int e = mrb_integer(exp);
    if (e < 0) goto raise;
    mpz_powm_i(ctx, &z, &a, e, &c);
  }

  /* Apply signed modulo adjustment for negative modulus */
  /* Ruby: result + m for non-zero result when m is negative */
  if (neg_mod && !zero_p(&z) && !uzero_p(&z)) {
    mpz_sub(ctx, &z, &z, &c);  /* z = z - |m| = z + m (since m is negative) */
  }

  if (mrb_integer_p(mod)) mpz_clear(ctx, &c);
  return bint_norm(mrb, bint_new(ctx, &z));

 raise:
  if (mrb_integer_p(mod)) mpz_clear(ctx, &c);
  mrb_raise(mrb, E_ARGUMENT_ERROR, "int.pow(n,m): n must be positive");
  /* not reached */
  return mrb_nil_value();
}

mrb_value
mrb_bint_to_s(mrb_state *mrb, mrb_value x, mrb_int base)
{
  mpz_t a;

  bint_as_mpz(RBIGINT(x), &a);
  if (zero_p(&a) || uzero_p(&a)) {
    return mrb_str_new_lit(mrb, "0");
  }
  size_t len = mpz_sizeinbase(&a, (int)base);
  if (sizeof(size_t) >= sizeof(mrb_int) && MRB_INT_MAX-2 < len) {
    mrb_raise(mrb, E_ARGUMENT_ERROR, "too long string from Integer");
  }
  mrb_value str = mrb_str_new(mrb, NULL, len+2);
  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_get_str(ctx, RSTRING_PTR(str), len, base, &a);
  RSTR_SET_LEN(RSTRING(str), strlen(RSTRING_PTR(str)));
  return str;
}

mrb_value
mrb_bint_and(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b, c;

  bint_as_mpz(RBIGINT(x), &a);
  if (mrb_integer_p(y)) {
    mrb_int z = mrb_integer(y);
    if (z == 0) return mrb_fixnum_value(0);
    if (z > 0 && (mp_dbl_limb)z < DIG_BASE) {
      z &= a.p[0];
      return mrb_int_value(mrb, z);
    }
    if (z == -1) return x;
  }

  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&a) || zero_p(&b)) return mrb_fixnum_value(0);

  MPZ_CTX_INIT(mrb, ctx, pool);
  mpz_init(ctx, &c);
  mpz_and(ctx, &c, &a, &b);
  return bint_norm(mrb, bint_new(ctx, &c));
}

mrb_value
mrb_bint_or(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b, c;

  bint_as_mpz(RBIGINT(x), &a);
  if (mrb_integer_p(y)) {
    mrb_int z = mrb_integer(y);
    if (z == 0) return x;
    if (z == -1) return y;
  }

  MPZ_CTX_INIT(mrb, ctx, pool);
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&a)) return y;
  if (zero_p(&b)) return x;
  mpz_init(ctx, &c);
  mpz_or(ctx, &c, &b, &a);
  return bint_norm(mrb, bint_new(ctx, &c));
}

mrb_value
mrb_bint_xor(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t a, b, c;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  if (mrb_integer_p(y) && a.sn > 0) {
    mrb_int z = mrb_integer(y);
    if (z == 0) return x;
    if (0 < z && (mp_dbl_limb)z < DIG_BASE) {
      mpz_init_set(ctx, &c, &a);
      if (a.sz == 0) {
        mpz_realloc(ctx, &c, 1);
        c.p[0] = (mp_limb)z;
      }
      else {
        c.p[0] ^= (mp_limb)z;
      }
      return bint_norm(mrb, bint_new(ctx, &c));
    }
  }
  y = mrb_as_bint(mrb, y);
  bint_as_mpz(RBIGINT(y), &b);
  if (zero_p(&a)) return y;
  if (zero_p(&b)) return x;
  mpz_init(ctx, &c);
  mpz_xor(ctx, &c, &a, &b);
  return bint_norm(mrb, bint_new(ctx, &c));
}

mrb_value
mrb_bint_neg(mrb_state *mrb, mrb_value x)
{
  mpz_t a, b;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  mpz_init(ctx, &b);
  mpz_neg(ctx, &b, &a);
  struct RBigint *b2 = bint_new(ctx, &b);
  /* no normalization */
  return mrb_obj_value(b2);
}

mrb_value
mrb_bint_rev(mrb_state *mrb, mrb_value x)
{
  mpz_t a, b;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  mpz_init(ctx, &b);
  mpz_neg(ctx, &b, &a);
  mpz_sub_int(ctx, &b, 1);
  return bint_norm(mrb, bint_new(ctx, &b));
}

mrb_value
mrb_bint_lshift(mrb_state *mrb, mrb_value x, mrb_int width)
{
  mpz_t a, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  mpz_init(ctx, &z);
  if (width < 0) {
    mpz_div_2exp(ctx, &z, &a, -width);
  }
  else {
    mpz_mul_2exp(ctx, &z, &a, width);
  }
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_rshift(mrb_state *mrb, mrb_value x, mrb_int width)
{
  mpz_t a, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  mpz_init(ctx, &z);
  if (width < 0) {
    mpz_mul_2exp(ctx, &z, &a, -width);
  }
  else {
    mpz_div_2exp(ctx, &z, &a, width);
  }
  return bint_norm(mrb, bint_new(ctx, &z));
}

void
mrb_bint_copy(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t b, temp;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(y), &b);
  mpz_init_set(ctx, &temp, &b);
  bint_set(ctx, RBIGINT(x), &temp);
}

size_t
mrb_bint_memsize(mrb_value x)
{
  mpz_t z;

  bint_as_mpz(RBIGINT(x), &z);
  return z.sz * sizeof(mp_limb);
}

mrb_value
mrb_bint_sqrt(mrb_state *mrb, mrb_value x)
{
  mpz_t a;

  bint_as_mpz(RBIGINT(x), &a);
  if (a.sn < 0) {
    mrb_raise(mrb, E_ARGUMENT_ERROR, "square root of negative number");
  }
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_t z;
  mpz_init(ctx, &z);
  mpz_sqrt(ctx, &z, &a);

  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_int
mrb_bint_sign(mrb_state *mrb, mrb_value bint)
{
  return RBIGINT_SIGN(RBIGINT(bint));
}

mrb_int
mrb_bint_size(mrb_state *mrb, mrb_value bint)
{
  mpz_t z;
  bint_as_mpz(RBIGINT(bint), &z);
  return z.sz * sizeof(mp_limb);
}

mrb_value
mrb_bint_from_bytes(mrb_state *mrb, const uint8_t *bytes, mrb_int len)
{
  mpz_t z;
  size_t limb_len = (len + sizeof(mp_limb) - 1) / sizeof(mp_limb);
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init_heap(ctx, &z, limb_len);
  memcpy(z.p, bytes, len);
  z.sn = (len > 0) ? 1 : 0;
  z.sz = limb_len;
  trim(&z);
  return bint_norm(mrb, bint_new(ctx, &z));
}

mrb_value
mrb_bint_hash(mrb_state *mrb, mrb_value x)
{
  mpz_t z;

  bint_as_mpz(RBIGINT(x), &z);
  uint32_t hash = mrb_byte_hash((uint8_t*)z.p, z.sz*sizeof(mp_limb));
  hash = mrb_byte_hash_step((uint8_t*)&z.sn, sizeof(z.sn), hash);
  return mrb_int_value(mrb, hash);
}

/* to be used only from mruby-sprintf */
mrb_value
mrb_bint_2comp(mrb_state *mrb, mrb_value x)
{
  mpz_t a, z;
  MPZ_CTX_INIT(mrb, ctx, pool);

  bint_as_mpz(RBIGINT(x), &a);
  mpz_init(ctx, &z);
  mrb_assert(a.sn < 0);
  size_t size = a.sz;
  mpz_realloc(ctx, &z, size);
  mp_limb *ds = a.p;
  mp_limb *dd = z.p;
  char carry = 1;
  for (size_t i=0; i<size; i++) {
    mp_limb xv = ds[i];
    make_2comp(xv, carry);
    dd[i] = xv;
  }
  z.sn = 1;

  struct RBigint *b2 = bint_new(ctx, &z);
  return mrb_obj_value(b2);
}

#ifdef MRB_USE_RATIONAL
void
mrb_bint_reduce(mrb_state *mrb, mrb_value *xp, mrb_value *yp)
{
  mpz_t r, x, y, a, b;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &r);
  mpz_init(ctx, &a); mpz_init(ctx, &b);

  bint_as_mpz(RBIGINT(*xp), &x);
  bint_as_mpz(RBIGINT(*yp), &y);

  mpz_gcd(ctx, &r, &x, &y);

  mpz_mdiv(ctx, &a, &x, &r);
  mpz_mdiv(ctx, &b, &y, &r);

  mpz_clear(ctx, &r);

  struct RBigint *b1 = bint_new(ctx, &a);
  struct RBigint *b2 = bint_new(ctx, &b);
  *xp = mrb_obj_value(b1);
  *yp = mrb_obj_value(b2);
}
#endif

mrb_value
mrb_bint_gcd(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t r, a, b;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &r);
  bint_as_mpz(RBIGINT(x), &a);
  bint_as_mpz(RBIGINT(y), &b);

  mpz_gcd(ctx, &r, &a, &b);

  struct RBigint *result = bint_new(ctx, &r);
  return bint_norm(mrb, result);
}

mrb_value
mrb_bint_lcm(mrb_state *mrb, mrb_value x, mrb_value y)
{
  mpz_t gcd_val, x_mpz, y_mpz, abs_x, abs_y, quotient, result_mpz;

  bint_as_mpz(RBIGINT(x), &x_mpz);
  bint_as_mpz(RBIGINT(y), &y_mpz);

  if (zero_p(&x_mpz) || zero_p(&y_mpz)) {
    return mrb_fixnum_value(0);
  }

  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_abs(ctx, &abs_x, &x_mpz);
  mpz_abs(ctx, &abs_y, &y_mpz);

  mpz_init(ctx, &gcd_val);
  mpz_gcd(ctx, &gcd_val, &abs_x, &abs_y);

  /* LCM = (abs_x / gcd) * abs_y — divide first to reduce intermediate size */
  mpz_init(ctx, &quotient);
  mpz_mdiv(ctx, &quotient, &abs_x, &gcd_val);
  mpz_init(ctx, &result_mpz);
  mpz_mul(ctx, &result_mpz, &quotient, &abs_y);

  mpz_clear(ctx, &gcd_val);
  mpz_clear(ctx, &abs_x);
  mpz_clear(ctx, &abs_y);
  mpz_clear(ctx, &quotient);

  struct RBigint *result = bint_new(ctx, &result_mpz);
  return bint_norm(mrb, result);
}

mrb_value
mrb_bint_abs(mrb_state *mrb, mrb_value x)
{
  mpz_t a, result_mpz;
  MPZ_CTX_INIT(mrb, ctx, pool);

  mpz_init(ctx, &result_mpz);
  bint_as_mpz(RBIGINT(x), &a);
  mpz_abs(ctx, &result_mpz, &a);

  struct RBigint *result = bint_new(ctx, &result_mpz);
  return mrb_obj_value(result);
}
