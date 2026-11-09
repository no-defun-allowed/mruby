/*
 ** pack.c - Array#pack, String#unpack
 */

#include <mruby.h>
#include <mruby/error.h>
#include <mruby/array.h>
#include <mruby/class.h>
#include <mruby/numeric.h>
#include <mruby/string.h>
#include <mruby/variable.h>
#include <mruby/endian.h>
#include <mruby/internal.h>

#include <ctype.h>
#include <string.h>

#define INT_OVERFLOW_P(n)  ((n) < MRB_INT_MIN || (n) > MRB_INT_MAX)
#define UINT_OVERFLOW_P(n) ((n) > MRB_INT_MAX)

struct tmpl {
  mrb_value str;
  int idx;
};

enum pack_dir {
  PACK_DIR_NONE,
  PACK_DIR_CHAR,      /* C */
  PACK_DIR_SHORT,     /* S */
  PACK_DIR_LONG,      /* L */
  PACK_DIR_STR,       /* A */
};

enum pack_type {
  PACK_TYPE_INTEGER,
  PACK_TYPE_STRING,
  PACK_TYPE_NONE
};

#define PACK_FLAG_a             0x01      /* null padding ("a") */
#define PACK_FLAG_SIGNED        0x02      /* signed ("!") */
#define PACK_FLAG_WIDTH         0x04      /* "count" is "width" */
#define PACK_FLAG_LSB           0x08      /* LSB / low nibble first */
#define PACK_FLAG_COUNT2        0x10      /* "count" is special... */
#define PACK_FLAG_LITTLEENDIAN  0x20      /* little endian actually */

/* template parsing optimization structures */
typedef struct {
  unsigned char dir;
  unsigned char type;
  unsigned char size;
  unsigned char base_flags;
} format_info_t;

static char format_names[] = {'A','a','C','c','N','n','V','v'};
static format_info_t format_infos[] = {
  {PACK_DIR_STR, PACK_TYPE_STRING, 0, PACK_FLAG_WIDTH | PACK_FLAG_COUNT2},
  {PACK_DIR_STR, PACK_TYPE_STRING, 0, PACK_FLAG_WIDTH | PACK_FLAG_COUNT2 | PACK_FLAG_a},
  {PACK_DIR_CHAR, PACK_TYPE_INTEGER, 1, 0},
  {PACK_DIR_CHAR, PACK_TYPE_INTEGER, 1, PACK_FLAG_SIGNED},
  {PACK_DIR_LONG, PACK_TYPE_INTEGER, 4, 0},
  {PACK_DIR_SHORT, PACK_TYPE_INTEGER, 2, 0},
  {PACK_DIR_LONG, PACK_TYPE_INTEGER, 4, PACK_FLAG_LITTLEENDIAN},
  {PACK_DIR_SHORT, PACK_TYPE_INTEGER, 2, PACK_FLAG_LITTLEENDIAN},
};

/* Format character lookup function */
static inline format_info_t pack_format_info(unsigned char c) {
  char *p = memchr(format_names, c, sizeof(format_names));
  if (p) return format_infos[p - format_names];
  else return (format_info_t){PACK_DIR_NONE, PACK_TYPE_NONE, 0, 0};
}

static mrb_value
str_len_ensure(mrb_state *mrb, mrb_value str, mrb_int len)
{
  mrb_int n = RSTRING_LEN(str);
  if (len < 0) {
    mrb_raise(mrb, E_RANGE_ERROR, "negative (or overflowed) integer");
  }
  if (len > n) {
    do {
      n *= 2;
    } while (len > n);
    str = mrb_str_resize(mrb, str, n);
  }
  return str;
}

static void
pack_int(mrb_state *mrb, mrb_value o, mrb_value str, mrb_int sidx, unsigned int bytes, unsigned int flags) {
  str = str_len_ensure(mrb, str, sidx + bytes);
  uint32_t n = (uint32_t)mrb_integer(o);
  char *dptr = RSTRING_PTR(str) + sidx;
  if (flags & PACK_FLAG_LITTLEENDIAN)
    for (int i = 0; i < bytes; i++)
      dptr[i] = (char)(n >> (8 * i));
  else
    for (int i = 0; i < bytes; i++)
      dptr[bytes - i - 1] = (char)(n >> (8 * i));
}

#ifndef MRB_INT64
static void
u32tostr(char *buf, size_t len, uint32_t n)
{
#ifdef MRB_NO_STDIO
  char *bufend = buf + len;
  char *p = bufend - 1;

  if (len < 1) {
    return;
  }

  *p-- = '\0';
  len--;

  if (n > 0) {
    for (; len > 0 && n > 0; len --, n /= 10) {
      *p -- = '0' + (n % 10);
    }
    p++;
  }
  else if (len > 0) {
    *p = '0';
    len--;
  }

  memmove(buf, p, bufend - p);
#else
  snprintf(buf, len, "%" PRIu32, n);
#endif /* MRB_NO_STDIO */
}
#endif /* MRB_INT64 */

static void
unpack_int(mrb_state *mrb, const void *src, mrb_value ary, unsigned int bytes, unsigned int flags)
{
  uint32_t ul = 0;
  const char *p = src;
  /* read the value */
  if (flags & PACK_FLAG_LITTLEENDIAN)
    for (int i = 0; i < bytes; i++)
      ul |= p[i] << (8 * i);
  else
    for (int i = 0; i < bytes; i++)
      ul |= p[bytes - i - 1] << (8 * i);
  if (flags & PACK_FLAG_SIGNED) {
    mrb_int mask = 1 << (8 * bytes - 1);
    /* If the value is positive, this sets the high bit
     * (relative to the data size) and then removes it.
     * If the value is negative, this clears the high bit
     * and then sets all bits from that bit to the MSB. */
    mrb_int n = (ul ^ mask) - mask;
    /* It always fits in a 32-bit (signed) integer. */
    mrb_ary_push(mrb, ary, mrb_int_value(mrb, n));
  } else {
#ifndef MRB_INT64
    if (UINT_OVERFLOW_P(ul)) {
#ifdef MRB_USE_BIGINT
      /* Overflow the MSB into a bignum. */
      int ai = mrb_gc_arena_save(mrb);
      mrb_ary_push(mrb, ary, mrb_bint_new_uint64(mrb, ul));
      mrb_gc_arena_restore(mrb, ai);
      return;
#else
      u32tostr(msg, sizeof(msg), ul);
      mrb_raisef(mrb, E_RANGE_ERROR, "cannot unpack to Integer: %s", msg);
#endif
    }
#endif
    mrb_ary_push(mrb, ary, mrb_int_value(mrb, ul));
  }
}

static int
pack_str(mrb_state *mrb, mrb_value src, mrb_value dst, mrb_int didx, int count, unsigned int flags)
{
  const char *sptr = RSTRING_PTR(src);
  mrb_int slen = RSTRING_LEN(src);
  mrb_int copylen, padlen;
  char pad;

  /* determine padding character based on format */
  if (flags & PACK_FLAG_a)
    pad = '\0';
  else
    pad = ' ';

  if (count == 0) {
    return 0;
  }
  else if (count == -1) {
    copylen = slen;
    padlen = 0;
  }
  else if (count < slen) {
    copylen = count;
    padlen = 0;
  }
  else {
    copylen = slen;
    padlen = count - slen;
  }

  /* pre-allocate exact buffer size */
  dst = str_len_ensure(mrb, dst, didx + copylen + padlen);
  char *dptr = RSTRING_PTR(dst) + didx;
  char *dptr0 = dptr;

  /* copy string data */
  if (copylen > 0) {
    memcpy(dptr, sptr, copylen);
    dptr += copylen;
  }

  /* bulk padding using memset instead of loop */
  if (padlen > 0) {
    memset(dptr, pad, padlen);
    dptr += padlen;
  }

  return (int)(dptr - dptr0);
}

#define CHECK_UNPACK_LEN(mrb, slen, ary) do {\
  if ((slen) <= 0) {\
    mrb_ary_push(mrb, ary, mrb_str_new(mrb, 0, 0));\
    return 0;\
  }\
} while (0)

static inline mrb_bool is_padding_char_a(unsigned char c) {
  switch (c) {
  case '\0':
  case ' ':
  case '\t':
  case '\n':
  case '\v':
  case '\f':
  case '\r': return 1;
  default: return 0;
  }
}

static int
unpack_str(mrb_state *mrb, const void *src, mrb_int slen, mrb_value ary, int count, unsigned int flags)
{
  CHECK_UNPACK_LEN(mrb, slen, ary);

  const char *sptr = (const char*)src;
  int copylen;

  if (count != -1 && count < slen) {
    slen = count;
  }
  copylen = (int)slen;

  if (!(flags & PACK_FLAG_a)) {  /* "A" format - trim spaces and nulls */
    while (copylen > 0 && is_padding_char_a(sptr[copylen - 1])) {
      copylen--;
    }
  }
  /* "a" format does no trimming */

  if (copylen < 0) copylen = 0;

  mrb_value dst = mrb_str_new(mrb, sptr, (mrb_int)copylen);
  mrb_ary_push(mrb, ary, dst);
  return (int)slen;
}

static void
check_x(mrb_state *mrb, mrb_int a, mrb_int count, char c)
{
  if (a < count) {
    mrb_raisef(mrb, E_ARGUMENT_ERROR, "%c outside of string", c);
  }
}

static void
prepare_tmpl(mrb_state *mrb, struct tmpl *tmpl)
{
  mrb_get_args(mrb, "S", &tmpl->str);
  tmpl->idx = 0;
}

static int
has_tmpl(const struct tmpl *tmpl)
{
  return (tmpl->idx < RSTRING_LEN(tmpl->str));
}

static enum pack_dir
read_tmpl(mrb_state *mrb, struct tmpl *tmpl, enum pack_type *typep, mrb_int *sizep, mrb_int *countp, unsigned int *flagsp)
{
  mrb_int t, tlen;
  int ch, size = 0;
  enum pack_dir dir;
  enum pack_type type;
  mrb_int count = 1;
  unsigned int flags = 0;
  const char *tptr;

  tptr = RSTRING_PTR(tmpl->str);
  tlen = RSTRING_LEN(tmpl->str);

 restart:
  if (tmpl->idx >= tlen) return PACK_DIR_NONE;
  t = tptr[tmpl->idx++];
  /* Handle whitespace - skip and restart */
  if (ISSPACE((char)t)) {
    goto restart;
  }

  if (t >= 0 && t < 256) {
    format_info_t info_val = pack_format_info((unsigned char)t);
    const format_info_t *info = &info_val;
    if (info->dir != PACK_DIR_NONE) {
      /* Valid format character found in lookup table */
      dir = info->dir;
      type = info->type;
      size = info->size;
      flags = info->base_flags;
    } else {
      char c = (char)t;
      mrb_value s = mrb_str_new(mrb, &c, 1);
      mrb_raisef(mrb, E_ARGUMENT_ERROR, "unknown unpack directive %!v", s);
    }
  }

  /* read suffix [0-9*!] */
  while (tmpl->idx < tlen) {
    ch = tptr[tmpl->idx];
    if (ISDIGIT(ch)) {
      char *e;
      mrb_int n;
      if (!mrb_read_int(tptr+tmpl->idx, tptr+tlen, &e, &n) || INT_MAX < n) {
        mrb_raise(mrb, E_RUNTIME_ERROR, "too big template length");
      }
      count = n;
      tmpl->idx = (int)(e - tptr);
      continue;
    }
    else if (ch == '*') {
      count = -1;
    }
    else if (ch == '!') {
      flags |= PACK_FLAG_SIGNED;
    }
    else {
      break;
    }
    tmpl->idx++;
  }

  *typep = type;
  *sizep = size;
  *countp = count;
  *flagsp = flags;
  return dir;
}

/*
 * call-seq:
 *   array.pack(template) -> string
 *
 * Packs the contents of array into a binary string according to the
 * directives in template. Directives are single characters from the
 * table below. Each directive may be followed by a number indicating
 * the number of times to repeat the directive.
 *
 * Template string directives:
 *   C - 8-bit unsigned integer (unsigned char)
 *   c - 8-bit signed integer (signed char)
 *   S - 16-bit unsigned integer, native endian (uint16_t)
 *   s - 16-bit signed integer, native endian (int16_t)
 *   L - 32-bit unsigned integer, native endian (uint32_t)
 *   l - 32-bit signed integer, native endian (int32_t)
 *   Q - 64-bit unsigned integer, native endian (uint64_t)
 *   q - 64-bit signed integer, native endian (int64_t)
 *   n - 16-bit unsigned integer, network byte order
 *   N - 32-bit unsigned integer, network byte order
 *   v - 16-bit unsigned integer, little endian
 *   V - 32-bit unsigned integer, little endian
 *   f - single precision float, native format
 *   d - double precision float, native format
 *   A - ASCII string, space padded
 *   a - ASCII string, null padded
 *   Z - null-terminated string
 *   H - hex string, high nibble first
 *   h - hex string, low nibble first
 *   x - null byte
 *   X - back up one byte
 *   @ - null fill to absolute position
 *
 *   [1, 2, 3].pack("CCC")     #=> "\x01\x02\x03"
 *   [1, 2].pack("S*")         #=> "\x01\x00\x02\x00" (little endian)
 *   ["hello"].pack("A10")     #=> "hello     "
 */
static mrb_value
mrb_pack_pack(mrb_state *mrb, mrb_value ary)
{
  mrb_value o;
  struct tmpl tmpl;
  enum pack_type type;
  mrb_int count;
  mrb_int size;
  unsigned int flags;
  enum pack_dir dir;

  prepare_tmpl(mrb, &tmpl);

  mrb_value result = mrb_str_new(mrb, NULL, 128);  /* allocate initial buffer */
  mrb_int aidx = 0;
  mrb_int ridx = 0;
  while (has_tmpl(&tmpl)) {
    dir = read_tmpl(mrb, &tmpl, &type, &size, &count, &flags);
    if (dir == PACK_DIR_NONE) break;

    if ((flags & PACK_FLAG_WIDTH) && aidx >= RARRAY_LEN(ary)) {
      mrb_raise(mrb, E_ARGUMENT_ERROR, "too few arguments");
    }
    for (; aidx < RARRAY_LEN(ary); aidx++) {
      if (count == 0 && !(flags & PACK_FLAG_WIDTH))
        break;
      
      o = RARRAY_PTR(ary)[aidx];
      if (type == PACK_TYPE_INTEGER) {
        o = mrb_ensure_int_type(mrb, o);
        pack_int(mrb, o, result, ridx, size, flags);
        ridx += size;
      }
      else if (type == PACK_TYPE_STRING) {
        if (!mrb_string_p(o)) {
          mrb_raisef(mrb, E_TYPE_ERROR, "can't convert %T into String", o);
        }
        ridx += pack_str(mrb, o, result, ridx, (int)count, flags);
      }
      if (flags & PACK_FLAG_COUNT2) {
        /* always consumes 1 entry */
        aidx++;
        break;
      }
      if (count > 0) {
        count--;
      }
    }
    if (ridx < 0) {
    overflow:
      mrb_raise(mrb, E_RANGE_ERROR, "negative (or overflowed) template size");
    }
  }

  mrb_str_resize(mrb, result, ridx);
  return result;
}

static mrb_value
pack_unpack(mrb_state *mrb, mrb_value str, mrb_bool single)
{
  struct tmpl tmpl;
  mrb_int count;
  unsigned int flags;
  enum pack_type type;
  mrb_int size;
  const unsigned char *sptr;

  prepare_tmpl(mrb, &tmpl);

  mrb_int srcidx = 0;
  mrb_int srclen = RSTRING_LEN(str);

  mrb_value result = mrb_ary_new(mrb);
  while (has_tmpl(&tmpl)) {
    enum pack_dir dir = read_tmpl(mrb, &tmpl, &type, &size, &count, &flags);

    if (dir == PACK_DIR_NONE) break;
    sptr = (const unsigned char*)RSTRING_PTR(str) + srcidx;
    if (dir == PACK_DIR_STR) {
      srcidx += unpack_str(mrb, sptr, srclen - srcidx, result, (int)count, flags);
      if (single) goto single_return;
      continue;
    }

    while (count != 0 && srcidx < srclen) {
      if (srclen - srcidx < size) {
        while (count-- > 0) {
          mrb_ary_push(mrb, result, mrb_nil_value());
        }
        break;
      }

      sptr = (const unsigned char*)RSTRING_PTR(str) + srcidx;
      unpack_int(mrb, sptr, result, size, flags);
      srcidx += size;
      if (count > 0) {
        count--;
      }
    }
  }
  if (single) {
  single_return:
    if (RARRAY_LEN(result) > 0) {
      return RARRAY_PTR(result)[0];
    }
    return mrb_nil_value();
  }
  return result;
}

/*
 * call-seq:
 *   string.unpack(template) -> array
 *
 * Unpacks the contents of string according to the template string,
 * returning an array of values. Template uses the same format as
 * Array#pack. See Array#pack for template string format.
 *
 *   "\x01\x02\x03".unpack("CCC")     #=> [1, 2, 3]
 *   "\x01\x00\x02\x00".unpack("S*")  #=> [1, 2] (little endian)
 *   "hello     ".unpack("A10")       #=> ["hello"]
 */
static mrb_value
mrb_pack_unpack(mrb_state *mrb, mrb_value str)
{
  return pack_unpack(mrb, str, FALSE);
}

/*
 * call-seq:
 *   string.unpack1(template) -> object
 *
 * Unpacks the first value from string according to the template string.
 * This is equivalent to string.unpack(template)[0] but more efficient
 * when only the first value is needed.
 *
 *   "\x01\x02\x03".unpack1("C")      #=> 1
 *   "\x01\x00\x02\x00".unpack1("S")  #=> 1 (little endian)
 *   "hello     ".unpack1("A10")      #=> "hello"
 */
static mrb_value
mrb_pack_unpack1(mrb_state *mrb, mrb_value str)
{
  return pack_unpack(mrb, str, TRUE);
}

void
mrb_mruby_mini_pack_gem_init(mrb_state *mrb)
{
  mrb_define_method_id(mrb, mrb->array_class, MRB_SYM(pack), mrb_pack_pack, MRB_ARGS_REQ(1));
  mrb_define_method_id(mrb, mrb->string_class, MRB_SYM(unpack), mrb_pack_unpack, MRB_ARGS_REQ(1));
  mrb_define_method_id(mrb, mrb->string_class, MRB_SYM(unpack1), mrb_pack_unpack1, MRB_ARGS_REQ(1));
}

void
mrb_mruby_mini_pack_gem_final(mrb_state *mrb)
{
}
