#include <mruby.h>

#ifndef MRB_NO_FLOAT

#include <string.h>
#include <math.h>
#include <stdint.h>

static double
mrb_pow10(int exp)
{
  return pow(10, (double)exp);
}

/*
** Parses a string representation of a floating-point number.
**
** @param str The input string to parse.
** @param endp A pointer to a char* that will be updated to point to the
**             character in str after the last character used in the
**             conversion. Can be NULL.
** @param fp A pointer to a double that will be set to the parsed
**           floating-point number.
** @return TRUE if a float is successfully parsed, FALSE otherwise.
*/
MRB_API mrb_bool
mrb_read_float(const char *str, char **endp, double *fp)
{
  const char *p = str;
  const char *a = p;
  uint64_t int_part = 0;
  uint64_t frac_part = 0;
  int frac_digits = 0;
  int sign = 1;
  int digits = 0;
  int overflow = 0;

  // Skip whitespace
  while (ISSPACE((unsigned char)*p)) p++;

  // Handle sign
  if (*p == '-') { sign = -1; p++; }
  else if (*p == '+') p++;

  // Parse integer part using integer arithmetic for better accuracy
  while (ISDIGIT(*p)) {
    if (int_part > (UINT64_MAX - 9) / 10) {
      overflow = 1;
      // Continue parsing to find the end
      while (ISDIGIT(*p)) p++;
      break;
    }
    int_part = int_part * 10 + (*p - '0');
    digits++;
    a = ++p;
  }

  // Parse fractional part
  if (*p == '.') {
    p++;
    while (ISDIGIT(*p) && frac_digits < 17) { // Limit precision to avoid overflow
      frac_part = frac_part * 10 + (*p - '0');
      frac_digits++;
      digits++;
      p++;
    }
    // Skip remaining fractional digits if any
    while (ISDIGIT(*p)) p++;
    a = p;
  }

  // If no digits were found, return FALSE
  if (digits == 0) {
    if (endp) *endp = (char*)str;
    *fp = 0.0;
    return FALSE;
  }

  double res;
  if (overflow) {
    // For overflow case, fall back to floating point parsing
    res = (double)int_part;
  }
  else {
    // Fast path: combine integer and fractional parts
    res = (double)int_part;
    if (frac_digits > 0) {
      res += (double)frac_part * mrb_pow10(-frac_digits);
    }
  }

  // Handle exponent
  if ((*p | 32) == 'e') {
    int e = 0;
    int exp_sign = 1;
    p++;
    if (*p == '-') { exp_sign = -1; p++; }
    else if (*p == '+') p++;

    // If no digits follow 'e', ignore the exponent part
    if (!ISDIGIT(*p)) goto done;

    while (ISDIGIT(*p)) {
      if (e < 10000) { // Prevent integer overflow
        e = e * 10 + (*p - '0');
      }
      p++;
    }
    e *= exp_sign;

    // Apply exponent directly - let mrb_pow10 handle overflow
    // Large exponents will return HUGE_VAL (infinity) or 0.0 as appropriate

    res *= mrb_pow10(e);
    a = p;
  }

  // Apply sign
  res *= sign;

  // Set endp
 done:
  if (endp) *endp = (char*)a;
  *fp = res;

  // strtod(3) stores ERANGE to errno for overflow/underflow
  // mruby does not require those checks
#if 0
  // Check for underflow after applying the exponent
  if (res != 0.0 && fabs(res) < DBL_MIN) {
    return FALSE;
  }

  // Check if the result is infinity or NaN
  if (isinf(res) || isnan(res)) {
    return FALSE;
  }
#endif
  return TRUE;
}

#endif
