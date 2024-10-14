#include "hashlittle.h"
#if defined(__BYTE_ORDER__)
#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define HASH_LITTLE_ENDIAN 1
#else
#define HASH_LITTLE_ENDIAN 0
#endif
#else
#if defined(__bg__)
#define HASH_LITTLE_ENDIAN 0
#else
#define HASH_LITTLE_ENDIAN 1
#endif
#endif
#define rot(x, k) (((x) << (k)) | ((x) >> (32 - (k))))
#define mix(a, b, c)                                                           \
  {                                                                            \
    a -= c;                                                                    \
    a ^= rot(c, 4);                                                            \
    c += b;                                                                    \
    b -= a;                                                                    \
    b ^= rot(a, 6);                                                            \
    a += c;                                                                    \
    c -= b;                                                                    \
    c ^= rot(b, 8);                                                            \
    b += a;                                                                    \
    a -= c;                                                                    \
    a ^= rot(c, 16);                                                           \
    c += b;                                                                    \
    b -= a;                                                                    \
    b ^= rot(a, 19);                                                           \
    a += c;                                                                    \
    c -= b;                                                                    \
    c ^= rot(b, 4);                                                            \
    b += a;                                                                    \
  }
#define final(a, b, c)                                                         \
  {                                                                            \
    c ^= b;                                                                    \
    c -= rot(b, 14);                                                           \
    a ^= c;                                                                    \
    a -= rot(c, 11);                                                           \
    b ^= a;                                                                    \
    b -= rot(a, 25);                                                           \
    c ^= b;                                                                    \
    c -= rot(b, 16);                                                           \
    a ^= c;                                                                    \
    a -= rot(c, 4);                                                            \
    b ^= a;                                                                    \
    b -= rot(a, 14);                                                           \
    c ^= b;                                                                    \
    c -= rot(b, 24);                                                           \
  }
uint32_t LAMMPS_NS::hashlittle(const void *key, size_t length,
                               uint32_t initval) {
#ifndef PURIFY_HATES_HASHLITTLE
  uint32_t a, b, c;
  union {
    const void *ptr;
    size_t i;
  } u;
  a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;
  u.ptr = key;
  if (HASH_LITTLE_ENDIAN && ((u.i & 0x3) == 0)) {
    const uint32_t *k = (const uint32_t *)key;
    while (length > 12) {
      a += k[0];
      b += k[1];
      c += k[2];
      mix(a, b, c);
      length -= 12;
      k += 3;
    }
#ifndef VALGRIND
    switch (length) {
    case 12:
      c += k[2];
      b += k[1];
      a += k[0];
      break;
    case 11:
      c += k[2] & 0xffffff;
      b += k[1];
      a += k[0];
      break;
    case 10:
      c += k[2] & 0xffff;
      b += k[1];
      a += k[0];
      break;
    case 9:
      c += k[2] & 0xff;
      b += k[1];
      a += k[0];
      break;
    case 8:
      b += k[1];
      a += k[0];
      break;
    case 7:
      b += k[1] & 0xffffff;
      a += k[0];
      break;
    case 6:
      b += k[1] & 0xffff;
      a += k[0];
      break;
    case 5:
      b += k[1] & 0xff;
      a += k[0];
      break;
    case 4:
      a += k[0];
      break;
    case 3:
      a += k[0] & 0xffffff;
      break;
    case 2:
      a += k[0] & 0xffff;
      break;
    case 1:
      a += k[0] & 0xff;
      break;
    case 0:
      return c;
    }
#else
    const uint8_t *k8 = (const uint8_t *)k;
    switch (length) {
    case 12:
      c += k[2];
      b += k[1];
      a += k[0];
      break;
    case 11:
      c += ((uint32_t)k8[10]) << 16;
    case 10:
      c += ((uint32_t)k8[9]) << 8;
    case 9:
      c += k8[8];
    case 8:
      b += k[1];
      a += k[0];
      break;
    case 7:
      b += ((uint32_t)k8[6]) << 16;
    case 6:
      b += ((uint32_t)k8[5]) << 8;
    case 5:
      b += k8[4];
    case 4:
      a += k[0];
      break;
    case 3:
      a += ((uint32_t)k8[2]) << 16;
    case 2:
      a += ((uint32_t)k8[1]) << 8;
    case 1:
      a += k8[0];
      break;
    case 0:
      return c;
    }
#endif
  } else if (HASH_LITTLE_ENDIAN && ((u.i & 0x1) == 0)) {
    const uint16_t *k = (const uint16_t *)key;
    const uint8_t *k8;
    while (length > 12) {
      a += k[0] + (((uint32_t)k[1]) << 16);
      b += k[2] + (((uint32_t)k[3]) << 16);
      c += k[4] + (((uint32_t)k[5]) << 16);
      mix(a, b, c);
      length -= 12;
      k += 6;
    }
    k8 = (const uint8_t *)k;
    switch (length) {
    case 12:
      c += k[4] + (((uint32_t)k[5]) << 16);
      b += k[2] + (((uint32_t)k[3]) << 16);
      a += k[0] + (((uint32_t)k[1]) << 16);
      break;
    case 11:
      c += ((uint32_t)k8[10]) << 16;
    case 10:
      c += k[4];
      b += k[2] + (((uint32_t)k[3]) << 16);
      a += k[0] + (((uint32_t)k[1]) << 16);
      break;
    case 9:
      c += k8[8];
    case 8:
      b += k[2] + (((uint32_t)k[3]) << 16);
      a += k[0] + (((uint32_t)k[1]) << 16);
      break;
    case 7:
      b += ((uint32_t)k8[6]) << 16;
    case 6:
      b += k[2];
      a += k[0] + (((uint32_t)k[1]) << 16);
      break;
    case 5:
      b += k8[4];
    case 4:
      a += k[0] + (((uint32_t)k[1]) << 16);
      break;
    case 3:
      a += ((uint32_t)k8[2]) << 16;
    case 2:
      a += k[0];
      break;
    case 1:
      a += k8[0];
      break;
    case 0:
      return c;
    }
  } else {
    const uint8_t *k = (const uint8_t *)key;
    while (length > 12) {
      a += k[0];
      a += ((uint32_t)k[1]) << 8;
      a += ((uint32_t)k[2]) << 16;
      a += ((uint32_t)k[3]) << 24;
      b += k[4];
      b += ((uint32_t)k[5]) << 8;
      b += ((uint32_t)k[6]) << 16;
      b += ((uint32_t)k[7]) << 24;
      c += k[8];
      c += ((uint32_t)k[9]) << 8;
      c += ((uint32_t)k[10]) << 16;
      c += ((uint32_t)k[11]) << 24;
      mix(a, b, c);
      length -= 12;
      k += 12;
    }
    switch (length) {
    case 12:
      c += ((uint32_t)k[11]) << 24;
    case 11:
      c += ((uint32_t)k[10]) << 16;
    case 10:
      c += ((uint32_t)k[9]) << 8;
    case 9:
      c += k[8];
    case 8:
      b += ((uint32_t)k[7]) << 24;
    case 7:
      b += ((uint32_t)k[6]) << 16;
    case 6:
      b += ((uint32_t)k[5]) << 8;
    case 5:
      b += k[4];
    case 4:
      a += ((uint32_t)k[3]) << 24;
    case 3:
      a += ((uint32_t)k[2]) << 16;
    case 2:
      a += ((uint32_t)k[1]) << 8;
    case 1:
      a += k[0];
      break;
    case 0:
      return c;
    }
  }
  final(a, b, c);
  return c;
#else
#define MAXINT_DIV_PHI 11400714819323198485U
  uint32_t h, rest, *p, bytes, num_bytes;
  char *byteptr;
  num_bytes = length;
  h = 0;
  for (p = (uint32_t *)key, bytes = num_bytes;
       bytes >= (uint32_t)sizeof(uint32_t); bytes -= sizeof(uint32_t), p++) {
    h = (h ^ (*p)) * MAXINT_DIV_PHI;
  }
  rest = 0;
  for (byteptr = (char *)p; bytes > 0; bytes--, byteptr++) {
    rest = (rest << 8) | (*byteptr);
  }
  if (rest)
    h = (h ^ rest) * MAXINT_DIV_PHI;
  return h;
#endif
}
