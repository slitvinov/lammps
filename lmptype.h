#ifndef LMP_LMPTYPE_H
#define LMP_LMPTYPE_H
#if __cplusplus < 201103L
#error LAMMPS requires a C++11 (or later) compliant compiler. Enable C++11 compatibility or upgrade the compiler.
#endif
#ifndef __STDC_LIMIT_MACROS
#define __STDC_LIMIT_MACROS
#endif
#ifndef __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif
#include <cinttypes>
#include <climits>
#include <cstdint>
#include <cstdlib>
#ifndef PRId64
#define PRId64 "ld"
#endif
namespace LAMMPS_NS {
#define SBBITS 30
#define HISTBITS 29
#define NEIGHMASK 0x1FFFFFFF
#define HISTMASK 0xDFFFFFFF
#define SPECIALMASK 0x3FFFFFFF
#if !defined(LAMMPS_SMALLSMALL) && !defined(LAMMPS_BIGBIG) &&                  \
    !defined(LAMMPS_SMALLBIG)
#define LAMMPS_SMALLBIG
#endif
#ifdef LAMMPS_LONGLONG_TO_LONG
#define MPI_LL MPI_LONG
#define ATOLL atoll
#else
#define MPI_LL MPI_LONG_LONG
#define ATOLL atol
#endif
#ifdef LAMMPS_SMALLBIG
typedef int smallint;
typedef int imageint;
typedef int tagint;
typedef int64_t bigint;
#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT64_MAX
#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_LL
#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%" PRId64
#define ATOTAGINT atoi
#define ATOBIGINT ATOLL
#define LAMMPS_TAGINT LAMMPS_INT
#define LAMMPS_TAGINT_2D LAMMPS_INT_2D
#define LAMMPS_BIGINT LAMMPS_INT64
#define LAMMPS_BIGINT_2D LAMMPS_INT64_2D
#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20
#endif
#ifdef LAMMPS_BIGBIG
typedef int smallint;
typedef int64_t imageint;
typedef int64_t tagint;
typedef int64_t bigint;
#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT64_MAX
#define MAXBIGINT INT64_MAX
#define MPI_LMP_TAGINT MPI_LL
#define MPI_LMP_IMAGEINT MPI_LL
#define MPI_LMP_BIGINT MPI_LL
#define TAGINT_FORMAT "%" PRId64
#define BIGINT_FORMAT "%" PRId64
#define ATOTAGINT ATOLL
#define ATOBIGINT ATOLL
#define LAMMPS_TAGINT LAMMPS_INT64
#define LAMMPS_TAGINT_2D LAMMPS_INT64_2D
#define LAMMPS_BIGINT LAMMPS_INT64
#define LAMMPS_BIGINT_2D LAMMPS_INT64_2D
#define IMGMASK 2097151
#define IMGMAX 1048576
#define IMGBITS 21
#define IMG2BITS 42
#endif
#ifdef LAMMPS_SMALLSMALL
typedef int smallint;
typedef int imageint;
typedef int tagint;
typedef int bigint;
#define MAXSMALLINT INT_MAX
#define MAXTAGINT INT_MAX
#define MAXBIGINT INT_MAX
#define MPI_LMP_TAGINT MPI_INT
#define MPI_LMP_IMAGEINT MPI_INT
#define MPI_LMP_BIGINT MPI_INT
#define TAGINT_FORMAT "%d"
#define BIGINT_FORMAT "%d"
#define ATOTAGINT atoi
#define ATOBIGINT atoi
#define LAMMPS_TAGINT LAMMPS_INT
#define LAMMPS_TAGINT_2D LAMMPS_INT_2D
#define LAMMPS_BIGINT LAMMPS_INT
#define LAMMPS_BIGINT_2D LAMMPS_INT_2D
#define IMGMASK 1023
#define IMGMAX 512
#define IMGBITS 10
#define IMG2BITS 20
#endif
union ubuf {
  double d;
  int64_t i;
  ubuf(const double &arg) : d(arg) {}
  ubuf(const int64_t &arg) : i(arg) {}
  ubuf(const int &arg) : i(arg) {}
};
} // namespace LAMMPS_NS
#ifdef _alignvar
#undef _alignvar
#endif
#ifdef _noalias
#undef _noalias
#endif
#ifdef _noopt
#undef _noopt
#endif
#if defined(__INTEL_COMPILER)
#define _alignvar(expr, val) __declspec(align(val)) expr
#elif defined(__GNUC__) || defined(__PGI) || defined(__INTEL_LLVM_COMPILER)
#define _alignvar(expr, val) expr __attribute((aligned(val)))
#else
#define _alignvar(expr, val) expr
#endif
#if defined(__INTEL_COMPILER) || (defined(__PGI) && !defined(__NVCOMPILER))
#define _noalias restrict
#elif defined(__GNUC__) || defined(__INTEL_LLVM_COMPILER) ||                   \
    defined(__NVCOMPILER)
#define _noalias __restrict
#else
#define _noalias
#endif
#if defined(__clang__)
#define _noopt __attribute__((optnone))
#elif defined(__INTEL_COMPILER) || defined(__INTEL_LLVM_COMPILER)
#define _noopt
#elif defined(__PGI)
#define _noopt
#elif defined(__GNUC__)
#if (__GNUC__ > 4) || ((__GNUC__ == 4) && (__GNUC_MINOR__ >= 9))
#if defined(_FORTIFY_SOURCE) && (_FORTIFY_SOURCE > 0)
#define _noopt __attribute__((optimize("no-var-tracking-assignments")))
#else
#define _noopt __attribute__((optimize("O0", "no-var-tracking-assignments")))
#endif
#else
#if defined(_FORTIFY_SOURCE) && (_FORTIFY_SOURCE > 0)
#define _noopt
#else
#define _noopt __attribute__((optimize("O0")))
#endif
#endif
#else
#define _noopt
#endif
#define LMP_UNUSED_PARAM(x) (void)(x)
#endif
