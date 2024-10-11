#ifndef LAMMPS_OMP_COMPAT
# if defined(__INTEL_LLVM_COMPILER)
#define LAMMPS_OMP_COMPAT 4
# elif defined(__INTEL_COMPILER)
# if __INTEL_COMPILER > 18
#define LAMMPS_OMP_COMPAT 4
# endif
# elif defined(__clang__)
# if __clang_major__ >= 10
#define LAMMPS_OMP_COMPAT 4
# endif
# elif defined(__GNUC__)
# if __GNUC__ >= 9
#define LAMMPS_OMP_COMPAT 4
# endif
# endif
#endif
#if LAMMPS_OMP_COMPAT == 4
#define LMP_SHARED(...) 
#define LMP_DEFAULT_NONE default(shared)
#else
#define LMP_SHARED(...) shared(__VA_ARGS__)
#define LMP_DEFAULT_NONE default(none)
#endif
