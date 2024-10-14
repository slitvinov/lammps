#include <cmath>
#include <direct.h>
#undef ATOBIGINT
#define ATOBIGINT _atoi64
#define pclose _pclose
#define strdup _strdup
#if !defined(__MINGW32__)
inline double pow(int i, int j) { return pow((double)i, j); }
inline double fabs(int i) { return fabs((double)i); }
inline double sqrt(int i) { return sqrt((double)i); }
#endif
inline double trunc(double x) { return x > 0 ? floor(x) : ceil(x); }
#ifndef S_IRWXU
#define S_IRWXU 0
#endif
#ifndef S_IRGRP
#define S_IRGRP 0
#endif
#ifndef S_IXGRP
#define S_IXGRP 0
#endif
inline int mkdir(const char *path, int) { return _mkdir(path); }
