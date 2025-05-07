#ifndef LMP_MEMORY_H
#define LMP_MEMORY_H
namespace LAMMPS_NS {
class Memory : protected Pointers {
public:
  Memory(class LAMMPS *);
  void *smalloc(bigint n, const char *);
  void *srealloc(void *, bigint n, const char *);
  void sfree(void *);
  void fail(const char *);
  template <typename TYPE> TYPE *create(TYPE *&array, int n, const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * n;
    array = (TYPE *)smalloc(nbytes, name);
    return array;
  }
  template <typename TYPE> TYPE **create(TYPE **&, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> TYPE *grow(TYPE *&array, int n, const char *name) {
    if (array == nullptr)
      return create(array, n, name);
    bigint nbytes = ((bigint)sizeof(TYPE)) * n;
    array = (TYPE *)srealloc(array, nbytes, name);
    return array;
  }
  template <typename TYPE> TYPE **grow(TYPE **&, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy(TYPE *&array) {
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE *create1d_offset(TYPE *&array, int nlo, int nhi, const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * (nhi - nlo + 1);
    array = (TYPE *)smalloc(nbytes, name);
    array -= nlo;
    return array;
  }
  template <typename TYPE>
  TYPE **create1d_offset(TYPE **&, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy1d_offset(TYPE *&array, int offset) {
    if (array)
      sfree(&array[offset]);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE **create(TYPE **&array, int n1, int n2, const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2;
    TYPE *data = (TYPE *)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1;
    array = (TYPE **)smalloc(nbytes, name);
    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }
  template <typename TYPE>
  TYPE ***create(TYPE ***&, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  TYPE **grow(TYPE **&array, int n1, int n2, const char *name) {
    if (array == nullptr)
      return create(array, n1, n2, name);
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2;
    TYPE *data = (TYPE *)srealloc(array[0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1;
    array = (TYPE **)srealloc(array, nbytes, name);
    bigint n = 0;
    for (int i = 0; i < n1; i++) {
      array[i] = &data[n];
      n += n2;
    }
    return array;
  }
  template <typename TYPE> TYPE ***grow(TYPE ***&, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy(TYPE **&array) {
    if (array == nullptr)
      return;
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ***create(TYPE ***&array, int n1, int n2, int n3, const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2 * n3;
    TYPE *data = (TYPE *)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1 * n2;
    TYPE **plane = (TYPE **)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE **)) * n1;
    array = (TYPE ***)smalloc(nbytes, name);
    int i, j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint)i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m + j] = &data[n];
        n += n3;
      }
    }
    return array;
  }
  template <typename TYPE>
  TYPE ****create(TYPE ****&, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  TYPE ***grow(TYPE ***&array, int n1, int n2, int n3, const char *name) {
    if (array == nullptr)
      return create(array, n1, n2, n3, name);
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2 * n3;
    TYPE *data = (TYPE *)srealloc(array[0][0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1 * n2;
    TYPE **plane = (TYPE **)srealloc(array[0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE **)) * n1;
    array = (TYPE ***)srealloc(array, nbytes, name);
    int i, j;
    bigint m;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m = ((bigint)i) * n2;
      array[i] = &plane[m];
      for (j = 0; j < n2; j++) {
        plane[m + j] = &data[n];
        n += n3;
      }
    }
    return array;
  }
  template <typename TYPE>
  TYPE ****grow(TYPE ****&, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy(TYPE ***&array) {
    if (array == nullptr)
      return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, int n2, int n3,
                          const char *name) {
    if (n1lo > n1hi || n2 <= 0 || n3 <= 0)
      return nullptr;
    int n1 = n1hi - n1lo + 1;
    create(array, n1, n2, n3, name);
    array -= n1lo;
    return array;
  }
  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****&, int, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy3d_offset(TYPE ***&array, int offset) {
    if (array == nullptr)
      return;
    sfree(&array[offset][0][0]);
    sfree(&array[offset][0]);
    sfree(&array[offset]);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ***create3d_offset(TYPE ***&array, int n1lo, int n1hi, int n2lo,
                          int n2hi, int n3lo, int n3hi, const char *name) {
    if (n1lo > n1hi || n2lo > n2hi || n3lo > n3hi)
      return nullptr;
    int n1 = n1hi - n1lo + 1;
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    create(array, n1, n2, n3, name);
    bigint m = ((bigint)n1) * n2;
    for (bigint i = 0; i < m; i++)
      array[0][i] -= n3lo;
    for (int i = 0; i < n1; i++)
      array[i] -= n2lo;
    array -= n1lo;
    return array;
  }
  template <typename TYPE>
  TYPE ****create3d_offset(TYPE ****&, int, int, int, int, int, int,
                           const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  void destroy3d_offset(TYPE ***&array, int n1_offset, int n2_offset,
                        int n3_offset) {
    if (array == nullptr)
      return;
    sfree(&array[n1_offset][n2_offset][n3_offset]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(&array[n1_offset]);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ***create3d_offset_last(TYPE ***&array, int n1lo, int n1hi, int n2lo,
                               int n2hi, int n3, const char *name) {
    if (n1lo > n1hi || n2lo > n2hi || n3 < 0)
      return nullptr;
    int n1 = n1hi - n1lo + 1;
    int n2 = n2hi - n2lo + 1;
    create(array, n1, n2, n3, name);
    for (int i = 0; i < n1; i++)
      array[i] -= n2lo;
    array -= n1lo;
    return array;
  }
  template <typename TYPE>
  TYPE ***create3d_offset_last(TYPE ****&, int, int, int, int, int,
                               const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  void destroy3d_offset_last(TYPE ***&array, int n1_offset, int n2_offset) {
    if (array == nullptr)
      return;
    sfree(&array[n1_offset][n2_offset][0]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(&array[n1_offset]);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ****create(TYPE ****&array, int n1, int n2, int n3, int n4,
                  const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2 * n3 * n4;
    TYPE *data = (TYPE *)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1 * n2 * n3;
    TYPE **cube = (TYPE **)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE **)) * n1 * n2;
    TYPE ***plane = (TYPE ***)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE ***)) * n1;
    array = (TYPE ****)smalloc(nbytes, name);
    bigint i, j, k;
    bigint m1, m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = i * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = i * n2 + j;
        m2 = i * n2 * n3 + j * n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = i * n2 * n3 + j * n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }
  template <typename TYPE>
  TYPE *****create(TYPE *****&, int, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  TYPE ****grow(TYPE ****&array, int n1, int n2, int n3, int n4,
                const char *name) {
    if (array == nullptr)
      return create(array, n1, n2, n3, n4, name);
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2 * n3 * n4;
    TYPE *data = (TYPE *)srealloc(array[0][0][0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1 * n2 * n3;
    TYPE **cube = (TYPE **)srealloc(array[0][0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE **)) * n1 * n2;
    TYPE ***plane = (TYPE ***)srealloc(array[0], nbytes, name);
    nbytes = ((bigint)sizeof(TYPE ***)) * n1;
    array = (TYPE ****)srealloc(array, nbytes, name);
    int i, j, k;
    bigint m1, m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint)i) * n2;
      array[i] = &plane[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint)i) * n2 + j;
        m2 = ((bigint)i) * n2 * n3 + j * n3;
        plane[m1] = &cube[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint)i) * n2 * n3 + j * n3 + k;
          cube[m1] = &data[n];
          n += n4;
        }
      }
    }
    return array;
  }
  template <typename TYPE>
  TYPE *****grow(TYPE *****&, int, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy(TYPE ****&array) {
    if (array == nullptr)
      return;
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ****create4d_offset(TYPE ****&array, int n1, int n2lo, int n2hi,
                           int n3lo, int n3hi, int n4lo, int n4hi,
                           const char *name) {
    if (n1 <= 0 || n2lo > n2hi || n3lo > n3hi || n4lo > n4hi)
      return nullptr;
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    int n4 = n4hi - n4lo + 1;
    create(array, n1, n2, n3, n4, name);
    bigint m = ((bigint)n1) * n2 * n3;
    for (bigint i = 0; i < m; i++)
      array[0][0][i] -= n4lo;
    m = ((bigint)n1) * n2;
    for (bigint i = 0; i < m; i++)
      array[0][i] -= n3lo;
    for (int i = 0; i < n1; i++)
      array[i] -= n2lo;
    return array;
  }
  template <typename TYPE>
  TYPE ****create4d_offset(TYPE *****&, int, int, int, int, int, int, int,
                           const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  void destroy4d_offset(TYPE ****&array, int n2_offset, int n3_offset,
                        int n4_offset) {
    if (array == nullptr)
      return;
    sfree(&array[0][n2_offset][n3_offset][n4_offset]);
    sfree(&array[0][n2_offset][n3_offset]);
    sfree(&array[0][n2_offset]);
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE ****create4d_offset_last(TYPE ****&array, int n1lo, int n1hi, int n2lo,
                                int n2hi, int n3lo, int n3hi, int n4,
                                const char *name) {
    if (n1lo > n1hi || n2lo > n2hi || n3lo > n3hi || n4 <= 0)
      return nullptr;
    int n1 = n1hi - n1lo + 1;
    int n2 = n2hi - n2lo + 1;
    int n3 = n3hi - n3lo + 1;
    create(array, n1, n2, n3, n4, name);
    bigint m = ((bigint)n1) * n2;
    for (bigint i = 0; i < m; i++)
      array[0][i] -= n3lo;
    for (int i = 0; i < n1; i++)
      array[i] -= n2lo;
    array -= n1lo;
    return array;
  }
  template <typename TYPE>
  TYPE ****create4d_offset_last(TYPE *****&, int, int, int, int, int, int, int,
                                const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE>
  void destroy4d_offset_last(TYPE ****&array, int n1_offset, int n2_offset,
                             int n3_offset) {
    if (array == nullptr)
      return;
    sfree(&array[n1_offset][n2_offset][n3_offset][0]);
    sfree(&array[n1_offset][n2_offset][n3_offset]);
    sfree(&array[n1_offset][n2_offset]);
    sfree(&array[n1_offset]);
    array = nullptr;
  }
  template <typename TYPE>
  TYPE *****create(TYPE *****&array, int n1, int n2, int n3, int n4, int n5,
                   const char *name) {
    bigint nbytes = ((bigint)sizeof(TYPE)) * n1 * n2 * n3 * n4 * n5;
    TYPE *data = (TYPE *)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE *)) * n1 * n2 * n3 * n4;
    TYPE **level4 = (TYPE **)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE **)) * n1 * n2 * n3;
    TYPE ***level3 = (TYPE ***)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE ***)) * n1 * n2;
    TYPE ****level2 = (TYPE ****)smalloc(nbytes, name);
    nbytes = ((bigint)sizeof(TYPE ****)) * n1;
    array = (TYPE *****)smalloc(nbytes, name);
    int i, j, k, l;
    bigint m1, m2;
    bigint n = 0;
    for (i = 0; i < n1; i++) {
      m2 = ((bigint)i) * n2;
      array[i] = &level2[m2];
      for (j = 0; j < n2; j++) {
        m1 = ((bigint)i) * n2 + j;
        m2 = ((bigint)i) * n2 * n3 + ((bigint)j) * n3;
        level2[m1] = &level3[m2];
        for (k = 0; k < n3; k++) {
          m1 = ((bigint)i) * n2 * n3 + ((bigint)j) * n3 + k;
          m2 = ((bigint)i) * n2 * n3 * n4 + ((bigint)j) * n3 * n4 +
               ((bigint)k) * n4;
          level3[m1] = &level4[m2];
          for (l = 0; l < n4; l++) {
            m1 = ((bigint)i) * n2 * n3 * n4 + ((bigint)j) * n3 * n4 +
                 ((bigint)k) * n4 + l;
            level4[m1] = &data[n];
            n += n5;
          }
        }
      }
    }
    return array;
  }
  template <typename TYPE>
  TYPE ******create(TYPE ******&, int, int, int, int, int, const char *name) {
    fail(name);
    return nullptr;
  }
  template <typename TYPE> void destroy(TYPE *****&array) {
    if (array == nullptr)
      return;
    sfree(array[0][0][0][0]);
    sfree(array[0][0][0]);
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }
  template <typename TYPE> double usage(TYPE *array, int n) {
    (void)array;
    double bytes = ((double)sizeof(TYPE)) * n;
    return bytes;
  }
  template <typename TYPE> double usage(TYPE **array, int n1, int n2) {
    (void)array;
    double bytes = ((double)sizeof(TYPE)) * n1 * n2;
    bytes += ((double)sizeof(TYPE *)) * n1;
    return bytes;
  }
  template <typename TYPE> double usage(TYPE ***array, int n1, int n2, int n3) {
    (void)array;
    double bytes = ((double)sizeof(TYPE)) * n1 * n2 * n3;
    bytes += ((double)sizeof(TYPE *)) * n1 * n2;
    bytes += ((double)sizeof(TYPE **)) * n1;
    return bytes;
  }
  template <typename TYPE>
  double usage(TYPE ****array, int n1, int n2, int n3, int n4) {
    (void)array;
    double bytes = ((double)sizeof(TYPE)) * n1 * n2 * n3 * n4;
    bytes += ((double)sizeof(TYPE *)) * n1 * n2 * n3;
    bytes += ((double)sizeof(TYPE **)) * n1 * n2;
    bytes += ((double)sizeof(TYPE ***)) * n1;
    return bytes;
  }
};
} // namespace LAMMPS_NS
#endif
