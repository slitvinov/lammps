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
  template <typename TYPE> void destroy(TYPE ***&array) {
    if (array == nullptr)
      return;
    sfree(array[0][0]);
    sfree(array[0]);
    sfree(array);
    array = nullptr;
  }
};
} // namespace LAMMPS_NS
#endif
