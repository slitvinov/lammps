#ifndef LMP_LAMMPSPLUGIN_H
#define LMP_LAMMPSPLUGIN_H 
#ifdef __cplusplus
extern "C" {
#endif
typedef void *(lammpsplugin_factory1) (void *);
typedef void *(lammpsplugin_factory2) (void *, int, char **);
typedef struct {
  const char *version;
  const char *style;
  const char *name;
  const char *info;
  const char *author;
  union {
    lammpsplugin_factory1 *v1;
    lammpsplugin_factory2 *v2;
  } creator;
  void *handle;
} lammpsplugin_t;
typedef void (*lammpsplugin_regfunc)(lammpsplugin_t *, void *);
typedef void (*lammpsplugin_initfunc)(void *, void *, void *);
void lammpsplugin_init(void *, void *, void *);
#ifdef __cplusplus
}
#endif
#endif
