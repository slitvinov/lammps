```
make
```

coverage

```
make 'CXXFLAGS = -coverage -O0 -g3' 'LDFLAGS = -coverage' -j `nproc` && mpiexec -n 2 ./main -i in.run && python -m gcovr --html-details cover.html
```
