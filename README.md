```
make -j 2
```

coverage

```
make 'CXXFLAGS = -coverage -O0 -g3' 'LDFLAGS = -coverage' -j 2
mpiexec -n 2 ./main -i in.run
python -m gcovr --html-details cover.html
```