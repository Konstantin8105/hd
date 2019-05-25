# sparse

[![codecov](https://codecov.io/gh/Konstantin8105/sparse/branch/master/graph/badge.svg)](https://codecov.io/gh/Konstantin8105/sparse)
[![Build Status](https://travis-ci.org/Konstantin8105/sparse.svg?branch=master)](https://travis-ci.org/Konstantin8105/sparse)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/sparse)](https://goreportcard.com/report/github.com/Konstantin8105/sparse)
[![GitHub license](https://img.shields.io/badge/license-LGPL%20v2.1-blue.svg)](https://github.com/Konstantin8105/sparse/blob/master/LICENSE)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/sparse?status.svg)](https://godoc.org/github.com/Konstantin8105/sparse)

**This package based on program CSparse from [SuiteSparse 5.3.0](http://faculty.cse.tamu.edu/davis/SuiteSparse/)**

### Comparation with `gonum.mat.LU`

```cmd
go test -v -cpuprofile cpu.prof -memprofile mem.prof -coverprofile=coverage.out -run=BenchmarkLU -bench=BenchmarkLU -benchmem
# goos: linux
# goarch: amd64
# pkg: github.com/Konstantin8105/sparse
# BenchmarkLU/Sparse:__30:__300-6         	   50000	     34007 ns/op	   51456 B/op	      32 allocs/op
# BenchmarkLU/Dense_:__30:__300-6         	   50000	     30012 ns/op	    9249 B/op	      11 allocs/op
# BenchmarkLU/Sparse:_100:_1070-6         	   10000	    125092 ns/op	  201458 B/op	      32 allocs/op
# BenchmarkLU/Dense_:_100:_1070-6         	    5000	    241981 ns/op	   85261 B/op	      11 allocs/op
# BenchmarkLU/Sparse:_300:_3270-6         	    5000	    354187 ns/op	  566515 B/op	      32 allocs/op
# BenchmarkLU/Dense_:_300:_3270-6         	     500	   2377480 ns/op	  746151 B/op	      22 allocs/op
# BenchmarkLU/Sparse:1000:10970-6         	    1000	   1152972 ns/op	 1794937 B/op	      32 allocs/op
# BenchmarkLU/Dense_:1000:10970-6         	      20	  86145264 ns/op	 8079531 B/op	      54 allocs/op
# BenchmarkLU/Sparse:3000:32970-6         	     500	   3320653 ns/op	 5383122 B/op	      32 allocs/op
# BenchmarkLU/Dense_:3000:32970-6         	       1	1328677635 ns/op	72346808 B/op	     156 allocs/op
```

Using sparse algorithm is effective in test case for square matrixes with size more 100.

### Example of comparing CSparse and CXSparse

```cmd
rm -rf /tmp/CXSparse
cp -R ./CXSparse/ /tmp/
sed -i 's/CS_ENTRY/double/g' /tmp/CXSparse/Source/*.c
sed -i 's/CS_INT/csi/g'   /tmp/CXSparse/Source/*.c
meld  /tmp/CXSparse/ ./CSparse/
```

### How updating package

* Check new version of `CSparse` is exist on [page](http://faculty.cse.tamu.edu/davis/SuiteSparse/).
* Download new version.
* Compare file-by-file with file in folder `CSparse` of that package.
* Inject changes into Go files of package.

> Note:
> CSparse haven`t any updates at the last few years, so
> that package is actual at the future.
>

### Just for information

**This package transpiled CSparse from C to Go by [c4go](https://github.com/Konstantin8105/c4go).**

### Profiling

```
go test -v -cpuprofile=cpu.prof -memprofile=mem.prof -run=Benchmark -bench=Benchmark -benchmem
go tool pprof cpu.prof
go tool pprof mem.prof
```

### Coverage

```
go test -coverprofile=coverage.out -run=TestLU
go tool cover -html=coverage.out
```

### Questions for CSparse

* Variables `css.lnz` and `css.unz` is not float type `double`, better to use integer type like `int`.
* 
