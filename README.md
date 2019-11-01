[![Coverage Status](https://coveralls.io/repos/github/Konstantin8105/hd/badge.svg?branch=master)](https://coveralls.io/github/Konstantin8105/hd?branch=master)
[![Build Status](https://travis-ci.org/Konstantin8105/hd.svg?branch=master)](https://travis-ci.org/Konstantin8105/hd)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/hd)](https://goreportcard.com/report/github.com/Konstantin8105/hd)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/hd?status.svg)](https://godoc.org/github.com/Konstantin8105/hd)
![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)

# hd

FEM(finite element method) for structural engineer

```
Specific of truss:
* Generally truss finite element and beam finite element with moment free on node is not the same.
for truss stiffiner matrix look like that
⎡ * 0 0 * 0 0 ⎤
⎢ 0 0 0 0 0 0 ⎥
⎢ 0 0 0 0 0 0 ⎥
⎢ * 0 0 * 0 0 ⎥
⎢ 0 0 0 0 0 0 ⎥
⎣ 0 0 0 0 0 0 ⎦
for beam with moment free on nodes stiffiner matrix look like that
⎡ * 0 0 * 0 0 ⎤
⎢ 0 * 0 0 * 0 ⎥
⎢ 0 0 0 0 0 0 ⎥
⎢ * 0 0 * 0 0 ⎥
⎢ 0 * 0 0 * 0 ⎥
⎣ 0 0 0 0 0 0 ⎦
So, if we have 2 free moment, then also add free by Y direction for bith beam points
```

#### Coverage

```
go test -coverprofile=coverage.out ./...
go tool cover -html=coverage.out
```

or on one line:

```
go test -coverprofile=coverage.out ./... ; go tool cover -html=coverage.out ; rm coverage.out
```

### Benchmark

In folder `mod` run:

```cmd
go test -v -run=Benchmark -bench=Benchmark -benchmem -cpuprofile cpu.prof -memprofile mem.prof
```

```
goos: linux
goarch: amd64
pkg: github.com/Konstantin8105/hd/mod
BenchmarkRun/____1-cases20-4         	    5184	    202061 ns/op	   95197 B/op	    1052 allocs/op
BenchmarkRun/____2-cases20-4         	    3676	    275535 ns/op	  134711 B/op	    1288 allocs/op
BenchmarkRun/____4-cases20-4         	    2658	    446752 ns/op	  226851 B/op	    1778 allocs/op
BenchmarkRun/____8-cases20-4         	    1460	    807288 ns/op	  440599 B/op	    2738 allocs/op
BenchmarkRun/___16-cases20-4         	     687	   1736805 ns/op	  928967 B/op	    4640 allocs/op
BenchmarkRun/___32-cases20-4         	     230	   5137632 ns/op	 2239788 B/op	    8456 allocs/op
BenchmarkRun/___64-cases20-4         	      57	  18934656 ns/op	 6173719 B/op	   16063 allocs/op
BenchmarkRun/__128-cases20-4         	      10	 105636607 ns/op	18726030 B/op	   31323 allocs/op
PASS
ok  	github.com/Konstantin8105/hd/mod	11.064s
```

