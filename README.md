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
