# pow
replace golang math.Pow to optimal

#### godoc

```
FUNCTIONS

func E2(x float64) float64
    E2 function replace `math.Pow(x,2.0)`

func E3(x float64) float64
    E3 function replace `math.Pow(x,3.0)`
```

#### benchmark

```
go test -bench=. -count=5 > bench.txt
benchstat bench.txt 
```
```
name          time/op
/math.Pow2-4   357ns ± 0%
/pow.E2-4     5.69ns ± 0%

/math.Pow3-4   359ns ± 1%
/pow.E3-4     8.86ns ± 1%
```

So, that package is more optimal at more 30 times.

#### example

```golang
func Example() {
	x := 2.0
	r2 := pow.E2(x) // math.Pow(x, 2.0)
	r3 := pow.E3(x) // math.Pow(x, 3.0)
	fmt.Fprintf(os.Stdout, "%.4f\n", r2)
	fmt.Fprintf(os.Stdout, "%.4f\n", r3)
	// Output:
	// 4.0000
	// 8.0000
}
```
