# pow
replace golang math.Pow to optimal

#### godoc

```
FUNCTIONS

func E2(x float64) float64
    E2 function replace `math.Pow(x,2.0)`

func E3(x float64) float64
    E3 function replace `math.Pow(x,3.0)`

func En(x float64, e int) float64
    En function replace `math.Pow(x,e)`
```

#### benchmark

```
go test -bench=. -count=5 > bench.txt
benchstat bench.txt 
```
```
name              time/op
/math.Pow2-5      41.1ns ± 3%
/pow.E2-5         0.60ns ± 1%
/pow.En2-5        7.10ns ± 0%

/math.Pow3-5      39.7ns ± 1%
/pow.E3-5         0.60ns ± 1%
/pow.En3-5        7.73ns ± 1%

/math.Pow(151)-5  41.0ns ± 0%
/pow.En(151)-5    17.6ns ± 1%
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
