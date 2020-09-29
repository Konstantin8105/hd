# pow
replace golang math.Pow to optimal

#### godoc

```
FUNCTIONS

func E2(x float64) float64
    E2 function replace `math.Pow(x,2.0)`

func E3(x float64) float64
    E3 function replace `math.Pow(x,3.0)`

func E4(x float64) float64
    E4 function replace `math.Pow(x,4.0)`

func En(x float64, e int) float64
    En function replace `math.Pow(x,e)`
```

#### benchmark

```
go test -bench=. -count=5 > bench.txt
benchstat bench.txt 
```
```
name                time/op
/math.Pow2-4        41.7ns ±11%
/pow.E2-4           0.60ns ± 9%
/pow.En2-4          6.88ns ± 1%

/math.Pow3-4        39.1ns ± 1%
/pow.E3-4           0.59ns ± 2%
/pow.En3-4          7.80ns ± 2%

/math.Pow4-4        42.8ns ± 2%
/pow.E4-4           0.58ns ± 0%
/pow.En4-4          7.76ns ± 3%

/math.Pow(x,_51)-4  46.3ns ± 0%
/pow.En(x,_51)-4    16.8ns ± 0%
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
