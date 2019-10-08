package pow

// E2 function replace `math.Pow(x,2.0)`
func E2(x float64) float64 {
	return x * x
}

// E3 function replace `math.Pow(x,3.0)`
func E3(x float64) float64 {
	return x * x * x
}

// En function replace `math.Pow(x,e)`
func En(x float64, e int) float64 {
	var isNegative bool
	if e < 0 {
		isNegative = true
		e = -e
	}
	var (
		i uint8
		n uint8   = uint8(e)
		r float64 = 1.0

		// variables for iterations
		tmp float64
		c   uint8
	)

	for n > 1 {
		tmp = x
		c = 1
		for c <= n/2 {
			tmp *= tmp
			c *= 2
		}
		n = n - c
		r *= tmp
	}

	for i = 0; i < n; i++ {
		r *= x
	}
	if isNegative {
		return 1 / r
	}
	return r
}
