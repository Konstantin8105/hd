package verif

import (
	"fmt"
	"os"
	"strings"

	"github.com/Konstantin8105/hd"
)

func Example() {
	tcs := []func() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) (tol []float64)){
		MSA21,
		MSA49,
		MSA413,
		MSA91,
	}

	for _, tc := range tcs {
		model, lc, name, isOk := tc()
		name = strings.ReplaceAll(name, "\n", " ")
		name = strings.ReplaceAll(name, "  ", " ")
		name = strings.TrimSpace(name)
		fmt.Fprintf(os.Stdout, "%s\n", name)
		err := hd.LinearStatic(nil, &model, &lc)
		if err != nil {
			fmt.Fprintf(os.Stdout, "%v\n", err)
			continue
		}
		tols := isOk(&lc)
		for i := range tols {
			fmt.Fprintf(os.Stdout, "case %2d: %7.2f%%\n", i, tols[i])
		}
	}
	// Output:
	// Book: William McGuire, Richard H.Gallagher, Ronald D.Ziemian Matrix Structural Analysis EXAMPLE 2.1 Page 21
	// case  0:    0.05%
	// case  1:    0.46%
	// case  2:    0.00%
	// case  3:   -0.01%
	// Book: William McGuire, Richard H.Gallagher, Ronald D.Ziemian Matrix Structural Analysis EXAMPLE 4.9 Page 80
	// case  0:   -1.79%
	// case  1:    0.44%
	// case  2:    0.44%
	// case  3:    0.00%
	// case  4:    0.00%
	// case  5:    0.06%
	// Book: William McGuire, Richard H.Gallagher, Ronald D.Ziemian Matrix Structural Analysis EXAMPLE 4.13 Page 84
	// case  0:    0.02%
	// case  1:    0.02%
	// case  2:    0.26%
	// Book: William McGuire, Richard H.Gallagher, Ronald D.Ziemian Matrix Structural Analysis EXAMPLE 9.1 Page 247
	// case  0:    5.05%
	// case  1:    5.05%
}
