package verif

import (
	"fmt"
	"os"
	"strings"

	"github.com/Konstantin8105/hd"
)

func Example() {
	tcs := []func() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) (tol []float64)){
		G,
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
	// Book: William McGuire, Richard H.Gallagher, Ronald D.Ziemian Matrix Structural Analysis EXAMPLE 9.1 Page 247
	// case  0:    5.05%
	// case  1:    5.05%
}
