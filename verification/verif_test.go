package verif

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"testing"

	"github.com/Konstantin8105/hd"
	"github.com/pmezard/go-difflib/difflib"
)

func Test(t *testing.T) {
	tcs := []struct {
		name string
		f    func() (model hd.Model, lc hd.LoadCase, name string, isOk func(lc *hd.LoadCase) string)
	}{
		{name: "MSA21", f: MSA21},
		{name: "MSA49", f: MSA49},
		{name: "MSA413", f: MSA413},
		{name: "MSA67", f: MSA67},
		{name: "MSA81", f: MSA81},
		{name: "MSA91", f: MSA91},
	}

	for _, tc := range tcs {
		t.Run(tc.name, func(t *testing.T) {
			var actual bytes.Buffer
			{
				model, lc, name, isOk := tc.f()
				fmt.Fprintf(&actual, "name : %s\n", name)
				err := hd.LinearStatic(nil, &model, &lc)
				if err != nil {
					t.Fatalf("%v\n%v\n", err, lc)
				}
				fmt.Fprintf(&actual, "%s\n", isOk(&lc))
			}

			// compare files
			if os.Getenv("UPDATE") != "" {
				err := ioutil.WriteFile("."+tc.name, actual.Bytes(), 0644)
				if err != nil {
					t.Fatalf("Cannot Update: %v", err)
				}
			}

			expect, err := ioutil.ReadFile("." + tc.name)
			if err != nil {
				t.Fatalf("Cannot read file : %v", err)
			}

			if !bytes.Equal(expect, actual.Bytes()) {
				// show a diff between files
				diff := difflib.UnifiedDiff{
					A:        difflib.SplitLines(string(expect)),
					B:        difflib.SplitLines(actual.String()),
					FromFile: "Original",
					ToFile:   "Current",
					Context:  30000,
				}
				text, _ := difflib.GetUnifiedDiffString(diff)
				t.Log(text)
				t.Errorf("result is not same")
			}
		})

	}
}
