package example_test

import (
	"bytes"
	"io/ioutil"
	"os"
	"testing"

	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/hd/example"
	"github.com/pmezard/go-difflib/difflib"
)

func TestModelString(t *testing.T) {
	ms := []struct {
		mf       func() (hd.Model, []hd.LoadCase, []hd.ModalCase)
		filename string
	}{{
		mf:       example.ConsoleBeam,
		filename: "beam",
	}, {
		mf:       example.Truss,
		filename: "truss",
	}, {
		mf:       example.GBeam,
		filename: "g-beam",
	}, {
		mf:       example.DoubleBeam,
		filename: "double-beam",
	}, {
		mf:       example.ModalTruss,
		filename: "truss-modal",
	}, {
		mf:       example.ModalTrussRotate,
		filename: "truss-modal-rotate",
	}, {
		mf:       example.ModalBeam,
		filename: "beam-modal",
	}, {
		mf:       example.ModalBeamRotate,
		filename: "beam-modal-rotate",
	}, {
		mf:       example.ModalBeam3mass,
		filename: "beam-modal-3mass",
	}, {
		mf:       example.BucklingBeam,
		filename: "beam-buckling",
	}, {
		mf: func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
			m, _, lc, _ := example.BeamDc()
			return func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
				return m, lc, nil
			}()
		},
		filename: "beam-dc-part1",
	}, {
		mf: func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
			_, m, _, lc := example.BeamDc()
			return func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
				return m, lc, nil
			}()
		},
		filename: "beam-dc-part2",
	}}
	for _, m := range ms {
		t.Run(m.filename, func(t *testing.T) {
			model, lcs, mcs := m.mf()
			var b bytes.Buffer
			if err := hd.Run(&b, &model, lcs, mcs); err != nil {
				t.Errorf("Cannot calculate : %v", err)
			}
			b.Reset()

			// compare files
			// TODO : String() for local case and modal in not added
			actual := []byte(model.String())

			if os.Getenv("UPDATE") != "" {
				err := ioutil.WriteFile("./testdata/"+m.filename, actual, 0644)
				if err != nil {
					t.Fatalf("Cannot Update: %v", err)
				}
			}

			expect, err := ioutil.ReadFile("./testdata/" + m.filename)
			if err != nil {
				t.Fatalf("Cannot read file : %v", err)
			}

			if !bytes.Equal(expect, actual) {
				// show a diff between files
				diff := difflib.UnifiedDiff{
					A:        difflib.SplitLines(string(expect)),
					B:        difflib.SplitLines(string(actual)),
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
