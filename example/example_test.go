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
		mf:       example.ModalBeam2,
		filename: "beam-modal2",
	}, {
		mf:       example.ModalBeamRotate2,
		filename: "beam-modal-rotate2",
	}, {
		mf:       example.ModalBeam,
		filename: "beam-modal",
	}, {
		mf:       example.ModalBeamRotate,
		filename: "beam-modal-rotate",
	}, {
		mf:       example.BeamWithBuckling,
		filename: "beam-with-buckling",
	}, {
		mf:       example.ModalBeam3mass,
		filename: "beam-modal-3mass",
	}, {
		mf:       example.FrameModal,
		filename: "beam-frame-modal",
	}, {
		mf:       example.BucklingBeam,
		filename: "beam-buckling",
	}, {
		mf:       example.Gframe,
		filename: "gframe",
	}, {
		mf: func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
			return example.G(true)
		},
		filename: "G-linear",
	}, {
		mf: func() (hd.Model, []hd.LoadCase, []hd.ModalCase) {
			return example.G(false)
		},
		filename: "G-nonlinear",
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
	}, {
		mf:       example.Pframe,
		filename: "pframe",
		// 	}, {
		// 		mf:       example.Frame,
		// 		filename: "frame",
	}, {
		mf:       example.EFESTS10bar,
		filename: "EFESTS10bar",
	},
	}
	for _, m := range ms {
		t.Run(m.filename, func(t *testing.T) {
			model, lcs, mcs := m.mf()
			var b bytes.Buffer
			if err := hd.Run(&b, &model, lcs, mcs); err != nil {
				t.Errorf("Cannot calculate : %v", err)
			}

			// compare files
			actual := []byte(b.String())
			b.Reset()

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
