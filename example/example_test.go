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
		m        hd.Model
		filename string
	}{{
		m:        example.ConsoleBeam(),
		filename: "beam",
	}, {
		m:        example.Truss(),
		filename: "truss",
	}, {
		m:        example.GBeam(),
		filename: "g-beam",
	}, {
		m:        example.DoubleBeam(),
		filename: "double-beam",
	}, {
		m:        example.ModalTruss(),
		filename: "truss-modal",
	}, {
		m:        example.ModalTrussRotate(),
		filename: "truss-modal-rotate",
	}, {
		m:        example.ModalBeam(),
		filename: "beam-modal",
	}, {
		m:        example.ModalBeamRotate(),
		filename: "beam-modal-rotate",
	}, {
		m:        example.ModalBeam3mass(),
		filename: "beam-modal-3mass",
	}, {
		m:        example.BucklingBeam(),
		filename: "beam-buckling",
	}, {
		m: func() hd.Model {
			m, _ := example.BeamDc()
			return m
		}(),
		filename: "beam-dc-part1",
	}, {
		m: func() hd.Model {
			_, m := example.BeamDc()
			return m
		}(),
		filename: "beam-dc-part2",
	}}
	for _, m := range ms {
		t.Run(m.filename, func(t *testing.T) {
			var b bytes.Buffer
			if err := m.m.Run(&b); err != nil {
				t.Errorf("Cannot calculate : %v", err)
			}
			b.Reset()

			// compare files
			actual := []byte(m.m.String())

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
