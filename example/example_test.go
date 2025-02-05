package example_test

import (
	"bytes"
	"testing"

	"github.com/Konstantin8105/compare"
	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/hd/example"
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
		mf:       example.TrussWithBuckling,
		filename: "truss-with-buckling",
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
			name := "./testdata/" + m.filename
			compare.Test(t, name, b.Bytes())
		})
	}
}
