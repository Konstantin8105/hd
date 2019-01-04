package mod

import (
	"bytes"
	"fmt"
	"math"
	"testing"

	"github.com/Konstantin8105/hd"
	"github.com/Konstantin8105/hd/example"
)

func TestSplitFail(t *testing.T) {
	m := example.ConsoleBeam()
	tcs := []struct {
		beamIndex, amounts int
	}{
		{-1, 5},
		{100, 5},
		{0, -1},
	}
	for i, tc := range tcs {
		t.Run(fmt.Sprintf("Case #%d", i), func(t *testing.T) {
			if err := SplitBeam(&m, tc.beamIndex, tc.amounts); err == nil {
				t.Errorf("Error for %v", tc)
			}
		})
	}
}

func TestSplit(t *testing.T) {
	models := []func() hd.Model{
		example.ConsoleBeam, example.DoubleBeam,
		example.GBeam,
		example.Truss,
		example.ModalBeam, example.ModalBeamRotate,
		example.ModalTruss, example.ModalTrussRotate,
	}
	for mIndex := range models {
		t.Run(fmt.Sprintf("Model%d", mIndex), func(t *testing.T) {
			m := models[mIndex]()
			var b bytes.Buffer
			if err := m.Run(&b); err != nil {
				t.Fatalf("Error : %v", err)
			}
			expectResult := m.LoadCases
			hz := m.ModalCases[0].Result[0].Hz

			for i := 1; i < 10; i++ {
				t.Run(fmt.Sprintf("Split%d", i), func(t *testing.T) {
					mLocal := models[mIndex]()

					// split each beams
					var b bytes.Buffer
					amountBeams := len(mLocal.Beams)
					for j := 0; j < amountBeams; j++ {
						if err := SplitBeam(&mLocal, j, i); err != nil {
							t.Fatalf("Cannot split %d: %v", i, err)
						}
					}

					// calculation
					if err := mLocal.Run(&b); err != nil {
						t.Fatalf("Error : %v", err)
					}

					// eps
					eps := 1e-9

					if err := Compare(expectResult, mLocal.LoadCases, eps); err != nil {
						t.Log(mLocal.String())
						t.Errorf("Result is not same: %v", err)
					}

					h := mLocal.ModalCases[0].Result[0].Hz
					diff := math.Abs((hz - h) / h)
					if diff > eps {
						t.Logf("Narural frequency: %15.5e != %15.5e", hz, h)
						t.Errorf("Diff in natural frequency is not ok : %15.5e", diff)
					}
				})
			}
		})
	}
}

func BenchmarkRun(b *testing.B) {
	minimalLoadCases := 20
	for ic := 1; ic <= 128; ic *= 2 {
		b.Run(fmt.Sprintf("%5d-cases%d", ic, minimalLoadCases), func(b *testing.B) {
			// prepare model
			m := example.ConsoleBeam()
			// add more finite elements
			err := SplitBeam(&m, 0, ic)
			if err != nil {
				panic(err)
			}
			// add more cases
			for i := 0; len(m.LoadCases) < minimalLoadCases; i++ {
				lc := m.LoadCases[0]
				lc.LoadNodes[0].Forces[1] += float64(i)
				m.LoadCases = append(m.LoadCases, lc)
			}
			for i := 2; i < len(m.Points); i += 2 {
				m.ModalCases[0].ModalMasses = append(m.ModalCases[0].ModalMasses,
					hd.ModalMass{N: i, Mass: 100})
			}
			var bb bytes.Buffer
			b.ResetTimer()
			// run benchmark
			for i := 0; i < b.N; i++ {
				err := m.Run(&bb)
				if err != nil {
					panic(err)
				}
				bb.Reset()
			}
		})
	}
}
