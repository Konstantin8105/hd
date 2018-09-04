package hd

import (
	"fmt"
	"math"
)

func (m *Model) checkInputData() error {
	et := ErrorTree{Name: "checkInputData"}
	// points
	if err := m.checkPoints(); err != nil {
		et.Add(err)
	}
	// beams
	if err := m.checkBeams(); err != nil {
		et.Add(err)
	}
	// supports
	if err := m.checkSupports(); err != nil {
		et.Add(err)
	}
	// load cases
	if err := m.checkLoadCases(); err != nil {
		et.Add(err)
	}
	// load
	if err := m.checkLoad(); err != nil {
		et.Add(err)
	}
	// error handling
	if et.IsError() {
		return et
	}
	return nil
}

type ErrorFunc struct {
	isError bool
	Err     error
}

func isNaN(f float64) ErrorFunc {
	return ErrorFunc{
		isError: math.IsNaN(f),
		Err:     fmt.Errorf("NaN float"),
	}
}

func isInf(f float64) ErrorFunc {
	return ErrorFunc{
		isError: math.IsInf(f, 0),
		Err:     fmt.Errorf("infinity float"),
	}
}

func isPositive(f float64) ErrorFunc {
	return ErrorFunc{
		isError: f < 0,
		Err:     fmt.Errorf("negative float"),
	}
}

func isTrue(b bool) ErrorFunc {
	return ErrorFunc{
		isError: b,
		Err:     fmt.Errorf("error case is true"),
	}
}

func isNotZero(f float64) ErrorFunc {
	return ErrorFunc{
		isError: f == 0,
		Err:     fmt.Errorf("zero float"),
	}
}

type ErrorPoint struct {
	PointIndex int
	CoordIndex int
	Err        error
}

func (e ErrorPoint) Error() string {
	return fmt.Sprintf("Error in point %d, coordinate %d : %v",
		e.PointIndex,
		e.CoordIndex,
		e.Err)
}

func isOk(errs ...ErrorFunc) (err error) {
	for _, e := range errs {
		if e.isError {
			return e.Err
		}
	}
	return nil
}

func (m *Model) checkPoints() error {
	et := ErrorTree{Name: "checkPoints"}
	for i := range m.Points {
		for j := 0; j < len(m.Points[i]); j++ {
			err := isOk(
				isNaN(m.Points[i][j]),
				isInf(m.Points[i][j]),
			)
			if err != nil {
				et.Add(ErrorPoint{
					PointIndex: i,
					CoordIndex: j,
					Err:        err,
				})
			}
		}
	}
	if et.IsError() {
		return et
	}
	return nil
}

type ErrorBeam struct {
	BeamIndex int
	Detail    string
	Err       error
}

func (e ErrorBeam) Error() string {
	return fmt.Sprintf("Error in beam %d, checking `%s` : %v",
		e.BeamIndex,
		e.Detail,
		e.Err)
}

func (m *Model) checkBeams() error {
	et := ErrorTree{Name: "checkBeams"}
	for i, b := range m.Beams {
		for j := 0; j < len(b.N); j++ {
			err := isOk(
				isTrue(b.N[j] < 0),
				isTrue(b.N[j] >= len(m.Points)),
			)
			if err != nil {
				et.Add(ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("outside index of point : %d", b.N[j]),
				})
			}
		}

		if m.distance(b.N[0], b.N[1]) <= 0 {
			et.Add(ErrorBeam{
				BeamIndex: i,
				Detail:    "distance between points",
				Err:       fmt.Errorf("is less or equal zero"),
			})
		}

		values := []struct {
			v    float64
			name string
		}{
			{b.A, "cross-section area"},
			{b.J, "moment of inertia"},
			{b.E, "modulus of elasticity"},
		}

		for _, value := range values {
			err := isOk(
				isNaN(value.v),
				isInf(value.v),
				isPositive(value.v),
				isNotZero(value.v),
			)
			if err != nil {
				et.Add(ErrorBeam{
					BeamIndex: i,
					Detail:    value.name,
					Err:       err,
				})
			}
		}
	}
	if et.IsError() {
		return et
	}
	return nil
}

func (m *Model) checkSupports() (err error) {
	if len(m.Points) != len(m.Supports) {
		return fmt.Errorf("Amount of supports is not same amount of points")
	}
	return nil
}

func (m *Model) checkLoadCases() (err error) {
	// always ok
	return nil
}

type ErrorLoad struct {
	LoadCase int
	LoadPos  int
	Err      error
}

func (e ErrorLoad) Error() string {
	return fmt.Sprintf("Error in load case %d, load position %d : %v",
		e.LoadCase,
		e.LoadPos,
		e.Err)
}

func (m *Model) checkLoad() (err error) {
	et := ErrorTree{Name: "checkLoad"}
	for i := range m.LoadCases {
		for j := range m.LoadCases[i].LoadNodes {
			ld := m.LoadCases[i].LoadNodes[j]
			err := isOk(
				isTrue(ld.N < 0),
				isTrue(ld.N >= len(m.Points)),
			)
			if err != nil {
				et.Add(ErrorLoad{
					LoadCase: i,
					LoadPos:  j,
					Err:      fmt.Errorf("outside index of point : %d", ld.N),
				})
			}
			for k := 0; k < 3; k++ {
				err := isOk(
					isNaN(ld.Forces[k]),
					isInf(ld.Forces[k]),
				)
				if err != nil {
					et.Add(ErrorLoad{
						LoadCase: i,
						LoadPos:  j,
						Err:      err,
					})
				}
			}
		}
	}
	if et.IsError() {
		return et
	}
	return nil
}
