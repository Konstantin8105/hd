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

	// TODO check loads, load cases

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
			if b.N[j] < 0 {
				et.Add(ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("negative index of point"),
				})
			}
			if b.N[j] >= len(m.Points) {
				et.Add(ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("outside index of point"),
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
			{b.G, "shear modulus"},
			{b.Density, "density"},
			{b.Cte, "coefficient of thermal expansion"},
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
