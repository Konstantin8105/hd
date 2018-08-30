package hd

import (
	"fmt"
	"math"
)

func (m *Model) checkInputData() (err error) {
	// points
	err = m.checkPoints()
	if err != nil {
		return err
	}
	// beams
	err = m.checkBeams()
	if err != nil {
		return err
	}
	// supports
	err = m.checkSupports()
	if err != nil {
		return err
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
		Err:     fmt.Errorf("Not acceptable NaN float"),
	}
}

func isInf(f float64) ErrorFunc {
	return ErrorFunc{
		isError: math.IsInf(f, 0),
		Err:     fmt.Errorf("Not acceptable infinity float"),
	}
}

func isPositive(f float64) ErrorFunc {
	return ErrorFunc{
		isError: f < 0,
		Err:     fmt.Errorf("Not acceptable negative float"),
	}
}

func isNotZero(f float64) ErrorFunc {
	return ErrorFunc{
		isError: f == 0,
		Err:     fmt.Errorf("Not acceptable zero float"),
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

func isOk(errs []ErrorFunc) (err error) {
	for _, e := range errs {
		if e.isError {
			return e.Err
		}
	}
	return nil
}

func (m *Model) checkPoints() (err error) {
	for i := range m.Points {
		for j := 0; j < len(m.Points[i]); j++ {
			err := isOk([]ErrorFunc{
				isNaN(m.Points[i][j]),
				isInf(m.Points[i][j]),
			})
			if err != nil {
				return ErrorPoint{
					PointIndex: i,
					CoordIndex: j,
					Err:        err,
				}
			}
		}
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

func (m *Model) checkBeams() (err error) {
	for i, b := range m.Beams {
		for j := 0; j < len(b.N); j++ {
			if b.N[j] < 0 {
				return ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("Not acceptable negative index of point"),
				}
			}
			if b.N[j] >= len(m.Points) {
				return ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("Not acceptable outside index of point"),
				}
			}
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
			err = isOk([]ErrorFunc{
				isNaN(value.v),
				isInf(value.v),
				isPositive(value.v),
				isNotZero(value.v),
			})
			if err != nil {
				return ErrorBeam{
					BeamIndex: i,
					Detail:    value.name,
					Err:       err,
				}
			}
		}
	}
	return nil
}

func (m *Model) checkSupports() (err error) {
	if len(m.Points) != len(m.Supports) {
		return fmt.Errorf("Amount of supports is not same amount of points")
	}
	return nil
}
