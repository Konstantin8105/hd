package hd

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"
)

func (m *Model) checkInputData() error {
	et := errors.Tree{Name: "checkInputData"}
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
	// pins
	if err := m.checkPins(); err != nil {
		et.Add(err)
	}
	// load cases
	// no test cases

	// load
	if err := m.checkLoad(); err != nil {
		et.Add(err)
	}
	// modal cases
	if err := m.checkModalCases(); err != nil {
		et.Add(err)
	}
	// modal
	if err := m.checkModal(); err != nil {
		et.Add(err)
	}
	// error handling
	if et.IsError() {
		return et
	}
	return nil
}

type errorFunc struct {
	isError bool
	Err     error
}

func isNaN(f float64) errorFunc {
	return errorFunc{
		isError: math.IsNaN(f),
		Err:     fmt.Errorf("NaN float"),
	}
}

func isInf(f float64) errorFunc {
	return errorFunc{
		isError: math.IsInf(f, 0),
		Err:     fmt.Errorf("infinity float"),
	}
}

func isPositive(f float64) errorFunc {
	return errorFunc{
		isError: f < 0,
		Err:     fmt.Errorf("negative float"),
	}
}

func isTrue(b bool) errorFunc {
	return errorFunc{
		isError: b,
		Err:     fmt.Errorf("error case is true"),
	}
}

func isNotZero(f float64) errorFunc {
	return errorFunc{
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

func isOk(errs ...errorFunc) (err error) {
	for _, e := range errs {
		if e.isError {
			return e.Err
		}
	}
	return nil
}

func (m *Model) checkPoints() error {
	et := errors.Tree{Name: "checkPoints"}
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
	et := errors.Tree{Name: "checkBeams"}
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

type ErrorPin struct {
	Beam int
	Err  error
}

func (e ErrorPin) Error() string {
	return fmt.Sprintf("Error in pin of beam %d: %v",
		e.Beam,
		e.Err)
}

func (m *Model) checkPins() (err error) {
	if len(m.Pins) == 0 {
		return nil
	}
	et := errors.Tree{Name: "checkPins"}
	if len(m.Beams) != len(m.Pins) {
		et.Add(fmt.Errorf("Amount of pins is not same beams, %d != %d",
			len(m.Beams), len(m.Pins)))
	}
	for beam := range m.Pins {
		// maximal allowable pins
		maxPins := 4
		pins := 0
		for j := 0; j < 6; j++ {
			if m.Pins[beam][j] {
				pins++
			}
		}
		if pins > maxPins {
			et.Add(ErrorPin{
				Beam: beam,
				Err: fmt.Errorf("Amount of pins of point is %d more %d",
					pins, maxPins),
			})
		}
		// both X is free
		if m.Pins[beam][0] && m.Pins[beam][3] {
			et.Add(ErrorPin{
				Beam: beam,
				Err:  fmt.Errorf("Not acceptable X direction is free"),
			})
		}
	}
	if et.IsError() {
		return et
	}
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
	et := errors.Tree{Name: "checkLoad"}
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

func (m *Model) checkModalCases() (err error) {
	et := errors.Tree{Name: "checkModalCases"}
	// empty modal mass is not acceptable
	for i := range m.ModalCases {
		err := isOk(
			isTrue(len(m.ModalCases[i].ModalMasses) == 0),
		)
		if err != nil {
			et.Add(ErrorModal{
				ModalCase: i,
				ModalPos:  0,
				Err:       fmt.Errorf("Modal case haven`t masses"),
			})
		}
	}
	if et.IsError() {
		return et
	}

	return nil
}

type ErrorModal struct {
	ModalCase int
	ModalPos  int
	Err       error
}

func (e ErrorModal) Error() string {
	return fmt.Sprintf("Error in modal case %d, mass position %d : %v",
		e.ModalCase,
		e.ModalPos,
		e.Err)
}

func (m *Model) checkModal() (err error) {
	et := errors.Tree{Name: "checkModal"}
	for i := range m.ModalCases {
		for j := range m.ModalCases[i].ModalMasses {
			ld := m.ModalCases[i].ModalMasses[j]
			err := isOk(
				isTrue(ld.N < 0),
				isTrue(ld.N >= len(m.Points)),
			)
			if err != nil {
				et.Add(ErrorModal{
					ModalCase: i,
					ModalPos:  j,
					Err:       fmt.Errorf("outside index of point : %d", ld.N),
				})
			}
			err = isOk(
				isNaN(ld.Mass),
				isInf(ld.Mass),
				isPositive(ld.Mass),
				isNotZero(ld.Mass),
			)
			if err != nil {
				et.Add(ErrorModal{
					ModalCase: i,
					ModalPos:  j,
					Err:       err,
				})
			}
			for s := 0; s < len(m.Supports); s++ {
				if m.Supports[s][0] || m.Supports[s][1] || m.Supports[s][2] {
					if ld.N == s {
						et.Add(ErrorModal{
							ModalCase: i,
							ModalPos:  j,
							Err:       fmt.Errorf("Modal mass on support"),
						})
					}
				}
			}
		}
	}
	if et.IsError() {
		return et
	}
	return nil
}
