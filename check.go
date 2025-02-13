package hd

import (
	"fmt"
	"math"

	goerrors "errors"

	"github.com/Konstantin8105/errors"
)

func (m *Model) checkInputData() error {
	et := errors.Tree{Name: "check input data"}
	// points
	if err := m.checkPoints(); err != nil {
		_ = et.Add(err)
	}
	// beams
	if err := m.checkBeams(); err != nil {
		_ = et.Add(err)
	}
	// supports
	if err := m.checkSupports(); err != nil {
		_ = et.Add(err)
	}
	// pins
	if err := m.checkPins(); err != nil {
		_ = et.Add(err)
	}

	if !et.IsError() {
		// graph checking
		if err := m.graphCheck(); err != nil {
			_ = et.Add(err)
		}
	}

	// error handling
	if et.IsError() {
		return et
	}
	return nil
}

// checking for calculate one structure on model.
// graph moving by beams with mark.
func (m *Model) graphCheck() (err error) {
	if len(m.Beams) == 0 {
		return
	}
	// all points mark are false
	mark := make([]bool, len(m.Points))
	beamMark := make([]bool, len(m.Beams))

	// walk by beam and mark all
	var presentBeam, nextBeam []int
	presentBeam = append(presentBeam, 0) // begin

	for iter := 0; iter < 10000000; iter++ {
		for _, pr := range presentBeam {
			// mark points
			mark[m.Beams[pr].N[0]] = true
			mark[m.Beams[pr].N[1]] = true
			// find near beams
			nodes := m.Beams[pr].N
			for index := range m.Beams {
				var add bool
				for i := range nodes {
					if m.Beams[index].N[i] == nodes[0] ||
						m.Beams[index].N[i] == nodes[1] {
						add = true
					}
				}
				if add {
					nextBeam = append(nextBeam, index)
				}
			}
		}
		for _, pr := range presentBeam {
			beamMark[pr] = true
		}
		presentBeam = presentBeam[:0]
		// remove marked beams
		for _, pr := range nextBeam {
			if beamMark[pr] == true {
				continue
			}
			presentBeam = append(presentBeam, pr)
		}
		// is ready to exit
		if len(presentBeam) == 0 {
			break
		}
	}

	// list of unmarked points
	var list []int
	for index := range mark {
		if !mark[index] {
			list = append(list, index)
		}
	}

	if len(list) == 0 {
		return nil
	}

	return ErrorAmountModel
}

// ErrorAmountModel error if calculation model is not single.
// Calculate calculation models separately.
var ErrorAmountModel = goerrors.New("Model isn`t single calculation model")

// ErrorLoad error in load data
type ErrorLoad struct {
	LoadPos int
	Err     error
}

func (e ErrorLoad) Error() string {
	return fmt.Sprintf("Error in load case, load position %d : %v",
		e.LoadPos,
		e.Err)
}

func (l LoadCase) checkInputData(m *Model) error {
	et := errors.Tree{Name: "check input data"}

	// load
	for j := range l.LoadNodes {
		ld := l.LoadNodes[j]
		err := isOk(
			isTrue(ld.N < 0),
			isTrue(ld.N >= len(m.Points)),
		)
		if err != nil {
			_ = et.Add(ErrorLoad{
				LoadPos: j,
				Err:     fmt.Errorf("outside index of point : %d", ld.N),
			})
		}
		for k := 0; k < 3; k++ {
			err := isOk(
				isNaN(ld.Forces[k]),
				isInf(ld.Forces[k]),
			)
			if err != nil {
				_ = et.Add(ErrorLoad{
					LoadPos: j,
					Err:     err,
				})
			}
		}
	}

	// error handling
	if et.IsError() {
		return et
	}

	return nil
}

func isNaN(f float64) error {
	if math.IsNaN(f) {
		return fmt.Errorf("NaN float")
	}
	return nil
}

func isInf(f float64) error {
	if math.IsInf(f, 0) {
		return fmt.Errorf("infinity float")
	}
	return nil
}

func isPositive(f float64) error {
	if f < 0 {
		return fmt.Errorf("negative float")
	}
	return nil
}

func isTrue(b bool) error {
	if b {
		return fmt.Errorf("error case is true")
	}
	return nil
}

func isNotZero(f float64) error {
	if f == 0 {
		return fmt.Errorf("zero float")
	}
	return nil
}

// ErrorPoint error in point data
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

func isOk(errs ...error) (err error) {
	for _, e := range errs {
		if e != nil {
			return e
		}
	}
	return nil
}

func (m *Model) checkPoints() error {
	et := errors.Tree{Name: "check points"}
	for i := range m.Points {
		for j := 0; j < len(m.Points[i]); j++ {
			err := isOk(
				isNaN(m.Points[i][j]),
				isInf(m.Points[i][j]),
			)
			if err != nil {
				_ = et.Add(ErrorPoint{
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

// ErrorBeam error in beam data
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

// ErrorBeamAmount error if calculation model have points, but haven`t beams.
var ErrorBeamAmount = goerrors.New("Zero beams with non-zero amount of points")

func (m *Model) checkBeams() error {
	et := errors.Tree{Name: "check beams"}
	// add len beam is zero but points more zero
	if len(m.Points) > 0 && len(m.Beams) == 0 {
		return ErrorBeamAmount
	}
	// check beam property
	for i, b := range m.Beams {
		for j := 0; j < len(b.N); j++ {
			err := isOk(
				isTrue(b.N[j] < 0),
				isTrue(b.N[j] >= len(m.Points)),
			)
			if err != nil {
				_ = et.Add(ErrorBeam{
					BeamIndex: i,
					Detail:    fmt.Sprintf("Point %d", j),
					Err:       fmt.Errorf("outside index of point : %d", b.N[j]),
				})
			}
		}

		if b.N[0] == b.N[1] {
			_ = et.Add(ErrorBeam{
				BeamIndex: i,
				Detail:    "Point indexes",
				Err:       fmt.Errorf("is same"),
			})
		}

		if m.distance(b.N[0], b.N[1]) <= 0 {
			_ = et.Add(ErrorBeam{
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
				_ = et.Add(ErrorBeam{
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

// ErrorPin error in pin data
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
	et := errors.Tree{Name: "checkPins"}
	if len(m.Beams) != len(m.Pins) {
		_ = et.Add(fmt.Errorf("Amount of pins is not same beams, %d != %d",
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
			_ = et.Add(ErrorPin{
				Beam: beam,
				Err: fmt.Errorf("Amount of pins of point is %d more %d",
					pins, maxPins),
			})
		}
		// both X is free
		if m.Pins[beam][0] && m.Pins[beam][3] {
			_ = et.Add(ErrorPin{
				Beam: beam,
				Err:  fmt.Errorf("Not acceptable X direction is free"),
			})
		}
	}

	// check truss elements mistake
	for beam := 0; beam < len(m.Pins); beam++ {
		if m.Pins[beam][2] && m.Pins[beam][5] &&
			!(m.Pins[beam][1] && m.Pins[beam][4]) {
			et.Add(ErrorPin{
				Beam: beam,
				Err:  fmt.Errorf("Pins [1] and [4] in direction Y must be true"),
			})
		}
	}

	if et.IsError() {
		return et
	}
	return nil
}

// ErrorModal error in modal case data
type ErrorModal struct {
	ModalPos int
	Err      error
}

func (e ErrorModal) Error() string {
	return fmt.Sprintf("Error in modal case, mass position %d : %v",
		e.ModalPos,
		e.Err)
}

func (mc *ModalCase) checkInputData(m *Model) (err error) {
	et := errors.Tree{Name: "check input data of modal case"}
	for j := range mc.ModalMasses {
		ld := mc.ModalMasses[j]
		err := isOk(
			isTrue(ld.N < 0),
			isTrue(ld.N >= len(m.Points)),
		)
		if err != nil {
			_ = et.Add(ErrorModal{
				ModalPos: j,
				Err:      fmt.Errorf("outside index of point : %d", ld.N),
			})
		}
		err = isOk(
			isNaN(ld.Mass),
			isInf(ld.Mass),
			isPositive(ld.Mass),
			isNotZero(ld.Mass),
		)
		if err != nil {
			_ = et.Add(ErrorModal{
				ModalPos: j,
				Err:      err,
			})
		}
		for s := 0; s < len(m.Supports); s++ {
			if m.Supports[s][0] || m.Supports[s][1] || m.Supports[s][2] {
				if ld.N == s {
					_ = et.Add(ErrorModal{
						ModalPos: j,
						Err:      fmt.Errorf("Modal mass on support"),
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
