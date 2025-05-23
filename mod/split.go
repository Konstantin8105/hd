package mod

import (
	"fmt"

	"github.com/Konstantin8105/hd"
)

// SplitBeam is split beam on small parts.
// Rules of splitting:
//
//   - beamIndex will beam connected to start beam point
//   - all new point add at the end of point list
//   - all new beam add at the end of beam list
func SplitBeam(m *hd.Model, beamIndex, amountIntermediatePoints int) (err error) {
	if amountIntermediatePoints < 1 {
		return fmt.Errorf("Not valid value of amount intermediate points : %d",
			amountIntermediatePoints)
	}
	if beamIndex < 0 || beamIndex >= len(m.Beams) {
		return fmt.Errorf("Not valid value of beam index : %d",
			beamIndex)
	}

	// create a new points
	lastPointIndex := len(m.Points)
	startPoint := m.Points[m.Beams[beamIndex].N[0]]
	endPoint := m.Points[m.Beams[beamIndex].N[1]]
	m.Points = append(m.Points, make([][2]float64, amountIntermediatePoints)...)
	for j := 0; j < 2; j++ {
		// j = 0 -  X coordinate
		// j = 1 -  Y coordinate
		Δ := (endPoint[j] - startPoint[j]) / float64(amountIntermediatePoints+1)
		for i := 0; i < amountIntermediatePoints; i++ {
			if startPoint[j] == endPoint[j] {
				m.Points[lastPointIndex+i][j] = startPoint[j]
				continue
			}
			delta := float64(i+1) * Δ
			m.Points[lastPointIndex+i][j] = startPoint[j] + delta
		}
	}

	// create a new beams
	lastBeamIndex := len(m.Beams)
	m.Beams = append(m.Beams, make([]hd.BeamProp, amountIntermediatePoints)...)
	for i := 0; i < amountIntermediatePoints; i++ {
		m.Beams[i+lastBeamIndex] = m.Beams[beamIndex]
		m.Beams[i+lastBeamIndex].N[0] = lastPointIndex + i
		m.Beams[i+lastBeamIndex].N[1] = lastPointIndex + i + 1
	}
	m.Beams[len(m.Beams)-1].N[1] = m.Beams[beamIndex].N[1]

	// change original beam
	m.Beams[beamIndex].N[1] = lastPointIndex

	// add free supports
	m.Supports = append(m.Supports, make([][3]bool, amountIntermediatePoints)...)

	// LoadNodes is not changed

	// Change pins for splitted beams
	if len(m.Pins) != 0 {
		m.Pins = append(m.Pins, make([][6]bool, amountIntermediatePoints)...)
		for i := 3; i < 6; i++ {
			m.Pins[len(m.Beams)-1][i] = m.Pins[beamIndex][i]
			m.Pins[beamIndex][i] = false
		}
	}

	return nil
}
