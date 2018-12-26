package mod

import (
	"fmt"
	"math"

	"github.com/Konstantin8105/errors"

	"github.com/Konstantin8105/hd"
)

// LoadUniform - convert uniform load on single beam to node load and
// return slice of node load or error if input data not valid.
// Use for selfweight - not-projection load.
//
//	uf is uniform load in global system direction
//	[0] - X , Unit: N
//	[1] - Y , Unit: N
//
//	Example of projection uniform load on beam with different location:
//	uf = [2]float64{ 0.0 , -1.0 }
//	+--+--+
//	|  |  |
//	V  V  V
//	•                •          +-----+-----+
//	 \               |          |     |     |
//	  \              |          V     V     V
//	   \             |          •-----------•
//	    \            |
//	     \           |
//	      •          •
//
//	Example of not-projection uniform load on beam with different location:
//	uf = [2]float64{ 0.0 , -1.0 }
//	|\
//	V \
//	•  \             •          +-----+-----+
//	 \ |\            | |        |     |     |
//	  \V \           | V        V     V     V
//	   \  |          |          •-----------•
//	    \ |          | |
//	     \V          | V
//	      •          •
//
// TODO (KI) : add specific for pin ends, because that function specific only for rigit ends.
func LoadUniform(m *hd.Model, beamIndex int, projection bool, uf [2]float64) (ln []hd.LoadNode, err error) {
	defer func() {
		if err != nil {
			err = errors.New("LoadUniform").Add(err)
		}
	}()

	// check input data
	et := errors.New("check input data")
	if m == nil {
		_ = et.Add(fmt.Errorf("Model is nil"))
	} else {
		if beamIndex >= len(m.Beams) || beamIndex < 0 {
			_ = et.Add(fmt.Errorf("index of beam is outside of model slice[0...%d]: %d", len(m.Beams), beamIndex))
		} else {
			// check point index is exist
			for i, node := range m.Beams[beamIndex].N {
				if node >= len(m.Points) || node < 0 {
					_ = et.Add(fmt.Errorf("index of node %d of beam is outside slice: [0...%d]", i, len(m.Points)))
				}
			}
		}
		for i := range uf {
			if math.IsNaN(uf[i]) {
				_ = et.Add(fmt.Errorf("load not valid load %d : NaN", i))
			}
			if math.IsInf(uf[i], 0) {
				_ = et.Add(fmt.Errorf("load not valid load %d : NaN", i))
			}
		}
	}

	if et.IsError() {
		return nil, et
	}

	// converting

	// calculate dx, dy is projection length of beam
	nodes := m.Beams[beamIndex].N
	dx := m.Points[nodes[0]][0] - m.Points[nodes[1]][0]
	dy := m.Points[nodes[0]][1] - m.Points[nodes[1]][1]

	// load on node:
	// M = (q * l * l) / 12.0

	// by X direction
	switch {
	case dx > 0:
		// location of beam node:
		// start - right
		// end   - left
		ln = append(ln, hd.LoadNode{
			N: nodes[0], // start node
			Forces: [3]float64{
				0.0,                       // X
				0.0,                       // Y
				-(uf[1] * dx * dx) / 12.0, // M
			},
		}, hd.LoadNode{
			N: nodes[1], // end node
			Forces: [3]float64{
				0.0,                      // X
				0.0,                      // Y
				(uf[1] * dx * dx) / 12.0, // M
			},
		})

	case dx < 0:
		// location of beam node:
		// start - left
		// end   - right
		ln = append(ln, hd.LoadNode{
			N: nodes[0], // start node
			Forces: [3]float64{
				0.0,                      // X
				0.0,                      // Y
				(uf[1] * dx * dx) / 12.0, // M
			},
		}, hd.LoadNode{
			N: nodes[1], // end node
			Forces: [3]float64{
				0.0,                       // X
				0.0,                       // Y
				-(uf[1] * dx * dx) / 12.0, // M
			},
		})
	}

	// by Y direction
	switch {
	case dy > 0:
		// location of beam node:
		// start - up
		// end   - down
		ln = append(ln, hd.LoadNode{
			N: nodes[0], // start node
			Forces: [3]float64{
				0.0,                      // X
				0.0,                      // Y
				(uf[0] * dx * dx) / 12.0, // M
			},
		}, hd.LoadNode{
			N: nodes[1], // end node
			Forces: [3]float64{
				0.0,                       // X
				0.0,                       // Y
				-(uf[0] * dx * dx) / 12.0, // M
			},
		})

	case dy < 0:
		// location of beam node:
		// start - down
		// end   - up
		ln = append(ln, hd.LoadNode{
			N: nodes[0], // start node
			Forces: [3]float64{
				0.0,                       // X
				0.0,                       // Y
				-(uf[0] * dx * dx) / 12.0, // M
			},
		}, hd.LoadNode{
			N: nodes[1], // end node
			Forces: [3]float64{
				0.0,                      // X
				0.0,                      // Y
				(uf[0] * dx * dx) / 12.0, // M
			},
		})
	}

	// load on node:
	// P = (q * l) / 2.0
	if projection {
		dx = math.Abs(dx)
		dy = math.Abs(dy)
	} else {
		dx = math.Sqrt(dx*dx + dy*dy)
		dy = dx
	}

	for i := range nodes {
		ln = append(ln,
			// load by X direction
			hd.LoadNode{
				N: nodes[i],
				Forces: [3]float64{
					(uf[0] * dy) / 2.0, // X
					0.0,                // Y
					0.0,                // M
				},
			},
			// load by Y direction
			hd.LoadNode{
				N: nodes[i],
				Forces: [3]float64{
					0.0,                // X
					(uf[1] * dx) / 2.0, // Y
					0.0,                // M
				},
			})
	}

	return
}
