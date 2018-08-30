package hd

// Model is structural calculation model
type Model struct {
	// Points is slice of point coordinate
	//
	// [0] - X coordinate
	// [1] - Y coordinate
	Points [][2]float64 `json:"points"`

	// Beams is slice of point index
	//
	// [0] - begin of beam
	// [1] - end of beam
	Beams [][2]int `json:"beams"`

	// Supports is slice of fixed supports.
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Supports [][3]bool `json:"supports"`
}
