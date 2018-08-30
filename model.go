package hd

// Model is structural calculation model
type Model struct {
	// Points is slice of point coordinate
	//
	// [0] - X coordinate
	// [1] - Y coordinate
	Points [][2]float64

	// Beams is slice of point index and beam property
	//
	Beams []BeamProp

	// Supports is slice of fixed supports.
	// Len of support must be same amount of Points
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Supports [][3]bool

	// LoadCases is slice of load cases
	LoadCases []LoadCase
}

// BeamProp is beam property
type BeamProp struct {
	// Start and end point index
	N1, N2 int

	// A cross-section area
	// Unit : sq. meter.
	A float64

	// J is moment inertia
	// Unit : meter^4
	J float64

	// E is modulus of elasticity
	// Unit : Pa
	E float64

	// G is shear modulus
	// Unit : Pa
	G float64

	// Density
	// Unit : N/m^3
	Density float64

	// Coefficient of thermal expansion
	// Unit: 1/deg.C
	Cte float64
}

// LoadCase is summary combination of loads
type LoadCase struct {
	// LoadNodes is nodal loads
	LoadNodes []LoadNode
}

// LoadNode is node load on specific point
type LoadNode struct {
	// N is point index
	N int

	// Forces is node load on each direction
	//
	// [0] - X
	// [1] - Y
	// [2] - M
	Forces [3]float64
}
