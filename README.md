![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)
[![Go Report Card](https://goreportcard.com/badge/github.com/Konstantin8105/hd)](https://goreportcard.com/report/github.com/Konstantin8105/hd)
[![GoDoc](https://godoc.org/github.com/Konstantin8105/hd?status.svg)](https://godoc.org/github.com/Konstantin8105/hd)

# hd

FEM(finite element method) for structural engineer

**TODO :**
```
#### For calculate the steel standalone flue gas stack


Input data:
* Lower elevation
* Outside diameter
* Thickness of shell
* Temperature (operating, ambient)
* Platforms (elevation, design)
* Ladders (elevations, design)
* Parameters of flanges(F.S.)
* Parameters of doors
* Parameters of brances(external ducts)
* Property of material - linear
* Property of concrete(thk, density)
* Corrosion
* Strakes


Loads and Check codes:
* Russian code: SNiP 2.01.07, SNiP II-23-81
* CICIND
* API560


Finite elements:
* Beam 2D


Type of loads:
* Concentrate loads
* Uniform load
* Uniformly varying load


Support:
* Fixed support


Analyze:
* Static analyze
* Natural frequency analyze (frequency, modes)
* Buckling analyze
* P-delta analyse(Non-linear)


GUI:
* Command line


Output data:
* Optimize the design(thk, outside diameter)
* Deformation
* Result of calculation


#### Advance


Notes:
* Detail design of field splices, connections,...
* Seismic loads


Loads and Check codes:
* UBC-97
* ASME STS-1
* Eurocode
* Checking fatigue


GUI:
* Visualize in 3D
* Lifting devices
```
