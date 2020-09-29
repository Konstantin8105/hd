package sparse

import (
	"bytes"
	"fmt"
	"io"
	"math"
	"os"
	"runtime/debug"

	"github.com/Konstantin8105/errors"

	"math/rand"
)

// Matrix - sparse matrix.
// Matrix in compressed-column or triplet fotmat.
//
// Name struct in CSparse : cs or cs_sparse
type Matrix struct { // struct cs_sparse
	nzmax int       // maximum number of entries
	m     int       // number of rows
	n     int       // number of columns
	p     []int     // column pointers (size n+1) or col indices (size nzmax)
	i     []int     // row indices, size nzmax
	x     []float64 // numerical values, size nzmax
	nz    int       // # of entries in triplet matrix, -1 for compressed-col
}

func (A *Matrix) Copy() (*Matrix, error) {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}

	if et.IsError() {
		et.Name = "Function Copy: check input data"
		return nil, et
	}

	// coping
	C := new(Matrix)
	C.nzmax = A.nzmax
	C.m = A.m
	C.n = A.n
	C.p = append([]int{}, A.p...)
	C.i = append([]int{}, A.i...)
	C.x = append([]float64{}, A.x...)
	C.nz = A.nz

	return C, nil
}

// Size return size of matrix
func (m *Matrix) Dims() (rows, columns int) {
	return m.m, m.n
}

type Triplet Matrix

func NewTriplet() (*Triplet, error) {
	m, err := cs_spalloc(0, 0, 1, true, tripletFormat)
	return (*Triplet)(m), err
}

// symbolic Cholesky, LU, or QR analysis
type css struct { // struct cs_symbolic
	pinv     []int // inverse row perm. for QR, fill red. perm for Chol
	q        []int // fill-reducing column permutation for LU and QR
	parent   []int // elimination tree for Cholesky and QR
	cp       []int // column pointers for Cholesky, row counts for QR
	leftmost []int // TODO
	m2       int   // # of rows for QR, after adding fictitious rows

	// TODO (KI): lnz, unz is "double" at the base, but I think it is "int"
	lnz int // # entries in L for LU or Cholesky; in V for QR
	unz int // # entries in U for LU; in R for QR
}

// numeric Cholesky, LU, or QR factorization
type csn struct { // struct cs_numeric
	L    *Matrix   // L for LU and Cholesky, V for QR
	U    *Matrix   // U for LU, R for QR, not used for Cholesky
	pinv []int     // partial pivoting for LU
	B    []float64 // beta [0..n-1] for QR
}

// cs_dmperm or cs_scc output
type csd struct { // struct cs_dmperm_results
	p  []int  // size m, row permutation
	q  []int  // size n, column permutation
	r  []int  // size nb+1, block k is rows R[k] to R[k+1]-1 in A(P,Q)
	s  []int  // size nb+1, block k is cols S[k] to S[k+1]-1 in A(P,Q)
	nb int    // # of blocks in fine dmperm decomposition
	rr [5]int // coarse row decomposition
	cc [5]int // coarse column decomposition
}

// Add - additon two sparse matrix with factors.
//
// Matrix A, B is sparse matrix in CSC format.
//
//	C = α*A + β*B
//
// Name function in CSparse : cs_add.
func Add(A *Matrix, B *Matrix, α float64, β float64) (*Matrix, error) {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if B == nil {
		_ = et.Add(fmt.Errorf("matrix B is nil"))
	}
	if B != nil && B.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix B is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil && B != nil {
		if A.m != B.m {
			_ = et.Add(fmt.Errorf("amount of rows in matrixes A and B is not same: %d != %d", A.m, B.m))
		}
		if A.n != B.n {
			_ = et.Add(fmt.Errorf("amount of columns in matrixes A and B is not same: %d != %d", A.n, B.n))
		}
	}
	if math.IsNaN(α) {
		_ = et.Add(fmt.Errorf("factor α is Nan value"))
	}
	if math.IsNaN(β) {
		_ = et.Add(fmt.Errorf("factor β is Nan value"))
	}
	if math.IsInf(α, 0) {
		_ = et.Add(fmt.Errorf("factor α is infinity value"))
	}
	if math.IsInf(β, 0) {
		_ = et.Add(fmt.Errorf("factor β is infinity value"))
	}

	if et.IsError() {
		et.Name = "Function Add: check input data"
		return nil, et
	}

	// internal : acceptable case B.x == nil && A.x == nil

	// TODO (KI) : if factors alpha, beta is zero, then more simplification

	// initialization of variables
	m, anz, n, bnz := A.m, A.p[A.n], B.n, B.p[B.n]

	// get workspace
	values := (A.x != nil && B.x != nil)
	w := make([]int, m)
	defer cs_free(w)

	var x []float64
	defer cs_free(x)
	if values {
		x = make([]float64, m)
	}

	// allocate result
	C, err := cs_spalloc(m, n, anz+bnz, true, cscFormat)
	if err != nil {
		return nil, err
	}
	Cp, Ci, Cx := C.p, C.i, C.x

	// calculation
	var nz int
	for j := 0; j < n; j++ {
		// column j of C starts here
		Cp[j] = nz
		// alpha*A(:,j)
		nz = cs_scatter(A, j, α, w, x, j+1, C, nz)
		// beta*B(:,j)
		nz = cs_scatter(B, j, β, w, x, j+1, C, nz)
		if values {
			for p := Cp[j]; p < nz; p++ {
				Cx[p] = x[Ci[p]]
			}
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	// remove extra space from C
	cs_sprealloc(C, 0)
	// success; free workspace, return C
	return C, nil
}

// cs_wclear - clear w
func cs_wclear(mark, lemax int, w []int, n int) int {
	if mark < 2 || mark+lemax < 0 {
		for k := 0; k < n; k++ {
			if w[k] != 0 {
				w[k] = 1
			}
		}
		mark = 2
	}
	// at this point, w [0..n-1] < mark holds
	return mark
}

type Order uint8

const (
	// natural ordering.
	AmdNatural Order = iota

	// matrix is square. Used for Cholesky or LU factorization of a matrix with
	// entries preliminary on the diagonal and a preliminary symmetric
	// nonzero pattern.
	//
	// amd(A+A')
	AmdChol

	// usually used for LU factorization of unsymmetric matrices.
	//
	// amd(S'*S)
	AmdLU

	// usually used for LU or QR factorization.
	//
	// amd(A'*A)
	AmdQR
)

func (o Order) String() (out string) {
	switch o {
	case AmdNatural:
		return "natural"
	case AmdChol:
		return "amd(A+A')"
	case AmdLU:
		return "amd(S'*S)"
	case AmdQR:
		return "amd(A'*A)"
	}
	return
}

// cs_amd - p = amd(A+A') if symmetric is true, or amd(A'A) otherwise
// order 0:natural, 1:Chol, 2:LU, 3:QR
func cs_amd(order Order, A *Matrix) []int {
	if !(A != nil && A.nz == -1) || order <= 0 || order > 3 {
		// check
		return nil
	}
	var C *Matrix
	// var A2 *Matrix
	// var AT []cs
	// var Cp []int
	// var Ci []int
	// var last []int
	// var W []int
	// var len []int
	// var nv []int
	// var next []int
	// var P []int
	// var head []int
	// var elen []int
	// var degree []int
	// var w []int
	// var hhead []int
	// var ATp []int
	// var ATi []int
	var d int     // int
	var dk int    // int
	var dext int  // int
	var lemax int // int
	var e int     // int
	// var elenk int
	var eln int // int
	var i int   // int
	var j int
	// var k int
	// var k1 int // int
	// var k2 int // int
	// var k3 int
	var jlast int // int
	var ln int    // int
	// var dense int
	// var nzmax int
	// var mindeg int
	var nvi int // int
	var nvj int // int
	// var nvk int
	var mark int //  int
	var wnvi int // int
	var ok bool  // int
	// var cnz int
	var nel int // int
	var p int
	var p1 int // int
	var p2 int
	var p3 int  // int
	var p4 int  // int
	var pj int  // int
	var pk int  // int
	var pk1 int // int
	var pk2 int // int
	var pn int  // int
	// var q int
	// var n int
	// var m int
	// var t int
	var h int // int
	// --- Construct matrix C -----------------------------------------------
	// compute A'
	AT, err := cs_transpose(A, false)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	var (
		m = A.m
		n = A.n
	)
	// find dense threshold
	dense := 16
	if val := int(10.0 * math.Sqrt(float64(n))); dense < val {
		dense = val
	}
	if n-2 < dense {
		dense = n - 2
	}

	switch {
	case order == 1 && n == m:
		// C = A+A'
		var err error
		C, err = Add(A, AT, 0, 0)
		if err != nil {
			return nil
		}

	case order == 2:
		// drop dense columns from AT
		ATp, ATi := AT.p, AT.i

		for p2, j = 0, 0; j < m; j++ {
			// column j of AT starts here
			p = ATp[j]
			// new column j starts here
			ATp[j] = p2
			if ATp[j+1]-p > dense {
				// skip dense col j
				continue
			}
			for ; p < ATp[j+1]; p++ {
				ATi[p2] = ATi[p]
				p2++
			}
		}

		// finalize AT
		ATp[m] = p2
		// A2 = AT'
		A2, err := cs_transpose(AT, false)
		if err != nil {
			fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
			return nil
		}
		// C=A'*A with no dense rows
		if A2 != nil {
			C, err = Multiply(AT, A2)
			if err != nil {
				fmt.Fprintf(os.Stderr, "%v", err) // TODO(KI) error handling
				return nil
			}
		} else {
			C = nil
		}
		cs_free(A2) // TODO (KI) : remove

	default:
		// C=A'*A
		C, err = Multiply(AT, A)
		if err != nil {
			fmt.Fprintf(os.Stderr, "%v", err) // TODO(KI) error handling
			return nil
		}
	}
	cs_free(AT)
	if C == nil {
		return nil
	}
	// drop diagonal entries
	Fkeep(C, func(i, j int, aij float64) bool {
		// cs_diag - keep off-diagonal entries; drop diagonal entries
		return (i != j)
	})

	Cp := C.p
	cnz := Cp[n]
	// allocate result
	P := make([]int, n+1)
	// get workspace
	W := make([]int, 8*(n+1))
	// add elbow room to C
	t := cnz + cnz/5 + 2*n
	if P == nil || W == nil || !cs_sprealloc(C, t) {
		return cs_idone(P, C, W, false)
	}
	var (
		len    = W[0:]
		nv     = W[1*(n+1):]
		next   = W[2*(n+1):]
		head   = W[3*(n+1):]
		elen   = W[4*(n+1):]
		degree = W[5*(n+1):]
		w      = W[6*(n+1):]
		hhead  = W[7*(n+1):]

		// use P as workspace for last
		last = P
	)

	// --- Initialize quotient graph ----------------------------------------
	for k := 0; k < n; k++ {
		len[k] = Cp[k+1] - Cp[k]
	}

	len[n] = 0
	nzmax := C.nzmax
	Ci := C.i

	for i := 0; i <= n; i++ {
		// degree list i is empty
		head[i] = -1
		last[i] = -1
		next[i] = -1
		// hash list i is empty
		hhead[i] = -1
		// node i is just one node
		nv[i] = 1
		// node i is alive
		w[i] = 1
		// Ek of node i is empty
		elen[i] = 0
		// degree of node i
		degree[i] = len[i]
	}
	// clear w
	mark = cs_wclear(0, 0, w, n)
	// n is a dead element
	elen[n] = -2
	// n is a root of assembly tree
	Cp[n] = -1
	// n is a dead element
	w[n] = 0

	// --- Initialize degree lists ------------------------------------------
	for i := 0; i < n; i++ {
		d := degree[i]
		switch {
		case d == 0:
			// node i is empty
			// element i is dead
			elen[i] = -2
			nel++
			// i is a root of assembly tree
			Cp[i] = -1
			w[i] = 0

		case d > dense:
			// node i is dense
			// absorb i into element n
			nv[i] = 0
			// node i is dead
			elen[i] = -1
			nel++
			Cp[i] = flip(n)
			nv[n]++

		default:
			if head[d] != -1 {
				last[head[d]] = i
			}
			// put node i in degree list d
			next[i] = head[d]
			head[d] = i
		}
	}

	var (
		mindeg int
		k      int
		elenk  int
		nvk    int
	)

	// while (selecting pivots) do
	for nel < n {

		// --- Select node of minimum approximate degree --------------------
		for k = -1; mindeg < n && (func() int {
			k = head[mindeg]
			return k
		}()) == -1; mindeg++ {
		}

		if next[k] != -1 {
			last[next[k]] = -1
		}
		// remove k from degree list
		head[mindeg] = next[k]
		// elenk = |Ek|
		elenk = elen[k]
		// # of nodes k represents
		nvk = nv[k]
		// nv[k] nodes of A eliminated
		nel += nvk
		if elenk > 0 && cnz+mindeg >= nzmax {

			// --- Garbage collection -------------------------------------------
			for j = 0; j < n; j++ {
				p = Cp[j]
				if p >= 0 {
					// j is a live node or element
					// save first entry of object
					Cp[j] = Ci[p]
					// first entry is now flip(j)
					Ci[p] = flip(j)
				}
			}

			// scan all of memory
			var q, p int
			for p = 0; p < cnz; {
				if (func() int {
					j = flip(Ci[p])
					p++
					return j
				}()) >= 0 {
					// found object j
					// restore first entry of object
					Ci[q] = Cp[j]
					// new pointer to object j
					Cp[j] = q
					q++
					for k3 := 0; k3 < len[j]-1; k3++ {
						Ci[q] = Ci[p]
						q++
						p++
					}
				}
			}

			// Ci [cnz...nzmax-1] now free
			cnz = q
		}
		// --- Construct new element ----------------------------------------
		dk = 0
		// flag k as in Lk
		nv[k] = -nvk
		p = Cp[k]
		// do in place if elen[k] == 0
		pk1 = cnz
		if elenk == 0 {
			pk1 = p
		}
		pk2 = pk1
		for k1 := 1; k1 <= elenk+1; k1++ {
			if k1 > elenk {
				// search the nodes in k
				e = k
				// list of nodes starts at Ci[pj]
				pj = p
				// length of list of nodes in k
				ln = len[k] - elenk
			} else {
				// search the nodes in e
				e = Ci[p]
				p++
				pj = Cp[e]
				// length of list of nodes in e
				ln = len[e]
			}
			for k2 := 1; k2 <= ln; k2++ {
				i := Ci[pj]
				pj++
				nvi = nv[i]
				if nvi <= 0 {
					// node i dead, or seen
					continue
				}
				// degree[Lk] += size of node i
				dk += nvi
				// negate nv[i] to denote i in Lk
				nv[i] = -nvi
				// place i in Lk
				Ci[pk2] = i
				pk2++
				if next[i] != -1 {
					last[next[i]] = last[i]
				}
				if last[i] != -1 {
					// remove i from degree list
					next[last[i]] = next[i]
				} else {
					head[degree[i]] = next[i]
				}
			}
			if e != k {
				// absorb e into k
				Cp[e] = flip(k)
				// e is now a dead element
				w[e] = 0
			}
		}
		if elenk != 0 {
			// Ci [cnz...nzmax] is free
			cnz = pk2
		}
		// external degree of k - |Lk\i|
		degree[k] = dk
		// element k is in Ci[pk1..pk2-1]
		Cp[k] = pk1
		len[k] = pk2 - pk1
		// k is now an element
		elen[k] = -2
		// --- Find set differences -----------------------------------------
		// clear w if necessary
		mark = cs_wclear(mark, lemax, w, n)

		// scan 1: find |Le\Lk|
		for pk = pk1; pk < pk2; pk++ {
			i = Ci[pk]
			eln = elen[i]
			if eln <= 0 {
				// skip if elen[i] empty
				continue
			}
			// nv [i] was negated
			nvi = -nv[i]
			wnvi = mark - nvi

			// scan Ei
			for p = Cp[i]; p <= Cp[i]+eln-1; p++ {
				e = Ci[p]
				switch {
				case w[e] >= mark:
					// decrement |Le\Lk|
					w[e] -= nvi
				case w[e] != 0:
					// ensure e is a live element
					// 1st time e seen in scan 1
					w[e] = degree[e] + wnvi
				}
			}

		}

		// --- Degree update ------------------------------------------------
		// scan2: degree update
		for pk = pk1; pk < pk2; pk++ {
			// consider node i in Lk
			i = Ci[pk]
			p1 = Cp[i]
			p2 = p1 + elen[i] - 1
			pn = p1

			// scan Ei
			h = 0
			d = 0
			p = p1
			for p = p1; p <= p2; p++ {
				e = Ci[p]
				if w[e] != 0 {
					// e is an unabsorbed element
					// dext = |Le\Lk|
					dext = w[e] - mark
					if dext > 0 {
						// sum up the set differences
						d += dext
						// keep e in Ei
						Ci[pn] = e
						pn++
						// compute the hash of node i
						h += e
					} else {
						// aggressive absorb. e->k
						Cp[e] = flip(k)
						// e is a dead element
						w[e] = 0
					}
				}
			}

			// elen[i] = |Ei|
			elen[i] = pn - p1 + 1
			p3 = pn
			p4 = p1 + len[i]

			// prune edges in Ai
			for p = p2 + 1; p < p4; p++ {
				j = Ci[p]
				nvj = nv[j]
				if nvj <= 0 {
					// node j dead or in Lk
					continue
				}
				// degree(i) += |j|
				d += nvj
				// place j in node list of i
				Ci[pn] = j
				pn++
				// compute hash for node i
				h += j
			}

			if d == 0 {
				// check for mass elimination
				// absorb i into k
				Cp[i] = flip(k)
				nvi = -nv[i]
				// |Lk| -= |i|
				dk -= nvi
				// |k| += nv[i]
				nvk += nvi
				nel += nvi
				nv[i] = 0
				// node i is dead
				elen[i] = -1
			} else {
				// update degree(i)
				degree[i] = func() int {
					if degree[i] < d {
						return degree[i]
					}
					return d
				}()
				// move first node to end
				Ci[pn] = Ci[p3]
				// move 1st el. to end of Ei
				Ci[p3] = Ci[p1]
				// add k as 1st element in of Ei
				Ci[p1] = k
				// new len of adj. list of node i
				len[i] = pn - p1 + 1
				// finalize hash of i
				h = func() int {
					if h < 0 {
						return -h
					}
					return h
				}() % n
				// place i in hash bucket
				next[i] = hhead[h]
				hhead[h] = i
				// save hash of i in last[i]
				last[i] = h
			}
		}

		// scan2 is done
		// finalize |Lk|
		degree[k] = dk
		lemax = func() int {
			if lemax > dk {
				return lemax
			}
			return dk
		}()
		// clear w
		mark = cs_wclear(mark+lemax, lemax, w, n)

		// --- Supernode detection ------------------------------------------
		for pk = pk1; pk < pk2; pk++ {
			i = Ci[pk]
			if nv[i] >= 0 {
				// skip if i is dead
				continue
			}
			// scan hash bucket of node i
			h = last[i]
			i = hhead[h]
			// hash bucket will be empty
			hhead[h] = -1
			for i != -1 && next[i] != -1 {
				ln = len[i]
				eln = elen[i]
				for p = Cp[i] + 1; p <= Cp[i]+ln-1; p++ {
					w[Ci[p]] = mark
				}
				jlast = i

				// compare i with all j
				for j = next[i]; j != -1; {
					ok = (len[j] == ln && elen[j] == eln)
					for p = Cp[j] + 1; ok && p <= Cp[j]+ln-1; p++ {
						if w[Ci[p]] != mark {
							// compare i and j
							ok = false
						}
					}
					if ok {
						// i and j are identical
						// absorb j into i
						Cp[j] = flip(i)
						nv[i] += nv[j]
						nv[j] = 0
						// node j is dead
						elen[j] = -1
						// delete j from hash bucket
						j = next[j]
						next[jlast] = j
					} else {
						// j and i are different
						jlast = j
						j = next[j]
					}
				}

				i = next[i]
				mark++
			}
		}

		// --- Finalize new element------------------------------------------
		// finalize Lk
		p = pk1
		for pk = pk1; pk < pk2; pk++ {
			i = Ci[pk]
			nvi = -nv[i]
			if nvi <= 0 {
				// skip if i is dead
				continue
			}
			// restore nv[i]
			nv[i] = nvi
			// compute external degree(i)
			d = degree[i] + dk - nvi
			d = func() int {
				if d < n-nel-nvi {
					return d
				}
				return n - nel - nvi
			}()
			if head[d] != -1 {
				last[head[d]] = i
			}
			// put i back in degree list
			next[i] = head[d]
			last[i] = -1
			head[d] = i
			// find new minimum degree
			mindeg = func() int {
				if mindeg < d {
					return mindeg
				}
				return d
			}()
			degree[i] = d
			// place i in Lk
			Ci[func() int {
				defer func() {
					p++
				}()
				return p
			}()] = i
		}

		// # nodes absorbed into k
		nv[k] = nvk
		if (func() int {
			len[k] = p - pk1
			return len[k]
		}()) == 0 {
			// length of adj list of element k
			// k is a root of the tree
			Cp[k] = -1
			// k is now a dead element
			w[k] = 0
		}
		if elenk != 0 {
			// free unused space in Lk
			cnz = p
		}
	}

	// --- Postordering -----------------------------------------------------
	// fix assembly tree
	for i = 0; i < n; i++ {
		Cp[i] = flip(Cp[i])
	}

	for j = 0; j <= n; j++ {
		head[j] = -1
	}

	// place unordered nodes in lists
	for j = n; j >= 0; j-- {
		if nv[j] > 0 {
			// skip if j is an element
			continue
		}
		// place j in list of its parent
		next[j] = head[Cp[j]]
		head[Cp[j]] = j
	}

	// place elements in lists
	for e = n; e >= 0; e-- {
		if nv[e] <= 0 {
			// skip unless e is an element
			continue
		}
		if Cp[e] != -1 {
			// place e in list of its parent
			next[e] = head[Cp[e]]
			head[Cp[e]] = e
		}
	}

	// postorder the assembly tree
	for i, k := 0, 0; i <= n; i++ {
		if Cp[i] == -1 {
			k = cs_tdfs(i, k, head, next, P, w)
		}
	}

	return cs_idone(P, C, W, true)
}

// cs_chol - L = chol (A, [pinv parent cp]), pinv is optional
func cs_chol(A *Matrix, S *css) *csn {
	var d float64
	var lki float64
	var Lx []float64
	var Cx []float64
	var top int
	var i int
	var p int
	var k int
	var Li []int
	var Lp []int
	var s []int
	var Cp []int
	var Ci []int
	if !(A != nil && A.nz == -1) || S == nil || S.cp == nil || S.parent == nil {
		return nil
	}
	n := A.n
	// allocate result
	N := new(csn)
	// get csi workspace
	c := make([]int, 2*n)
	// get double workspace
	var (
		x      = make([]float64, n)
		cp     = S.cp
		pinv   = S.pinv
		parent = S.parent
	)
	C := A
	if pinv != nil {
		C = cs_symperm(A, pinv, true)
	}
	// E is alias for A, or a copy E=A(p,p)
	var E *Matrix
	if pinv != nil {
		E = C
	}
	if N == nil || c == nil || x == nil || C == nil {
		return cs_ndone(N, E, c, x, false)
	}
	s = c[n:]
	Cp = C.p
	Ci = C.i
	Cx = C.x
	L, err := cs_spalloc(n, n, int(cp[n]), true, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}

	// allocate result
	N.L = L
	if L == nil {
		return (cs_ndone(N, E, c, x, false))
	}
	Lp = L.p
	Li = L.i
	Lx = L.x
	for k = 0; k < n; k++ {
		c[k] = cp[k]
		Lp[k] = c[k]
	}

	// compute L(k,:) for L*L' = C
	for k = 0; k < n; k++ {
		// --- Nonzero pattern of L(k,:) ------------------------------------
		// find pattern of L(k,:)
		top = cs_ereach(C, k, parent, s, c)
		// x (0:k) is now zero
		x[k] = 0

		// x = full(triu(C(:,k)))
		for p = Cp[k]; p < Cp[k+1]; p++ {
			if Ci[p] <= k {
				x[Ci[p]] = Cx[p]
			}
		}

		// d = C(k,k)
		d = x[k]
		// clear x for k+1st iteration
		x[k] = 0
		for ; top < n; top++ {
			// --- Triangular solve ---------------------------------------------
			// solve L(0:k-1,0:k-1) * x = C(:,k)
			// s [top..n-1] is pattern of L(k,:)
			i = s[top]
			// L(k,i) = x (i) / L(i,i)
			lki = x[i] / Lx[Lp[i]] // TODO (KI): check divided by zero
			// clear x for k+1st iteration
			x[i] = 0
			for p = Lp[i] + 1; p < c[i]; p++ {
				x[Li[p]] -= Lx[p] * lki
			}
			// d = d - L(k,i)*L(k,i)
			d -= lki * lki
			p = c[i]
			c[i]++
			// store L(k,i) in column i
			Li[p] = k
			Lx[p] = lki
		}
		// --- Compute L(k,k) -----------------------------------------------
		if d <= 0 {
			// not pos def
			return cs_ndone(N, E, c, x, false)
		}
		p = c[k]
		c[k]++
		// store L(k,k) = sqrt (d) in column k
		Li[p] = k
		Lx[p] = math.Sqrt(d)
	}

	// finalize L
	Lp[n] = cp[n]
	// success: free E,s,x; return N
	return cs_ndone(N, E, c, x, true)
}

// cs_cholsol - x=A\b where A is symmetric positive definite; b overwritten with solution
func cs_cholsol(order Order, A *Matrix, b []float64) (result bool) {
	if !(A != nil && A.nz == -1) || b == nil {
		// check inputs
		return false
	}
	n := A.n
	// ordering and symbolic analysis
	S := cs_schol(order, A)
	// numeric Cholesky factorization
	N := cs_chol(A, S)
	// get workspace
	x := make([]float64, n)
	ok := (S != nil && N != nil && x != nil)
	if ok {
		// x = P*b
		cs_ipvec(S.pinv, b, x, n)
		// x = L\x
		cs_lsolve(N.L, x)
		// x = L'\x
		cs_ltsolve(N.L, x)
		// b = P'*x
		cs_pvec(S.pinv, x, b, n)
	}
	cs_free(x)
	cs_free(S)
	cs_free(N)
	return ok
}

// Compress - compress triplet matrix T to compressed sparse column(CSC) format.
//
// Name function in CSparse : cs_compress.
func Compress(T *Triplet) (_ *Matrix, err error) {
	// check input data
	et := errors.New("")
	if T == nil {
		_ = et.Add(fmt.Errorf("matrix T is nil"))
	}
	if T != nil && T.nz == -1 {
		_ = et.Add(fmt.Errorf("matrix T is not triplet format"))
	}

	if et.IsError() {
		et.Name = "Function Add: check input data"
		return nil, et
	}

	defer func() {
		if r := recover(); r != nil {
			err = errors.New("Recovery").Add(fmt.Errorf("%v", r)).
				Add(fmt.Errorf("Stack:\n%s", debug.Stack()))
		}
	}()

	m, n, Ti, Tj, Tx, nz := T.m, T.n, T.i, T.p, T.x, T.nz

	// allocate result
	C, err := cs_spalloc(m, n, nz, Tx != nil, cscFormat)
	if err != nil {
		return nil, err
	}

	// get workspace
	w := make([]int, n)
	defer cs_free(w)

	// initialization
	Cp, Ci, Cx := C.p, C.i, C.x

	// Example:
	//
	//	input matrix:
	//	[ 1 2 ]
	//	[ 3 4 ]
	//	triplets:
	//	[0 0 1]
	//	[0 1 2]
	//	[1 0 3]
	//	[1 1 4]
	//	matrix T:
	//	p:[]int    {0, 1, 0, 1}
	//	i:[]int    {0, 0, 1, 1}
	//	x:[]float64{1, 2, 3, 4}
	//	m:2
	//	n:2
	//	nz:4
	//	nzmax:4

	// column counts
	// amount non-zero values in column
	for k := 0; k < nz; k++ {
		w[Tj[k]]++
	}

	// column pointers
	// example:
	// Cp: [0 0 0]
	// w : [2 2]
	_, err = cs_cumsum(Cp, w)
	if err != nil {
		return nil, err
	}
	// example:
	// Cp: [0 2 4] - index in slice `i` for each column
	// w : [0 2]

	// calculation
	for k := 0; k < nz; k++ {
		// A(i,j) is the pth entry in C
		p := w[Tj[k]]
		w[Tj[k]]++
		Ci[p] = Ti[k]
		if Cx != nil {
			Cx[p] = Tx[k]
		}
	}

	//	Result:
	//	p:[]int    {0, 2, 4}
	//	i:[]int    {0, 1, 0, 1}
	//	x:[]float64{1, 3, 2, 4}
	//	m:2
	//	n:2
	//	nz:-1 // indicate CSC format
	//	nzmax:4

	// success
	return C, nil
}

// init_ata - column counts of LL'=A or LL'=A'A, given parent & post ordering
func init_ata(AT *Matrix, post []int, w []int, head *[]int, next *[]int) {
	var (
		m   = AT.n
		n   = AT.m
		ATp = AT.p
		ATi = AT.i
	)
	*head = w[4*n:]
	*next = w[5*n+1:]

	// invert post
	for k := 0; k < n; k++ {
		w[post[k]] = k
	}

	for i := 0; i < m; i++ {
		k := n
		for p := ATp[i]; p < ATp[i+1]; p++ {
			if !(k < w[ATi[p]]) {
				k = w[ATi[p]]
			}
		}

		// place row i in linked list k
		(*next)[i] = (*head)[k]
		(*head)[k] = i
	}
}

// cs_counts -
func cs_counts(A *Matrix, parent []int, post []int, ata bool) []int {
	var J int
	var q int
	var jleaf int
	var head []int
	var next []int
	if !(A != nil && A.nz == -1) || parent == nil || post == nil {
		// check inputs
		return nil
	}
	m := A.m
	n := A.n
	s := 4 * n
	if ata {
		s += (n + m + 1)
	}
	// allocate result
	colcount := make([]int, n)
	delta := colcount
	// get workspace
	w := make([]int, s)
	// AT = A'
	AT, err := cs_transpose(A, false)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	if AT == nil || colcount == nil {
		return cs_idone(colcount, AT, w, false)
	}
	var (
		ancestor = w
		maxfirst = w[n:]
		prevleaf = w[2*n:]
		first    = w[3*n:]
	)

	// clear workspace w [0..s-1]
	for k := 0; k < s; k++ {
		w[k] = -1
	}

	// find first [j]
	for k := 0; k < n; k++ {
		j := post[k]
		// delta[j]=1 if j is a leaf
		delta[j] = 0
		if first[j] == -1 {
			delta[j] = 1
		}
		for ; j != -1 && first[j] == -1; j = parent[j] {
			first[j] = k
		}
	}

	var (
		ATp = AT.p
		ATi = AT.i
	)
	if ata {
		init_ata(AT, post, w, &head, &next)
	}

	// each node in its own set
	for i := 0; i < n; i++ {
		ancestor[i] = i
	}

	for k := 0; k < n; k++ {
		// j is the kth node in postordered etree
		j := post[k]
		if parent[j] != -1 {
			// j is not a root
			delta[parent[j]]--
		}

		// J=j for LL'=A case
		for J = func() int {
			if ata {
				return head[k]
			}
			return j
		}(); J != -1; J = func() int {
			if ata {
				return next[J]
			}
			return -1
		}() {
			for p := ATp[J]; p < ATp[J+1]; p++ {
				i := ATi[p]
				q = cs_leaf(i, j, first, maxfirst, prevleaf, ancestor, &jleaf)
				if jleaf >= 1 {
					// A(i,j) is in skeleton
					delta[j]++
				}
				if jleaf == 2 {
					// account for overlap in q
					delta[q]--
				}
			}
		}

		if parent[j] != -1 {
			ancestor[j] = parent[j]
		}
	}

	// sum up delta's of each child
	for j := 0; j < n; j++ {
		if parent[j] != -1 {
			colcount[parent[j]] += colcount[j]
		}
	}

	// success: free workspace
	return (cs_idone(colcount, AT, w, true))
}

// cs_cumsum - p [0..n] = cumulative sum of c [0..n-1], and then copy p [0..n-1] into c
// n is len of vector c.
//
// Example:
//
//	input data:
//	p =  [0 0 0 0 0]
//	c =  [8 8 8 6]
//	n =  len(c) = 4
//	output data:
//	p =  [0 8 16 24 30]
//	c =  [0 8 16 24]
//	Return int is : 30
//	Error      is : nil
//
// Overflow checking: Let`s take square dense matrix with rows, columns `g`.
// So, vector `c` have length `g` with value `g` in input data.
// Last value of vector `p` in output is `g`*`g`.
// If we suppose type of `p` is `int`, then limit of matrix size `g`:
//
// * for 32-bit machine is math.Sqrt(float64(math.MaxInt32-1)) = 4e+4.
// * for 64-bit machine is math.Sqrt(float64(math.MaxInt64-1)) = 3e+9.
//
func cs_cumsum(p []int, c []int) (int, error) {
	// check input data
	et := errors.New("")
	if p == nil {
		_ = et.Add(fmt.Errorf("Vector p is nil"))
	}
	if c == nil {
		_ = et.Add(fmt.Errorf("Vector c is nil"))
	}
	if len(p) != len(c)+1 {
		_ = et.Add(fmt.Errorf("length of `p` is not length of `c` + 1: %d != %d + 1", len(p), len(c)))
	}

	if et.IsError() {
		et.Name = "Function cs_cumsum: check input data"
		return -1, et
	}

	// calculation
	nz := 0
	for i := range c {
		p[i] = nz
		nz += c[i]
		// also copy p[0..n-1] back into c[0..n-1]
		c[i] = p[i]

		// this is usually happen for int overflow
		if nz < 0 {
			return -1, fmt.Errorf("Function cs_cumsum: value overflow: %d", nz)
		}
	}
	p[len(c)] = nz // add last summ

	// return sum (c [0..n-1])
	return nz, nil
}

// cs_dfs - depth-first-search of the graph of a matrix, starting at node j
func cs_dfs(j int, G *Matrix, top int, xi []int, pstack []int, pinv []int) int {
	if !(G != nil && G.nz == -1) || xi == nil || pstack == nil {
		// check inputs
		return -1
	}
	// initialization
	Gp, Gi := G.p, G.i

	// initialize the recursion stack
	xi[0] = j
	var head int
	for head >= 0 {
		// get j from the top of the recursion stack
		j = xi[head]
		jnew := j
		if pinv != nil {
			jnew = pinv[j]
		}
		if !marked(Gp, j) {
			// mark node j as visited
			Gp[j] = mark(Gp, j)

			if jnew < 0 {
				pstack[head] = 0
			} else {
				pstack[head] = unflip(Gp[jnew])
			}
		}
		// node j done if no unvisited neighbors
		done := true

		var p2 int
		if jnew < 0 {
			p2 = 0
		} else {
			p2 = unflip(Gp[jnew+1])
		}

		// examine all neighbors of j
		for p := pstack[head]; p < p2; p++ {
			// consider neighbor node i
			i := Gi[p]
			if marked(Gp, i) {
				// skip visited node i
				continue
			}
			// pause depth-first search of node j
			pstack[head] = p
			// start dfs at node i
			head++
			xi[head] = i
			// node j is not done
			done = false
			// break, to start dfs (i)
			break
		}

		if done {
			// depth-first search at node j is done
			// remove j from the recursion stack
			head--
			// and place in the output stack
			top--
			xi[top] = j
		}
	}
	return top
}

// cs_bfs - breadth-first search for coarse decomposition (C0,C1,R1 or R0,R3,C3)
func cs_bfs(A *Matrix, n int,
	wi []int,
	wj []int,
	queue []int,
	imatch []int,
	jmatch []int,
	mark int) bool {

	// check input data
	// if wi == nil || wj == nil || queue == nil || imatch == nil || jmatch == nil {
	// 	return false
	// }

	var head int
	var tail int
	var j2 int
	var C *Matrix

	// place all unmatched nodes in queue
	for j := 0; j < n; j++ {
		if imatch[j] >= 0 {
			// skip j if matched
			continue
		}
		// j in set C0 (R0 if transpose)
		wj[j] = 0
		// place unmatched col j in queue
		queue[tail] = j
		tail++
	}

	if tail == 0 {
		// quick return if no unmatched nodes
		return true
	}
	C = A
	if mark != 1 {
		var err error
		C, err = cs_transpose(A, false)
		if err != nil {
			fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
			return false
		}
	}
	if C == nil {
		// bfs of C=A' to find R3,C3 from R0
		return false
	}
	Ap := C.p
	Ai := C.i
	for head < tail {
		// while queue is not empty
		// get the head of the queue
		j := queue[head]
		head++
		for p := Ap[j]; p < Ap[j+1]; p++ {
			i := Ai[p]
			if wi[i] >= 0 {
				// skip if i is marked
				continue
			}
			// i in set R1 (C3 if transpose)
			wi[i] = mark
			// traverse alternating path to j2
			j2 = jmatch[i]
			if wj[j2] >= 0 {
				// skip j2 if it is marked
				continue
			}
			// j2 in set C1 (R3 if transpose)
			wj[j2] = mark
			// add j2 to queue
			queue[tail] = j2
			tail++
		}
	}
	if mark != 1 {
		// free A' if it was created
		cs_free(C)
	}
	return true
}

// cs_matched - collect matched rows and columns into p and q
func cs_matched(n int,
	wj []int,
	imatch []int,
	p []int,
	q []int,
	cc *[5]int,
	rr *[5]int,
	set int,
	mark int) {

	kc := cc[set]
	kr := rr[set-1]

	for j := 0; j < n; j++ {
		if wj[j] != mark {
			// skip if j is not in C set
			continue
		}
		p[kr] = imatch[j]
		kr++
		q[kc] = j
		kc++
	}

	cc[set+1] = kc
	rr[set] = kr
}

// cs_unmatched - collect unmatched rows into the permutation vector p
func cs_unmatched(m int, wi []int, p []int, rr *[5]int, set int) {
	kr := rr[set]
	for i := 0; i < m; i++ {
		if wi[i] == 0 {
			p[kr] = i
			kr++
		}
	}
	rr[set+1] = kr
}

// cs_dmperm - Given A, compute coarse and then fine dmperm
func cs_dmperm(A *Matrix, seed int) *csd {
	var cnz int
	// var nc int
	// var pinv []int
	// var Cp []int
	// var Ci []int
	// var ps []int
	// var rs []int
	// var nb1 int
	// var nb2 int
	// var C []cs
	// var scc []csd
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	// --- Maximum matching -------------------------------------------------
	m := A.m
	n := A.n
	// allocate result
	D := cs_dalloc(m, n)
	if D == nil {
		return nil
	}
	p := D.p
	q := D.q
	r := D.r
	s := D.s
	cc := &D.cc
	rr := &D.rr
	// max transversal
	jmatch := cs_maxtrans(A, seed)
	// imatch = inverse of jmatch
	imatch := jmatch[m:] // jmatch+m
	if jmatch == nil {
		return cs_ddone(D, nil, jmatch, false)
	}
	// --- Coarse decomposition ---------------------------------------------
	// use r and s as workspace
	wi := r
	wj := s

	// unmark all cols for bfs
	for j := 0; j < n; j++ {
		wj[j] = -1
	}

	// unmark all rows for bfs
	for i := 0; i < m; i++ {
		wi[i] = -1
	}

	// find C1, R1 from C0
	cs_bfs(A, n, wi, wj, q, imatch, jmatch, 1)
	// find R3, C3 from R0
	ok := cs_bfs(A, m, wj, wi, p, jmatch, imatch, 3)
	if !ok {
		return (cs_ddone(D, nil, jmatch, false))
	}
	// unmatched set C0
	cs_unmatched(n, wj, q, cc, 0)
	// set R1 and C1
	cs_matched(n, wj, imatch, p, q, cc, rr, 1, 1)
	// set R2 and C2
	cs_matched(n, wj, imatch, p, q, cc, rr, 2, -1)
	// set R3 and C3
	cs_matched(n, wj, imatch, p, q, cc, rr, 3, 3)
	// unmatched set R0
	cs_unmatched(m, wi, p, rr, 3)
	cs_free(jmatch)
	// --- Fine decomposition -----------------------------------------------
	// pinv=p'
	pinv := cs_pinv(p, m)
	if pinv == nil {
		return (cs_ddone(D, nil, nil, false))
	}
	// C=A(p,q) (it will hold A(R2,C2))
	C := cs_permute(A, pinv, q, false)
	cs_free(pinv)
	if C == nil {
		return (cs_ddone(D, nil, nil, false))
	}
	Cp := C.p
	// delete cols C0, C1, and C3 from C
	nc := cc[3] - cc[2]
	if cc[2] > 0 {
		for j := cc[2]; j <= cc[3]; j++ {
			Cp[j-cc[2]] = Cp[j]
		}
	}
	C.n = nc
	if rr[2]-rr[1] < m {
		// delete rows R0, R1, and R3 from C
		Fkeep(C, func(i, j int, aij float64) bool {
			// cs_rprune - return 1 if row i is in R2
			return (i >= rr[1] && i < rr[2])
		})

		cnz = Cp[nc]
		Ci := C.i
		if rr[1] > 0 {
			for k := 0; k < cnz; k++ {
				Ci[k] -= rr[1]
			}
		}
	}
	C.m = nc
	// find strongly connected components of C
	scc := cs_scc(C)
	if scc == nil {
		return (cs_ddone(D, C, nil, false))
	}
	// --- Combine coarse and fine decompositions ---------------------------
	// C(ps,ps) is the permuted matrix
	ps := scc.p
	// kth block is rs[k]..rs[k+1]-1
	rs := scc.r
	// # of blocks of A(R2,C2)
	nb1 := scc.nb
	for k := 0; k < nc; k++ {
		wj[k] = q[ps[k]+cc[2]]
	}
	for k := 0; k < nc; k++ {
		q[k+cc[2]] = wj[k]
	}
	for k := 0; k < nc; k++ {
		wi[k] = p[ps[k]+rr[1]]
	}
	for k := 0; k < nc; k++ {
		p[k+rr[1]] = wi[k]
	}
	// create the fine block partitions
	nb2 := 0
	s[0] = 0
	r[0] = s[0]
	if cc[2] > 0 {
		// leading coarse block A (R1, [C0 C1])
		nb2++
	}

	// coarse block A (R2,C2)
	for k := 0; k < nb1; k++ {
		// A (R2,C2) splits into nb1 fine blocks
		r[nb2] = rs[k] + rr[1]
		s[nb2] = rs[k] + cc[2]
		nb2++
	}

	if rr[2] < m {
		// trailing coarse block A ([R3 R0], C3)
		r[nb2] = rr[2]
		s[nb2] = cc[3]
		nb2++
	}
	r[nb2] = m
	s[nb2] = n
	D.nb = nb2
	cs_free(scc)
	return (cs_ddone(D, C, nil, true))
}

// cs_droptol -
func cs_droptol(A *Matrix, tol float64) (int, error) {
	// keep all large entries
	return Fkeep(A, func(i, j int, aij float64) bool {
		return (math.Abs(aij) > tol)
	})
}

// Dupl - remove duplicate entries from A
//
// Name function in CSparse: cs_dupl
func Dupl(A *Matrix) error {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}

	if et.IsError() {
		et.Name = "Function Dupl: check input data"
		return et
	}

	// initialization
	m, n, Ap, Ai, Ax := A.m, A.n, A.p, A.i, A.x

	// get workspace
	w := make([]int, m)
	defer cs_free(w)

	// row i not yet seen
	for i := 0; i < m; i++ {
		w[i] = -1
	}

	// Example of input:
	// Ap =  [0     4      6]
	// Ai =  [0 1 0 1  0 1]
	// Ax =  [1 3 1 10 2 4]

	var nz int = 0
	for j := 0; j < n; j++ {
		// column j will start at q
		q := nz

		// Example of w:
		// step j = 0 : [-1 -1]
		// step j = 1 : [ 0  1]
		for p := Ap[j]; p < Ap[j+1]; p++ {
			// A(i,j) is nonzero
			i := Ai[p]
			if w[i] >= q {
				// A(i,j) is a duplicate
				Ax[w[i]] += Ax[p]
				continue
			}

			// record where row i occurs
			w[i] = nz
			// keep A(i,j)
			Ai[nz] = i
			Ax[nz] = Ax[p]
			nz++
		}
		// record start of column j
		Ap[j] = q
	}
	// finalize A
	Ap[n] = nz
	// remove extra space from A
	cs_sprealloc(A, 0)

	// Example of output:
	// Ap =  [0    2    4]
	// Ai =  [0 1  0 1 ]
	// Ax =  [2 13 2 4 ]

	return nil
}

// Entry - add an entry to a triplet matrix; return 1 if ok, 0 otherwise
//
// Name function in CSparse : cs_entry.
func Entry(T *Triplet, i, j int, x float64) error {
	// check input data
	et := errors.New("")
	if T == nil {
		_ = et.Add(fmt.Errorf("matrix T is nil"))
	}
	if T != nil && T.nz < 0 {
		_ = et.Add(fmt.Errorf("matrix T is not triplet format"))
	}
	if i < 0 {
		_ = et.Add(fmt.Errorf("index `i` is less zero"))
	}
	if j < 0 {
		_ = et.Add(fmt.Errorf("index `j` is less zero"))
	}
	if math.IsNaN(x) {
		_ = et.Add(fmt.Errorf("value `x` is Nan value"))
	}
	if math.IsInf(x, 0) {
		_ = et.Add(fmt.Errorf("value `x` is infinity value"))
	}
	// memory reallocation
	if T != nil {
		if T.nz >= T.nzmax && !cs_sprealloc((*Matrix)(T), 2*T.nzmax) { // TODO (KI) : add error handling
			_ = et.Add(fmt.Errorf("cannot allocate new vector"))
		}
	}

	if et.IsError() {
		et.Name = "Function Entry: check input data"
		return et
	}

	// calculate amount rows and columns
	if T.m < i+1 {
		T.m = i + 1
	}
	if T.n < j+1 {
		T.n = j + 1
	}

	// do not add zero entry
	if T.x != nil && x == 0.0 {
		return nil
	}

	// add value at the end of vectors
	if T.x != nil {
		T.x[T.nz] = x
	}
	T.i[T.nz] = i
	T.p[T.nz] = j
	T.nz++

	return nil
}

// cs_ereach - find nonzero pattern of Cholesky L(k,1:k-1) using etree and triu(A(:,k))
func cs_ereach(A *Matrix, k int, parent []int, s []int, w []int) int {
	if !(A != nil && A.nz == -1) || parent == nil || s == nil || w == nil {
		// check inputs
		return -1
	}
	// initialization
	n, Ap, Ai := A.n, A.p, A.i
	top := n

	// mark node k as visited
	w[k] = mark(w, k)

	for p := Ap[k]; p < Ap[k+1]; p++ {
		// A(i,k) is nonzero
		i := Ai[p]
		if i > k {
			// only use upper triangular part of A
			continue
		}

		// traverse up etree
		var len int
		for len = 0; !marked(w, i); i = parent[i] {
			// L(k,i) is nonzero
			s[len] = i
			len++

			// mark i as visited
			w[i] = mark(w, i)
		}

		for len > 0 {
			// push path onto stack
			top--
			len--
			s[top] = s[len]
		}
	}

	// unmark all nodes
	for p := top; p < n; p++ {
		w[s[p]] = mark(w, s[p])
	}

	// unmark node k
	w[k] = mark(w, k)

	// s [top..n-1] contains pattern of L(k,:)
	return top
}

// cs_etree - compute the etree of A (using triu(A), or A'A without forming A'A
func cs_etree(A *Matrix, ata bool) []int {
	var i int
	var k int
	var p int
	var m int
	var n int
	var inext int
	var Ap []int
	var Ai []int
	var w []int
	var parent []int
	var ancestor []int
	var prev []int
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	m = A.m
	n = A.n
	Ap = A.p
	Ai = A.i
	// allocate result
	parent = make([]int, n)
	// get workspace
	w = make([]int, (n + func() int {
		if ata {
			return m
		}
		return 0
	}()))
	if w == nil || parent == nil {
		return (cs_idone(parent, nil, w, false))
	}
	ancestor = w
	prev = w[n:]
	if ata {
		for i = 0; i < m; i++ {
			prev[i] = -1
		}
	}
	for k = 0; k < n; k++ {
		// node k has no parent yet
		parent[k] = -1
		// nor does k have an ancestor
		ancestor[k] = -1
		for p = Ap[k]; p < Ap[k+1]; p++ {
			i = int(func() int {
				if ata {
					return (prev[Ai[p]])
				}
				return (Ai[p])
			}())

			// traverse from i to k
			for ; i != -1 && i < k; i = inext {
				// inext = ancestor of i
				inext = ancestor[i]
				// path compression
				ancestor[i] = k
				if inext == -1 {
					// no anc., parent is k
					parent[i] = k
				}
			}

			if ata {
				prev[Ai[p]] = k
			}
		}
	}
	return (cs_idone(parent, nil, w, true))
}

// Fkeep - drop entries for which fkeep(A(i,j)) is false; return nz if OK, else -1 and error
// Name function of CSparse: cs_fkeep
func Fkeep(A *Matrix, fkeep func(i int, j int, x float64) bool) (_ int, err error) {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if fkeep == nil {
		_ = et.Add(fmt.Errorf("function fkeep is nil"))
	}

	if et.IsError() {
		et.Name = "Function Add: check input data"
		return -1, et
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// calculation
	var nz int = 0
	for j := 0; j < n; j++ {
		// get current location of col j
		p := Ap[j]
		// record new location of col j
		Ap[j] = nz
		for ; p < Ap[j+1]; p++ {
			row := Ai[p] // row
			col := j     // column
			val := 1.0   // value, default
			if Ax != nil {
				val = Ax[p] // value
			}
			if fkeep(row, col, val) {
				if Ax != nil {
					// keep A(i,j)
					Ax[nz] = Ax[p]
				}
				Ai[nz] = Ai[p]
				nz++
			}
		}
	}
	// finalize A
	Ap[n] = nz
	// remove extra space from A
	cs_sprealloc(A, 0)
	return nz, nil
}

// Gaxpy - calculate by next formula.
//
// Matrix A is sparse matrix in CSC format.
//
//	y = A*x+y
//
// Name function in CSparse : cs_gaxpy.
func Gaxpy(A *Matrix, x []float64, y []float64) error {
	// check input data
	etGaxpy := errors.New("")
	if A == nil {
		_ = etGaxpy.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = etGaxpy.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if x == nil {
		_ = etGaxpy.Add(fmt.Errorf("Vector x is nil"))
	}
	if y == nil {
		_ = etGaxpy.Add(fmt.Errorf("Vector y is nil"))
	}
	if A != nil {
		if x != nil && A.n != len(x) {
			_ = etGaxpy.Add(fmt.Errorf("Amount of columns in matrix is not equal length of vector x: %d != %d", A.n, len(x)))
		}
		if y != nil && A.m != len(y) {
			_ = etGaxpy.Add(fmt.Errorf("Amount of rows in matrix is not equal length of vector y: %d != %d", A.m, len(y)))
		}
	}

	if etGaxpy.IsError() {
		etGaxpy.Name = "Function Gaxpy: check input data"
		return etGaxpy
	}

	// initialization
	n, Ap, Ai, Ax := A.n, A.p, A.i, A.x

	// calculation
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			y[Ai[p]] += Ax[p] * x[j]
		}
	}
	return nil
}

// cs_happly - apply the ith Householder vector to x
func cs_happly(V *Matrix, i int, beta float64, x []float64) error {
	if !(V != nil && int(V.nz) == -1) || x == nil {
		// check inputs
		return fmt.Errorf("Not correct input data")
	}

	// initialization
	Vp, Vi, Vx := V.p, V.i, V.x
	tau := 0.0

	// tau = v'*x
	for p := Vp[i]; p < Vp[i+1]; p++ {
		tau += Vx[p] * x[Vi[p]]
	}

	// tau = beta*(v'*x)
	tau *= beta

	// x = x - v*tau
	for p := Vp[i]; p < Vp[i+1]; p++ {
		x[Vi[p]] -= Vx[p] * tau
	}

	return nil
}

// cs_house - create a Householder reflection [v,beta,s]=house(x), overwrite x with v,
// * where (I-beta*v*v')*x = s*e1.  See Algo 5.1.1, Golub & Van Loan, 3rd ed.
func cs_house(x []float64, beta []float64, n int) float64 {
	var s float64
	var sigma float64
	var i int
	if x == nil || beta == nil {
		// check inputs
		return float64((-1))
	}
	for i = 1; i < n; i++ {
		sigma += x[i] * x[i]
	}
	if sigma == 0 {
		// s = |x(0)|
		s = math.Abs(x[0])
		beta[0] = float64(func() int {
			if x[0] <= 0 {
				return 2
			}
			return 0
		}())
		x[0] = 1
	} else {
		// s = norm (x)
		s = math.Sqrt(x[0]*x[0] + sigma)
		x[0] = func() float64 {
			if x[0] <= 0 {
				return (x[0] - s)
			}
			return (-sigma / (x[0] + s))
		}()
		beta[0] = -1 / (s * x[0])
	}
	return (s)
}

// cs_ipvec - x(p) = b, for dense vectors x and b; p=NULL denotes identity
func cs_ipvec(p []int, b []float64, x []float64, n int) bool {
	if x == nil || b == nil {
		// check inputs
		return false
	}
	for k := 0; k < n; k++ {
		if p != nil {
			x[p[k]] = b[k]
			continue
		}
		x[k] = b[k]
	}
	return true
}

// cs_leaf - consider A(i,j), node j in ith row subtree and return lca(jprev,j)
func cs_leaf(i int, j int, first []int, maxfirst []int, prevleaf []int, ancestor []int, jleaf *int) int {
	if first == nil || maxfirst == nil || prevleaf == nil || ancestor == nil || jleaf == nil {
		return -1
	}
	*jleaf = 0
	if i <= j || first[j] <= maxfirst[i] {
		// j not a leaf
		return -1
	}
	// update max first[j] seen so far
	maxfirst[i] = first[j]
	// jprev = previous leaf of ith subtree
	jprev := prevleaf[i]
	prevleaf[i] = j
	// j is first or subsequent leaf
	*jleaf = 2
	if jprev == -1 {
		*jleaf = 1
	}
	if *jleaf == 1 {
		// if 1st leaf, q = root of ith subtree
		return int((i))
	}
	var q int
	for q = jprev; q != ancestor[q]; q = ancestor[q] {
	}
	var sparent int
	for s := jprev; s != q; s = sparent {
		// path compression
		sparent = ancestor[s]
		ancestor[s] = q
	}
	// q = least common ancester (jprev,j)
	return q
}

// Load - load a triplet matrix from a file.
//
// Name function in CSparse : cs_load.
func Load(f io.Reader) (T *Triplet, err error) {
	defer func() {
		if err != nil {
			et := errors.New("Function Load")
			_ = et.Add(err)
			err = et
		}
	}()
	if f == nil {
		// check inputs
		err = fmt.Errorf("reader is nil")
		return
	}
	// allocate result
	if T, err = NewTriplet(); err != nil {
		return nil, err
	}
	for {
		var i, j int
		var x float64

		n, err := fmt.Fscanf(f, "%d %d %f\n", &i, &j, &x)
		if err == io.EOF {
			break
		}
		if n != 3 {
			return nil, fmt.Errorf("scan more then 3 variables")
		}
		if err != nil {
			return nil, fmt.Errorf("cannot scan: %v", err)
		}
		if err := Entry(T, i, j, x); err != nil {
			return nil, err
		}
	}
	return T, nil
}

// cs_lsolve - solve Lx=b where x and b are dense.  x=b on input, solution on output.
func cs_lsolve(L *Matrix, x []float64) bool {
	if !(L != nil && L.nz == -1) || x == nil {
		// check inputs
		return false
	}
	// initialization
	n, Lp, Li, Lx := L.n, L.p, L.i, L.x

	// calculation
	for j := 0; j < n; j++ {
		// (KI) : first entry in CSC matrix column is diagonal
		x[j] /= Lx[Lp[j]]
		for p := Lp[j] + 1; p < Lp[j+1]; p++ {
			x[Li[p]] -= Lx[p] * x[j]
		}
	}
	return true
}

// cs_ltsolve - solve L'x=b where x and b are dense.  x=b on input, solution on output.
//	[ L11 L21T ] [x1]   [ b1 ]
//	[  0  L22T ] [x2] = [ b2 ]
func cs_ltsolve(L *Matrix, x []float64) bool {
	if !(L != nil && int(L.nz) == -1) || x == nil {
		// check inputs
		return false
	}
	// initialization
	n, Lp, Li, Lx := L.n, L.p, L.i, L.x

	// calculation
	for j := n - 1; j >= 0; j-- {
		for p := Lp[j] + 1; p < Lp[j+1]; p++ {
			x[j] -= Lx[p] * x[Li[p]]
		}
		x[j] /= Lx[Lp[j]]
	}
	return true
}

// cs_lu - [L,U,pinv]=lu(A, [q lnz unz]). lnz and unz can be guess
func cs_lu(A *Matrix, S *css, tol float64) *csn {
	var col int
	if !(A != nil && A.nz == -1) || S == nil {
		// check inputs
		return nil
	}
	n := A.n
	q := S.q
	lnz := S.lnz
	unz := S.unz
	// get double workspace
	x := make([]float64, n)
	// get csi workspace
	xi := make([]int, 2*n)
	// allocate result
	N := new(csn)
	if x == nil || xi == nil || N == nil {
		return (cs_ndone(N, nil, xi, x, false))
	}
	L, err := cs_spalloc(n, n, lnz, true, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	N.L = L
	// allocate result L
	U, err := cs_spalloc(n, n, unz, true, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	N.U = U
	// allocate result U
	pinv := make([]int, n)
	N.pinv = pinv
	// allocate result pinv
	if L == nil || U == nil || pinv == nil {
		return (cs_ndone(N, nil, xi, x, false))
	}
	Lp, Up := L.p, U.p

	// clear workspace
	for i := 0; i < n; i++ {
		x[i] = 0
	}

	// no rows pivotal yet
	for i := 0; i < n; i++ {
		pinv[i] = -1
	}

	// no cols of L yet
	for k := 0; k <= n; k++ {
		Lp[k] = 0
	}

	unz = 0
	lnz = unz

	// compute L(:,k) and U(:,k)
	for k := 0; k < n; k++ {
		// --- Triangular solve ---------------------------------------------
		// L(:,k) starts here
		Lp[k] = lnz
		// U(:,k) starts here
		Up[k] = unz
		if (lnz+n > L.nzmax && !cs_sprealloc(L, 2*L.nzmax+n)) ||
			(unz+n > U.nzmax && !cs_sprealloc(U, 2*U.nzmax+n)) {
			return cs_ndone(N, nil, xi, x, false)
		}
		Li, Lx, Ui, Ux := L.i, L.x, U.i, U.x

		col = k
		if q != nil {
			col = q[k]
		}

		// x = L\A(:,col)
		top := cs_spsolve(L, A, col, xi, x, pinv, true)
		// --- Find pivot ---------------------------------------------------
		ipiv := -1
		a := -1.0
		for p := top; p < n; p++ {
			// x(i) is nonzero
			i := xi[p]
			if pinv[i] < 0 {
				t := math.Abs(x[i])
				if t > a {
					// row i is not yet pivotal
					// largest pivot candidate so far
					a = t
					ipiv = i
				}
			} else {
				// x(i) is the entry U(pinv[i],k)
				Ui[unz] = pinv[i]
				Ux[unz] = x[i]
				unz++
			}
		}
		if ipiv == -1 || a <= 0 {
			return (cs_ndone(N, nil, xi, x, false))
		}
		if pinv[col] < 0 && math.Abs(x[col]) >= a*tol {
			// tol=1 for  partial pivoting; tol<1 gives preference to diagonal
			ipiv = col
		}
		// --- Divide by pivot ----------------------------------------------
		// the chosen pivot
		pivot := x[ipiv]
		// last entry in U(:,k) is U(k,k)
		Ui[unz] = k
		Ux[unz] = pivot
		unz++
		// ipiv is the kth pivot row
		pinv[ipiv] = k
		// first entry in L(:,k) is L(k,k) = 1
		Li[lnz] = ipiv
		Lx[lnz] = 1
		lnz++

		// L(k+1:n,k) = x / pivot
		for p := top; p < n; p++ {
			i := xi[p]
			if pinv[i] < 0 {
				// x(i) is an entry in L(:,k)
				// save unpermuted row in L
				Li[lnz] = i
				// scale pivot column
				Lx[lnz] = x[i] / pivot
				lnz++
			}
			// x [0..n-1] = 0 for next k
			x[i] = 0
		}

	}

	// --- Finalize L and U -------------------------------------------------
	Lp[n] = lnz
	Up[n] = unz
	// fix row indices of L for final pinv
	Li := L.i
	for p := 0; p < lnz; p++ {
		Li[p] = pinv[Li[p]]
	}
	// remove extra space from L and U
	cs_sprealloc(L, 0)
	cs_sprealloc(U, 0)
	// success
	return (cs_ndone(N, nil, xi, x, true))
}

// cs_lusol - x=A\b where A is unsymmetric; b overwritten with solution
func cs_lusol(order Order, A *Matrix, b []float64, tol float64) bool {
	if !(A != nil && A.nz == -1) || b == nil {
		// check inputs
		return false
	}
	n := A.n
	// ordering and symbolic analysis
	S := cs_sqr(order, A, false)
	// numeric LU factorization
	N := cs_lu(A, S, tol)
	// get workspace
	x := make([]float64, n)
	ok := (S != nil && N != nil && x != nil)
	if ok {
		// x = b(p)
		cs_ipvec(N.pinv, b, x, n)
		// x = L\x
		cs_lsolve(N.L, x)
		// x = U\x
		cs_usolve(N.U, x)
		// b(q) = x
		cs_ipvec(S.q, x, b, n)
	}
	cs_free(x)
	cs_free(S)
	cs_free(N)
	return ok
}

// cs_free - wrapper for free
//
// Names function in CSparse: cs_free, cs_spfree, cs_nfree, cs_sfree, cs_dfree.
func cs_free(p interface{}) {
	if p == nil {
		// if input is nil, then return
		return
	}

	// TODO (KI): remove "return" for reused memory or debugging
	return

	switch v := p.(type) {
	case []float64:
		if v == nil || (v != nil && cap(v) == 0) {
			return
		}
		// TODO : only for debugging
		for i := range v {
			v[i] = -12121212
		}
		// TODO (KI) : fmt.Fprintf(os.Stdout, "Type : %8d %T\n", cap(v), v)

	case []int:
		if v == nil || (v != nil && cap(v) == 0) {
			return
		}
		// TODO : only for debugging
		for i := range v {
			v[i] = -12121212
		}
		// TODO (KI) : fmt.Fprintf(os.Stdout, "Type : %8d %T\n", cap(v), v)

	case LU:
		cs_free(v.s)
		cs_free(v.n)

	case *Matrix:
		if v != nil {
			cs_free(v.p)
			cs_free(v.i)
			cs_free(v.x)
		}

	case *csn:
		if v != nil {
			cs_free(v.L)
			cs_free(v.U)
			cs_free(v.pinv)
			cs_free(v.B)
		}

	case *css:
		if v != nil {
			cs_free(v.pinv)
			cs_free(v.q)
			cs_free(v.parent)
			cs_free(v.cp)
			cs_free(v.leftmost)
		}

	case *csd:
		if v != nil {
			cs_free(v.p)
			cs_free(v.q)
			cs_free(v.r)
			cs_free(v.s)
			// cs_free(v.rr) // ignore type `int[5]`
			// cs_free(v.cc) // ignore type `int[5]`
		}

		// default:
		// TODO (KI) : fmt.Fprintf(os.Stdout,"add memory reusing for type : %T\n", p)
	}
}

// cs_realloc - wrapper for realloc
func cs_realloc(p interface{}, n int, ok *bool) interface{} {

	realloc := func(length, n int) bool {
		// reallocation critetia
		return 2*n < length || // slice is too big
			length < n // slice is too small
	}

	switch v := p.(type) {
	case []int:
		// reallocate memory
		if realloc(len(v), n) {
			arr := make([]int, n)
			copy(arr, v)
			v, arr = arr, v
			cs_free(arr)
		}
		*ok = true
		return v

	case []float64:
		// reallocate memory
		if realloc(len(v), n) {
			arr := make([]float64, n)
			copy(arr, v)
			v, arr = arr, v
			cs_free(arr)
		}
		*ok = true
		return v
	}

	return nil
}

// cs_augment - find an augmenting path starting at column k and extend the match if found
func cs_augment(k int,
	A *Matrix, jmatch []int,
	cheap []int,
	w []int,
	js []int,
	is []int,
	ps []int) {

	var found bool = false
	var p int
	var i int = -1
	Ap := A.p
	Ai := A.i
	var head int
	var j int
	// start with just node k in jstack
	js[0] = k
	for head >= 0 {
		// --- Start (or continue) depth-first-search at node j -------------
		// get j from top of jstack
		j = js[head]
		if w[j] != k {
			// 1st time j visited for kth path
			// mark j as visited for kth path
			w[j] = k
			for p = cheap[j]; p < Ap[j+1] && !found; p++ {
				// try a cheap assignment (i,j)
				i = Ai[p]
				found = (jmatch[i] == -1)
			}
			// start here next time j is traversed
			cheap[j] = p
			if found {
				// column j matched with row i
				is[head] = i
				// end of augmenting path
				break
			}
			// no cheap match: start dfs for j
			ps[head] = Ap[j]
		}

		// --- Depth-first-search of neighbors of j -------------------------
		for p = ps[head]; p < Ap[j+1]; p++ {
			// consider row i
			i = Ai[p]
			if w[jmatch[i]] == k {
				// skip jmatch [i] if marked
				continue
			}
			// pause dfs of node j
			ps[head] = p + 1
			// i will be matched with j if found
			is[head] = i
			// start dfs at column jmatch [i]
			head++
			js[head] = jmatch[i]
			break
		}

		if p == Ap[j+1] {
			// node j is done; pop from stack
			head--
		}
	}
	if found {
		// augment the match if path found:
		for p = head; p >= 0; p-- {
			jmatch[is[p]] = js[p]
		}
	}
}

// cs_maxtrans - find a maximum transveral
//[jmatch [0..m-1]; imatch [0..n-1]]
func cs_maxtrans(A *Matrix, seed int) []int {
	var i int
	var j int
	var k int
	// var n int
	// var m int
	var p int
	// var n2 int
	var m2 int
	// var Ap []int
	// var jimatch []int
	// var w []int
	// var cheap []int
	// var js []int
	// var is []int
	// var ps []int
	// var Ai []int
	// var Cp []int
	var jmatch []int
	var imatch []int
	var q []int
	// var C []cs
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	m := A.m
	Ap := A.p
	Ai := A.i
	jimatch := make([]int, m+n) // cs_calloc(m+n, uint(0)).([]int)
	// allocate result
	w := jimatch
	if jimatch == nil {
		return nil
	}

	// count nonempty rows and columns
	n2 := 0
	for j, k = 0, 0; j < n; j++ {
		if Ap[j] < Ap[j+1] {
			n2++
		}
		for p = Ap[j]; p < Ap[j+1]; p++ {
			w[Ai[p]] = 1
			// count entries already on diagonal
			if j == Ai[p] {
				k++
			}
		}
	}

	if k == func() int {
		if m < n {
			return m
		}
		return n
	}() {
		// quick return if diagonal zero-free
		jmatch = jimatch
		imatch = jimatch[m:]
		for i = 0; i < k; i++ {
			jmatch[i] = i
		}
		for ; i < m; i++ {
			jmatch[i] = -1
		}
		for j = 0; j < k; j++ {
			imatch[j] = j
		}
		for ; j < n; j++ {
			imatch[j] = -1
		}
		return (cs_idone(jimatch, nil, nil, true))
	}
	for i = 0; i < m; i++ {
		m2 += w[i]
	}
	// transpose if needed
	C := A
	if m2 < n2 {
		var err error
		C, err = cs_transpose(A, false)
		if err != nil {
			fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
			return nil
		}
	}
	if C == nil {
		return (cs_idone(jimatch, func() *Matrix {
			if m2 < n2 {
				return C
			}
			return nil
		}(), nil, false))
	}
	n = C.n
	m = C.m
	Cp := C.p
	jmatch = func() []int {
		if m2 < n2 {
			return jimatch[n:]
		}
		return jimatch
	}()
	imatch = func() []int {
		if m2 < n2 {
			return jimatch
		}
		return jimatch[m:]
	}()
	// get workspace
	w = make([]int, 5*n)
	var (
		cheap = w[n:]
		js    = w[2*n:]
		is    = w[3*n:]
		ps    = w[4*n:]
	)

	// for cheap assignment
	for j = 0; j < n; j++ {
		cheap[j] = Cp[j]
	}

	// all columns unflagged
	for j = 0; j < n; j++ {
		w[j] = -1
	}

	// nothing matched yet
	for i = 0; i < m; i++ {
		jmatch[i] = -1
	}

	// q = random permutation
	q = cs_randperm(n, seed)

	// augment, starting at column q[k]
	for k = 0; k < n; k++ {
		cs_augment(func() int {
			if q != nil {
				return q[k]
			}
			return k
		}(), C, jmatch, cheap, w, js, is, ps)
	}

	cs_free(q)

	// find row match
	for j = 0; j < n; j++ {
		imatch[j] = -1
	}

	for i = 0; i < m; i++ {
		if jmatch[i] >= 0 {
			imatch[jmatch[i]] = i
		}
	}
	return (cs_idone(jimatch, func() *Matrix {
		if m2 < n2 {
			return C
		}
		return nil
	}(), w, true))
}

// Multiply - C = A*B
//
// Name function in CSparse : cs_multiply.
func Multiply(A *Matrix, B *Matrix) (*Matrix, error) {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}
	if B == nil {
		_ = et.Add(fmt.Errorf("matrix B is nil"))
	}
	if B != nil && B.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix B is not CSC(Compressed Sparse Column) format"))
	}
	if A != nil && B != nil {
		if A.n != B.m {
			_ = et.Add(fmt.Errorf("amount of columns in matrix A is not same amount of column matrix B: %d != %d", A.n, B.m))
		}
	}

	if et.IsError() {
		et.Name = "Function Add: check input data"
		return nil, et
	}

	// initialization
	m, anz, n := A.m, A.p[A.n], B.n
	Bp, Bi, Bx := B.p, B.i, B.x
	bnz := Bp[n]

	// get workspace
	w := make([]int, m)
	defer cs_free(w)

	// get workspace
	values := (A.x != nil && Bx != nil)
	var x []float64
	defer cs_free(x)
	if values {
		x = make([]float64, m)
	}

	// allocate result
	C, err := cs_spalloc(m, n, anz+bnz, values, cscFormat)
	if err != nil {
		return nil, err
	}
	Cp := C.p

	// calculation
	var nz int
	for j := 0; j < n; j++ {
		if nz+m > C.nzmax && !cs_sprealloc(C, 2*C.nzmax+m) {
			// out of memory
			return nil, fmt.Errorf("Out of memory") ///TODO(KI) error handling
		}
		// C->i and C->x may be reallocated
		Ci, Cx := C.i, C.x

		// column j of C starts here
		Cp[j] = nz
		for p := Bp[j]; p < Bp[j+1]; p++ {
			nz = cs_scatter(A, Bi[p], func() float64 {
				if Bx != nil {
					return Bx[p]
				}
				return 1
			}(), w, x, j+1, C, nz)
		}
		if values {
			for p := Cp[j]; p < nz; p++ {
				Cx[p] = x[Ci[p]]
			}
		}
	}
	// finalize the last column of C
	Cp[n] = nz

	// remove extra space from C
	cs_sprealloc(C, 0)

	// success
	return C, nil
}

// Norm - 1-norm of a sparse matrix = max (sum (abs (A))), largest column sum
//
// Name function in CSparse : cs_norm.
func Norm(A *Matrix) float64 {
	if !(A != nil && A.nz == -1) || A.x == nil {
		// check inputs
		return -1
	}

	// initialization
	n, Ap, Ax := A.n, A.p, A.x

	var norm float64
	for j := 0; j < n; j++ {
		s := 0.0
		for p := Ap[j]; p < Ap[j+1]; p++ {
			s += math.Abs(Ax[p])
		}
		if norm < s {
			norm = s
		}
	}
	return norm
}

// cs_permute - C = A(p,q) where p and q are permutations of 0..m-1 and 0..n-1.
func cs_permute(A *Matrix, pinv []int, q []int, values bool) *Matrix {
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}

	// initialization
	m, n, Ap, Ai, Ax := A.m, A.n, A.p, A.i, A.x

	// alloc result
	C, err := cs_spalloc(m, n, Ap[n], values && Ax != nil, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}

	// initialization
	Cp, Ci, Cx := C.p, C.i, C.x

	// calculation
	nz := 0
	for k := 0; k < n; k++ {
		// column k of C is column q[k] of A
		Cp[k] = nz
		j := k
		if q != nil {
			j = q[k]
		}
		for t := Ap[j]; t < Ap[j+1]; t++ {
			if Cx != nil {
				// row i of A is row pinv[i] of C
				Cx[nz] = Ax[t]
			}
			Ci[nz] = Ai[t]
			if pinv != nil {
				Ci[nz] = pinv[Ai[t]]
			}
			nz++
		}
	}
	// finalize the last column of C
	Cp[n] = nz
	return cs_done(C, nil, nil, true)
}

// cs_pinv - pinv = p', or p = pinv'
func cs_pinv(p []int, n int) []int {
	// var pinv []int
	if p == nil {
		// p = NULL denotes identity
		return nil
	}
	// allocate result
	pinv := make([]int, n)

	// invert the permutation
	for k := 0; k < n; k++ {
		pinv[p[k]] = k
	}

	// return result
	return pinv
}

// cs_post - post order a forest
func cs_post(parent []int, n int) []int {
	// var j int
	var k int
	// var post []int
	// var w []int
	// var head []int
	// var next []int
	// var stack []int
	if parent == nil {
		// check inputs
		return nil
	}
	// allocate result
	post := make([]int, n)
	// get workspace
	w := make([]int, 3*n)
	defer cs_free(w)
	// if w == nil || post == nil {
	// 	return (cs_idone(post, nil, w, false))
	// }
	var (
		head  = w
		next  = w[n:]
		stack = w[2*n:]
	)
	// empty linked lists
	for j := 0; j < n; j++ {
		head[j] = -1
	}

	// traverse nodes in reverse order
	for j := n - 1; j >= 0; j-- {
		if parent[j] == -1 {
			// j is a root
			continue
		}
		// add j to list of its parent
		next[j] = head[parent[j]]
		head[parent[j]] = j
	}

	for j := 0; j < n; j++ {
		if parent[j] != -1 {
			// skip j if it is not a root
			continue
		}
		k = cs_tdfs(j, k, head, next, post, stack)
	}
	// success; free w, return post
	return post //(cs_idone(post, nil, w, true))
}

// Print - print a sparse matrix.
//
//	if brief is true, then print shortly
//
// Name function in CSparse : cs_print.
func (A *Matrix) Print(out io.Writer, brief bool) error {
	if A == nil {
		return fmt.Errorf("Matrix is nil")
	}
	if A.nz != -1 {
		return fmt.Errorf("A is not Matrix type")
	}

	// initialization
	m, n, Ap, Ai, Ax, nzmax := A.m, A.n, A.p, A.i, A.x, A.nzmax

	// print in buffer
	var buf bytes.Buffer
	defer func() {
		if out != nil {
			fmt.Fprintf(out, "%s", buf.String())
		}
	}()

	fmt.Fprintf(&buf, "Sparse\n")

	fmt.Fprintf(&buf, "%d-by-%d, nzmax: %d nnz: %d, 1-norm: %10e\n", m, n, nzmax, Ap[n], Norm(A))
	for j := 0; j < n; j++ {
		fmt.Fprintf(&buf, "    col %d : locations %d to %d\n", j, Ap[j], Ap[j+1]-1)
		for p := Ap[j]; p < Ap[j+1]; p++ {
			fmt.Fprintf(&buf, "      %d : %10e\n", Ai[p], func() float64 {
				if Ax != nil {
					return Ax[p]
				}
				return 1
			}())
			if brief && p > 20 {
				fmt.Fprintf(&buf, "  ...\n")
				return nil
			}
		}
	}
	return nil
}

// Print triplets of matrix
//
//	if brief is true, then print shortly
func (A *Triplet) Print(out io.Writer, brief bool) error {
	if A == nil {
		return fmt.Errorf("Matrix is nil")
	}

	m, n, Ap, Ai, Ax, nzmax, nz := A.m, A.n, A.p, A.i, A.x, A.nzmax, A.nz

	// print in buffer
	var buf bytes.Buffer
	defer func() {
		if out != nil {
			fmt.Fprintf(out, "%s", buf.String())
		}
	}()

	fmt.Fprintf(&buf, "Sparse\n")
	fmt.Fprintf(&buf, "triplet: %d-by-%d, nzmax: %d nnz: %d\n", m, n, nzmax, nz)
	for p := 0; p < nz; p++ {
		fmt.Fprintf(&buf, "    %d %d : %10e\n", Ai[p], Ap[p], func() float64 {
			if Ax != nil {
				return Ax[p]
			}
			return 1
		}())
		if brief && p > 20 {
			fmt.Fprintf(&buf, "  ...\n")
			return nil
		}
	}

	return nil
}

// cs_pvec - x = b(p), for dense vectors x and b; p=NULL denotes identity
func cs_pvec(p []int, b []float64, x []float64, n int) bool {
	if x == nil || b == nil {
		// check inputs
		return false
	}
	for k := 0; k < n; k++ {
		if p != nil {
			x[k] = b[p[k]]
			continue
		}
		x[k] = b[k]
	}
	return true
}

// cs_qr - sparse QR factorization [V,beta,pinv,R] = qr (A)
func cs_qr(A *Matrix, S *css) *csn {
	var Rx []float64
	var Vx []float64
	var Ax []float64
	var x []float64
	var Beta []float64
	var i int
	var k int
	var p int
	// var m int
	var n int
	var vnz int
	var p1 int
	var top int
	var m2 int
	var len int
	var col int
	var rnz int
	var s []int
	var leftmost []int
	var Ap []int
	var Ai []int
	var parent []int
	var Rp []int
	var Ri []int
	var Vp []int
	var Vi []int
	var w []int
	var pinv []int
	var q []int
	var N *csn
	if !(A != nil && A.nz == -1) || S == nil {
		return nil
	}
	// m = A.m
	n = A.n
	Ap = A.p
	Ai = A.i
	Ax = A.x
	q = S.q
	parent = S.parent
	pinv = S.pinv
	m2 = S.m2
	vnz = S.lnz
	rnz = S.unz
	leftmost = S.leftmost
	// get csi workspace
	w = make([]int, m2+n) // cs_malloc(m2+n, uint(0)).([]int)
	// get double workspace
	x = make([]float64, m2) // cs_malloc(int(m2), uint(8)).([]float64)
	// allocate result
	N = new(csn) // cs_calloc(1, uint(32)).([]csn)
	if w == nil || x == nil || N == nil {
		return (cs_ndone(N, nil, w, x, false))
	}
	// s is size n
	s = w[m2:]

	// clear workspace x
	for k = 0; k < m2; k++ {
		x[k] = 0
	}

	V, err := cs_spalloc(int(m2), n, int(vnz), true, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	// allocate result V
	N.L = V
	R, err := cs_spalloc(int(m2), n, int(rnz), true, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	// allocate result R
	N.U = R
	Beta = make([]float64, n) // cs_malloc(n, uint(8)).([]float64)
	// allocate result Beta
	N.B = Beta
	if R == nil || V == nil || Beta == nil {
		return (cs_ndone(N, nil, w, x, false))
	}
	Rp = R.p
	Ri = R.i
	Rx = R.x
	Vp = V.p
	Vi = V.i
	Vx = V.x

	// clear w, to mark nodes
	for i = 0; i < m2; i++ {
		w[i] = -1
	}

	rnz = 0
	vnz = 0

	// compute V and R
	for k = 0; k < n; k++ {
		// R(:,k) starts here
		Rp[k] = rnz
		p1 = vnz
		// V(:,k) starts here
		Vp[k] = p1
		// add V(k,k) to pattern of V
		w[k] = k
		Vi[func() int {
			defer func() {
				vnz++
			}()
			return vnz
		}()] = k
		top = n
		col = int(func() int {
			if q != nil {
				return (q[k])
			}
			return (k)
		}())

		// find R(:,k) pattern
		for p = Ap[col]; p < Ap[col+1]; p++ {
			// i = min(find(A(i,q)))
			i = leftmost[Ai[p]]

			// traverse up to k
			for len = 0; w[i] != k; i = parent[i] {
				s[func() int {
					defer func() {
						len++
					}()
					return len
				}()] = i
				w[i] = k
			}

			for len > 0 {
				// push path on stack
				s[func() int {
					top--
					return top
				}()] = s[func() int {
					len--
					return len
				}()]
			}
			// i = permuted row of A(:,col)
			i = pinv[Ai[p]]
			// x (i) = A(:,col)
			x[i] = Ax[p]
			if i > k && w[i] < k {
				// pattern of V(:,k) = x (k+1:m)
				// add i to pattern of V(:,k)
				Vi[func() int {
					defer func() {
						vnz++
					}()
					return vnz
				}()] = i
				w[i] = k
			}
		}

		// for each i in pattern of R(:,k)
		for p = top; p < n; p++ {
			// R(i,k) is nonzero
			i = s[p]
			// apply (V(i),Beta(i)) to x
			cs_happly(V, int(i), Beta[i], x)
			// R(i,k) = x(i)
			Ri[rnz] = i
			Rx[func() int {
				defer func() {
					rnz++
				}()
				return rnz
			}()] = x[i]
			x[i] = 0
			if parent[i] == k {
				vnz = cs_scatter(V, int(i), 0, w, nil, int(k), V, int(vnz))
			}
		}

		// gather V(:,k) = x
		for p = p1; p < vnz; p++ {
			Vx[p] = x[Vi[p]]
			x[Vi[p]] = 0
		}

		// R(k,k) = norm (x)
		Ri[rnz] = k
		// [v,beta]=house(x)
		Rx[func() int {
			defer func() {
				rnz++
			}()
			return rnz
		}()] = cs_house(Vx[p1:], Beta[k:], vnz-p1)
	}

	// finalize R
	Rp[n] = rnz
	// finalize V
	Vp[n] = vnz
	// success
	return cs_ndone(N, nil, w, x, true)
}

// cs_qrsol - x=A\b where A can be rectangular; b overwritten with solution
func cs_qrsol(order Order, A *Matrix, b []float64) bool {
	var x []float64
	var S *css
	var N *csn
	var AT *Matrix
	var k int
	var ok bool
	if !(A != nil && A.nz == -1) || b == nil {
		// check inputs
		return false
	}
	n := A.n
	m := A.m
	if m >= n {
		// ordering and symbolic analysis
		S = cs_sqr(order, A, true)
		// numeric QR factorization
		N = cs_qr(A, S)
		// get workspace
		if S != nil {
			x = make([]float64, S.m2)
		} else {
			x = make([]float64, 1)
		}
		ok = (S != nil && N != nil && x != nil)
		if ok {
			// x(0:m-1) = b(p(0:m-1)
			cs_ipvec(S.pinv, b, x, m)

			// apply Householder refl. to x
			for k = 0; k < n; k++ {
				cs_happly(N.L, int(k), N.B[k], x)
			}

			// x = R\x
			cs_usolve(N.U, x)
			// b(q(0:n-1)) = x(0:n-1)
			cs_ipvec(S.q, x, b, n)
		}
	} else {
		// Ax=b is underdetermined
		var err error
		AT, err = cs_transpose(A, true)
		if err != nil {
			fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
			return false
		}
		// ordering and symbolic analysis
		S = cs_sqr(order, AT, true)
		// numeric QR factorization of A'
		N = cs_qr(AT, S)
		// get workspace
		if S != nil {
			x = make([]float64, S.m2)
		} else {
			x = make([]float64, 1)
		}
		ok = (AT != nil && S != nil && N != nil && x != nil)
		if ok {
			// x(q(0:m-1)) = b(0:m-1)
			cs_pvec(S.q, b, x, m)
			// x = R'\x
			cs_utsolve(N.U, x)

			// apply Householder refl. to x
			for k = m - 1; k >= 0; k-- {
				cs_happly(N.L, int(k), N.B[k], x)
			}

			// b(0:n-1) = x(p(0:n-1))
			cs_pvec(S.pinv, x, b, n)
		}
	}
	cs_free(x)
	cs_free(S)
	cs_free(N)
	cs_free(AT)
	return ok
}

// cs_randperm - return a random permutation vector, the identity perm, or p = n-1:-1:0.
// * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
// * p = random permutation.
func cs_randperm(n int, seed int) []int {
	if seed == 0 { // || n < 1 {
		// return p = NULL (identity)
		return nil
	}
	// allocate result
	p := make([]int, n)
	if p == nil {
		// out of memory
		return nil
	}
	for k := 0; k < n; k++ {
		p[k] = n - k - 1
	}
	if seed == -1 {
		// return reverse permutation
		return (p)
	}
	// get new random number seed
	rand.Seed(int64(seed))
	for k := 0; k < n; k++ {
		// j = rand integer in range k to n-1
		j := k + (rand.Int() % (n - k))
		// swap p[k] and p[j]
		p[k], p[j] = p[j], p[k]
	}
	return p
}

// cs_reach -
// * xi [top...n-1] = nodes reachable from graph of G*P' via nodes in B(:,k).
// * xi [n...2n-1] used as workspace
func cs_reach(G *Matrix, B *Matrix, k int, xi []int, pinv []int) int {
	if !(G != nil && int(G.nz) == -1) || !(B != nil && int(B.nz) == -1) || xi == nil {
		// check inputs
		return -1
	}

	// initialization
	n, Bp, Bi, Gp := G.n, B.p, B.i, G.p
	top := n

	for p := Bp[k]; p < Bp[k+1]; p++ {
		if !marked(Gp, Bi[p]) {
			// start a dfs at unmarked node i
			top = cs_dfs(Bi[p], G, top, xi, xi[n:], pinv)
		}
	}

	// restore G
	for p := top; p < n; p++ {
		Gp[xi[p]] = mark(Gp, xi[p])
	}

	return top
}

// cs_scatter - x = x + beta * A(:,j), where x is a dense vector and A(:,j) is sparse
func cs_scatter(A *Matrix, j int, beta float64, w []int, x []float64, mark int, C *Matrix, nz int) int {
	if !(A != nil && A.nz == -1) || w == nil || !(C != nil && C.nz == -1) {
		// check inputs
		return -1
	}
	Ap, Ai, Ax, Ci := A.p, A.i, A.x, C.i

	for p := Ap[j]; p < Ap[j+1]; p++ {
		// A(i,j) is nonzero
		i := Ai[p]
		if w[i] < mark {
			// i is new entry in column j
			w[i] = mark
			// add i to pattern of C(:,j)
			Ci[nz] = i
			nz++
			if x != nil {
				// x(i) = beta*A(i,j)
				x[i] = beta * Ax[p]
			}
			continue
		}
		if x != nil {
			// i exists in C(:,j) already
			x[i] += beta * Ax[p]
		}
	}
	return nz
}

// cs_scc - find the strongly connected components of a square matrix
// matrix A temporarily modified, then restored
func cs_scc(A *Matrix) *csd {
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	Ap := A.p
	// allocate result
	D := cs_dalloc(n, 0)
	// AT = A'
	AT, err := cs_transpose(A, false)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	// get workspace
	xi := make([]int, 2*n+1)
	if D == nil || AT == nil || xi == nil {
		return cs_ddone(D, AT, xi, false)
	}
	Blk := xi
	pstack := xi[n:]
	rcopy := pstack
	p := D.p
	r := D.r
	ATp := AT.p
	top := n

	// first dfs(A) to find finish times (xi)
	for i := 0; i < n; i++ {
		if !marked(Ap, i) {
			top = cs_dfs(i, A, top, xi, pstack, nil)
		}
	}

	// restore A; unmark all nodes
	for i := 0; i < n; i++ {
		Ap[i] = mark(Ap, i)
	}

	top = n
	nb := n

	// dfs(A') to find strongly connnected comp
	for k := 0; k < n; k++ {
		// get i in reverse order of finish times
		i := xi[k]
		if marked(ATp, i) {
			// skip node i if already ordered
			continue
		}
		// node i is the start of a component in p
		r[nb] = top
		nb--
		top = cs_dfs(i, AT, top, p, pstack, nil)
	}

	// first block starts at zero; shift r up
	r[nb] = 0
	for k := nb; k <= n; k++ {
		r[k-nb] = r[k]
	}
	nb = n - nb
	// nb = # of strongly connected components
	D.nb = nb

	// sort each block in natural order
	for b := 0; b < nb; b++ {
		for k := r[b]; k < r[b+1]; k++ {
			Blk[p[k]] = b
		}
	}

	for b := 0; b <= nb; b++ {
		rcopy[b] = r[b]
	}
	for i := 0; i < n; i++ {
		p[func() int {
			tempVar := &rcopy[Blk[i]]
			defer func() {
				*tempVar++
			}()
			return *tempVar
		}()] = i
	}
	return (cs_ddone(D, AT, xi, true))
}

// cs_schol - ordering and symbolic analysis for a Cholesky factorization
func cs_schol(order Order, A *Matrix) (result *css) {
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	// allocate result S
	S := new(css)
	// P = amd(A+A'), or natural
	P := cs_amd(order, A)
	// find inverse permutation
	S.pinv = cs_pinv(P, n)
	cs_free(P)
	if order != 0 && S.pinv == nil {
		cs_free(S)
		return nil
	}
	// C = spones(triu(A(P,P)))
	C := cs_symperm(A, S.pinv, false)
	// find etree of C
	S.parent = cs_etree(C, false)
	// postorder the etree
	post := cs_post(S.parent, n)
	// find column counts of chol(C)
	c := cs_counts(C, S.parent, post, false)
	cs_free(post)
	cs_free(C)
	// allocate result S->cp
	S.cp = make([]int, n+1)
	var err error
	S.lnz, err = cs_cumsum(S.cp, c)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	// find column pointers for L
	S.unz = S.lnz
	cs_free(c)
	return (func() *css {
		if S.lnz >= 0 {
			return S
		}
		cs_free(S)
		return nil
	}())
}

// cs_spsolve - solve Gx=b(:,k), where G is either upper (lo=0) or lower (lo=1) triangular
func cs_spsolve(G *Matrix, B *Matrix, k int, xi []int, x []float64, pinv []int, lo bool) int {
	if !(G != nil && G.nz == -1) || !(B != nil && B.nz == -1) || xi == nil || x == nil {
		return -1
	}

	// initialization
	Gp, Gi, Gx, n, Bp, Bi, Bx := G.p, G.i, G.x, G.n, B.p, B.i, B.x

	// xi[top..n-1]=Reach(B(:,k))
	top := cs_reach(G, B, k, xi, pinv)

	// clear x
	for p := top; p < n; p++ {
		x[xi[p]] = 0
	}

	// scatter B
	for p := Bp[k]; p < Bp[k+1]; p++ {
		x[Bi[p]] = Bx[p]
	}

	for px := top; px < n; px++ {
		// x(j) is nonzero
		j := xi[px]
		// j maps to col J of G
		J := j
		if pinv != nil {
			J = pinv[j]
		}
		if J < 0 {
			// column J is empty
			continue
		}
		// x(j) /= G(j,j)
		x[j] /= Gx[func() int {
			if lo {
				return Gp[J]
			}
			return (Gp[J+1] - 1)
		}()]
		// lo: L(j,j) 1st entry
		p := Gp[J]
		if lo {
			p += 1
		}
		// up: U(j,j) last entry
		q := Gp[J+1] - 1
		if lo {
			q += 1
		}
		for ; p < q; p++ {
			// x(i) -= G(i,j) * x(j)
			x[Gi[p]] -= Gx[p] * x[j]
		}
	}
	// return top of stack
	return top
}

// cs_vcount - compute nnz(V) = S->lnz, S->pinv, S->leftmost, S->m2 from A and S->parent
func cs_vcount(A *Matrix, S *css) bool {

	// if A == nil || S == nil {
	// 	return false
	// }

	var i int
	var k int
	var p int
	var pa int
	var n int = A.n
	var m int = A.m
	var Ap []int = A.p
	var Ai []int = A.i
	var next []int
	var head []int
	var tail []int
	var nque []int
	var pinv []int
	var leftmost []int
	var w []int
	var parent []int = S.parent
	pinv = make([]int, m+n)
	// allocate pinv,
	S.pinv = pinv
	leftmost = make([]int, m)
	// and leftmost
	S.leftmost = leftmost
	// get workspace
	w = make([]int, m+3*n)
	if pinv == nil || w == nil || leftmost == nil {
		// pinv and leftmost freed later
		cs_free(w)
		// out of memory
		return false
	}
	next = w
	head = w[m:]
	tail = w[m+n:]
	nque = w[m+2*n:]

	// queue k is empty
	for k = 0; k < n; k++ {
		head[k] = -1
	}

	for k = 0; k < n; k++ {
		tail[k] = -1
	}
	for k = 0; k < n; k++ {
		nque[k] = 0
	}
	for i = 0; i < m; i++ {
		leftmost[i] = -1
	}
	for k = n - 1; k >= 0; k-- {
		for p = Ap[k]; p < Ap[k+1]; p++ {
			// leftmost[i] = min(find(A(i,:)))
			leftmost[Ai[p]] = k
		}
	}

	// scan rows in reverse order
	for i = m - 1; i >= 0; i-- {
		// row i is not yet ordered
		pinv[i] = -1
		k = leftmost[i]
		if k == -1 {
			// row i is empty
			continue
		}
		if func() int {
			tempVar := &nque[k]
			defer func() {
				*tempVar++
			}()
			return *tempVar
		}() == 0 {
			// first row in queue k
			tail[k] = i
		}
		// put i at head of queue k
		next[i] = head[k]
		head[k] = i
	}

	S.lnz = 0
	S.m2 = m

	// find row permutation and nnz(V)
	for k = 0; k < n; k++ {
		// remove row i from queue k
		i = head[k]
		// count V(k,k) as nonzero
		S.lnz++
		if i < 0 {
			// add a fictitious row
			i = func() int {
				tempVar := &S.m2
				defer func() {
					*tempVar++
				}()
				return *tempVar
			}()
		}
		// associate row i with V(:,k)
		pinv[i] = k
		if func() int {
			tempVar := &nque[k]
			*tempVar--
			return *tempVar
		}() <= 0 {
			// skip if V(k+1:m,k) is empty
			continue
		}
		// nque [k] is nnz (V(k+1:m,k))
		S.lnz += nque[k]
		if (func() int {
			pa = parent[k]
			return pa
		}()) != -1 {
			if nque[pa] == 0 {
				// move all rows to parent of k
				tail[pa] = tail[k]
			}
			next[tail[k]] = head[pa]
			head[pa] = next[i]
			nque[pa] += nque[k]
		}
	}

	for i = 0; i < m; i++ {
		if pinv[i] < 0 {
			pinv[i] = k
			k++
		}
	}
	cs_free(w)
	return true
}

// cs_sqr - symbolic ordering and analysis for QR or LU
func cs_sqr(order Order, A *Matrix, qr bool) *css {
	var ok bool = true
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n := A.n
	// allocate result S
	S := new(css)

	// fill-reducing ordering
	S.q = cs_amd(order, A)
	if order != 0 && S.q == nil {
		cs_free(S)
		return nil
	}
	if qr {
		C := A
		if order != 0 {
			C = cs_permute(A, nil, S.q, false)
		}
		// QR symbolic analysis
		// etree of C'*C, where C=A(:,q)
		S.parent = cs_etree(C, true)
		post := cs_post(S.parent, n)
		// col counts chol(C'*C)
		S.cp = cs_counts(C, S.parent, post, true)
		cs_free(post)
		ok = (C != nil && S.parent != nil && S.cp != nil && cs_vcount(C, S))
		if ok {
			S.unz = 0
			for k := 0; k < n; k++ {
				S.unz += S.cp[k]
			}
		}
		if order != 0 {
			cs_free(C)
		}
	} else {
		// for LU factorization only, guess nnz(L) and nnz(U)
		// (KI) : Amount of memory allocated for matrix L and U
		// (KI) : by default is can be changed.
		// (KI) : Before use formula:  4*A.p[n] + n
		// (KI) : Acceptable any more 1.
		S.unz = 4*A.p[n] + n
		S.lnz = S.unz
	}
	// return result S
	if ok {
		return S
	}
	cs_free(S)
	return nil
}

// cs_symperm - C = A(p,p) where A and C are symmetric the upper part stored; pinv not p
func cs_symperm(A *Matrix, pinv []int, values bool) *Matrix {
	var i int
	var j int
	var p int
	var q int
	var i2 int
	var j2 int
	var n int
	var Ap []int
	var Ai []int
	var Cp []int
	var Ci []int
	var w []int
	var Cx []float64
	var Ax []float64
	if !(A != nil && A.nz == -1) {
		// check inputs
		return nil
	}
	n = A.n
	Ap = A.p
	Ai = A.i
	Ax = A.x
	// alloc result
	C, err := cs_spalloc(n, n, Ap[n], values && Ax != nil, cscFormat)
	if err != nil {
		fmt.Fprintf(os.Stderr, err.Error()) // TODO (KI) error hanling
		return nil
	}
	// get workspace
	w = make([]int, n)
	if C == nil {
		// out of memory
		return (cs_done(C, w, nil, false))
	}
	Cp = C.p
	Ci = C.i
	Cx = C.x

	// count entries in each column of C
	for j = 0; j < n; j++ {
		// column j of A is column j2 of C
		j2 = int(func() int {
			if pinv != nil {
				return (pinv[j])
			}
			return (j)
		}())
		for p = Ap[j]; p < Ap[j+1]; p++ {
			i = Ai[p]
			if i > j {
				// skip lower triangular part of A
				continue
			}
			// row i of A is row i2 of C
			i2 = int(func() int {
				if pinv != nil {
					return (pinv[i])
				}
				return (i)
			}())
			// column count of C
			w[func() int {
				if i2 > j2 {
					return (i2)
				}
				return (j2)
			}()]++
		}
	}

	// compute column pointers of C
	_, err = cs_cumsum(Cp, w)
	if err != nil {
		fmt.Fprintf(os.Stderr, "%v\n", err) // TODO (KI): add error hangling
		return nil
	}
	for j = 0; j < n; j++ {
		// column j of A is column j2 of C
		j2 = j
		if pinv != nil {
			j2 = pinv[j]
		}
		for p = Ap[j]; p < Ap[j+1]; p++ {
			i = Ai[p]
			if i > j {
				// skip lower triangular part of A
				continue
			}
			// row i of A is row i2 of C
			i2 = i
			if pinv != nil {
				i2 = pinv[i]
			}

			Ci[(func() int {
				q = func() int {
					tempVar := &w[func() int {
						if i2 > j2 {
							return (i2)
						}
						return (j2)
					}()]
					defer func() {
						*tempVar++
					}()
					return *tempVar
				}()
				return q
			}())] = int(func() int {
				if i2 < j2 {
					return (i2)
				}
				return (j2)
			}())
			if Cx != nil {
				Cx[q] = Ax[p]
			}
		}
	}
	// success; free workspace, return C
	return (cs_done(C, w, nil, true))
}

// cs_tdfs - depth-first search and postorder of a tree rooted at node j
func cs_tdfs(j, k int, head []int, next []int, post []int, stack []int) int {
	if head == nil || next == nil || post == nil || stack == nil {
		// check inputs
		return -1
	}
	// place j on the stack
	stack[0] = j
	top := 0
	for top >= 0 {
		// while (stack is not empty)
		// p = top of stack
		p := stack[top]
		// i = youngest child of p
		i := head[p]
		if i == -1 {
			// p has no unordered children left
			top--
			// node p is the kth postordered node
			post[k] = p
			k++
			continue
		}

		// remove i from children of p
		head[p] = next[i]
		// start dfs on child node i
		top++
		stack[top] = i
	}
	return k
}

// Transpose - C = A'
//
// Name function in CSparse : cs_transpose.
func Transpose(A *Matrix) (*Matrix, error) {
	return cs_transpose(A, true)
}

// cs_transpose - C = A'.
// if values == true, then initialize vector x in Matrix
func cs_transpose(A *Matrix, values bool) (*Matrix, error) {
	// check input data
	et := errors.New("")
	if A == nil {
		_ = et.Add(fmt.Errorf("matrix A is nil"))
	}
	if A != nil && A.nz != -1 {
		_ = et.Add(fmt.Errorf("matrix A is not CSC(Compressed Sparse Column) format"))
	}

	if et.IsError() {
		et.Name = "Function Transpose: check input data"
		return nil, et
	}

	// initialization
	m, n, Ap, Ai, Ax := A.m, A.n, A.p, A.i, A.x

	// allocate result
	C, err := cs_spalloc(n, m, Ap[n], values && Ax != nil, cscFormat)
	if err != nil {
		return nil, err
	}
	// get workspace
	w := make([]int, m)
	defer cs_free(w)

	// initialization
	Cp, Ci, Cx := C.p, C.i, C.x

	// row counts
	for p := 0; p < Ap[n]; p++ {
		w[Ai[p]]++
	}

	// row pointers
	_, err = cs_cumsum(Cp, w)
	if err != nil {
		return nil, err
	}
	for j := 0; j < n; j++ {
		for p := Ap[j]; p < Ap[j+1]; p++ {
			// place A(i,j) as entry C(j,i)
			q := w[Ai[p]]
			Ci[q] = j
			w[Ai[p]]++
			if Cx != nil {
				Cx[q] = Ax[p]
			}
		}
	}
	// success
	return C, nil
}

// cs_updown - sparse Cholesky update/downdate, L*L' + sigma*w*w' (sigma = +1 or -1)
//
// Name function in CSparse : cs_updown.
func cs_updown(L *Matrix, sigma int, C *Matrix, parent []int) int {
	var n int
	var p int
	var f int
	var j int
	var Lp []int
	var Li []int
	var Cp []int
	var Ci []int
	var Lx []float64
	var Cx []float64
	var alpha float64
	var beta float64 = 1
	var delta float64
	var gamma float64
	var w1 float64
	var w2 float64
	var w []float64
	var beta2 float64 = 1
	if !(L != nil && L.nz == -1) || !(C != nil && int(C.nz) == -1) || parent == nil {
		// check inputs
		return 0
	}
	Lp = L.p
	Li = L.i
	Lx = L.x
	n = L.n
	Cp = C.p
	Ci = C.i
	Cx = C.x
	p = Cp[0]
	if p >= Cp[1] {
		// return if C empty
		return 1
	}
	// get workspace
	w = make([]float64, n)
	f = Ci[p]
	for ; p < Cp[1]; p++ {
		// f = min (find (C))
		f = int(func() int {
			if f < Ci[p] {
				return f
			}
			return Ci[p]
		}())
	}

	// clear workspace w
	for j = f; j != -1; j = parent[j] {
		w[j] = 0
	}

	// w = C
	for p = Cp[0]; p < Cp[1]; p++ {
		w[Ci[p]] = Cx[p]
	}

	// walk path f up to root
	for j = f; j != -1; j = parent[j] {
		p = Lp[j]
		// alpha = w(j) / L(j,j)
		alpha = w[j] / Lx[p]
		beta2 = beta*beta + float64(sigma)*alpha*alpha
		if beta2 <= 0 {
			// not positive definite
			break
		}
		beta2 = math.Sqrt(beta2)
		delta = func() float64 {
			if sigma > 0 {
				return (beta / beta2)
			}
			return (beta2 / beta)
		}()
		gamma = float64(sigma) * alpha / (beta2 * beta)
		Lx[p] = delta*Lx[p] + func() float64 {
			if sigma > 0 {
				return (gamma * w[j])
			}
			return 0
		}()
		beta = beta2
		for p += 1; p < Lp[j+1]; p++ {
			w1 = w[Li[p]]
			w2 = w1 - alpha*Lx[p]
			w[Li[p]] = w2
			Lx[p] = delta*Lx[p] + gamma*func() float64 {
				if sigma > 0 {
					return w1
				}
				return w2
			}()
		}
	}

	cs_free(w)
	if beta2 > 0 {
		return 1
	}
	return 0
}

// cs_usolve - solve Ux=b where x and b are dense.  x=b on input, solution on output.
//	[ U11 U12 ] [x1]   [ b1 ]
//	[  0  U22 ] [x2] = [ b2 ]
func cs_usolve(U *Matrix, x []float64) bool {
	if !(U != nil && U.nz == -1) || x == nil {
		// check inputs
		return false
	}

	// initialization
	n, Up, Ui, Ux := U.n, U.p, U.i, U.x

	// calculation
	for j := n - 1; j >= 0; j-- {
		x[j] /= Ux[Up[j+1]-1]
		for p := Up[j]; p < Up[j+1]-1; p++ {
			x[Ui[p]] -= Ux[p] * x[j]
		}
	}
	return true
}

type matrixFormat bool

const (
	tripletFormat matrixFormat = false
	cscFormat     matrixFormat = true
)

// cs_spalloc - allocate a sparse matrix (triplet form or compressed-column form)
func cs_spalloc(m, n, nzmax int, values bool, mf matrixFormat) (*Matrix, error) {
	// check input data
	et := errors.New("")
	if m < 0 {
		_ = et.Add(fmt.Errorf("Value m is less zero : %d", m))
	}
	if n < 0 {
		_ = et.Add(fmt.Errorf("Value n is less zero : %d", n))
	}
	if nzmax < 0 {
		_ = et.Add(fmt.Errorf("Value nzmax is less zero : %d", nzmax))
	}

	if et.IsError() {
		et.Name = "Function cs_spalloc: check input data"
		return nil, et
	}

	// create a new struct
	A := new(Matrix)

	// define dimensions and nzmax
	A.m, A.n = m, n
	if nzmax < 1 {
		nzmax = 1
	}
	A.nzmax = nzmax

	// allocate triplet or comp.col
	switch mf {
	case tripletFormat:
		A.p = make([]int, nzmax)
		A.nz = 0
	case cscFormat:
		A.p = make([]int, n+1)
		A.nz = -1
	}
	A.i = make([]int, nzmax)
	A.x = nil
	if values {
		A.x = make([]float64, nzmax)
	}
	return A, nil
}

// cs_sprealloc - change the max # of entries sparse matrix
func cs_sprealloc(A *Matrix, nzmax int) (result bool) {
	var oki bool
	var okj bool = true
	var okx bool = true
	if A == nil {
		return false
	}
	if nzmax <= 0 {
		nzmax = func() int {
			if A != nil && A.nz == -1 {
				return A.p[A.n]
			}
			return A.nz
		}()
	}
	if nzmax < 1 {
		nzmax = 1
	}

	A.i = cs_realloc(A.i, nzmax, &oki).([]int)
	if A != nil && A.nz >= 0 {
		A.p = cs_realloc(A.p, nzmax, &okj).([]int)
	}
	if A.x != nil {
		A.x = cs_realloc(A.x, nzmax, &okx).([]float64)
	}
	ok := oki && okj && okx
	if ok {
		A.nzmax = nzmax
	}
	return ok
}

// cs_dalloc - allocate a cs_dmperm or cs_scc result
func cs_dalloc(m, n int) *csd {
	// if m < 1 || n < 1 { // TODO (KI) error handling
	// 	return nil
	// }
	D := new(csd)

	D.p = make([]int, m)
	D.r = make([]int, m+6)
	D.q = make([]int, n)
	D.s = make([]int, n+6)

	return D
}

// cs_done - free workspace and return a sparse matrix result
func cs_done(C *Matrix, w []int, x []float64, ok bool) *Matrix {
	cs_free(w)
	cs_free(x)
	// return result if OK, else free it
	if ok {
		return C
	}
	return nil
}

// cs_idone - free workspace and return csi array result
func cs_idone(p []int, C *Matrix, w interface{}, ok bool) []int {
	cs_free(C)
	cs_free(w)
	// return result, or free it
	if ok {
		return p
	}
	cs_free(p)
	return nil
}

// cs_ndone - free workspace and return a numeric factorization (Cholesky, LU, or QR)
func cs_ndone(N *csn, C *Matrix, w interface{}, x interface{}, ok bool) *csn {
	cs_free(C)
	cs_free(w)
	cs_free(x)
	// return result if OK, else free it
	if ok {
		return N
	}
	cs_free(N)
	return nil
}

// cs_ddone - free workspace and return a csd result
func cs_ddone(D *csd, C *Matrix, w interface{}, ok bool) *csd {
	cs_free(C)
	cs_free(w)
	// return result if OK, else free it
	if ok {
		return D
	}
	cs_free(D)
	return nil
}

// cs_utsolve - solve U'x=b where x and b are dense.  x=b on input, solution on output.
func cs_utsolve(U *Matrix, x []float64) bool {
	if !(U != nil && U.nz == -1) || x == nil {
		// check inputs
		return false
	}

	// initialization
	n, Up, Ui, Ux := U.n, U.p, U.i, U.x

	// calculation
	for j := 0; j < n; j++ {
		for p := Up[j]; p < Up[j+1]-1; p++ {
			x[j] -= Ux[p] * x[Ui[p]]
		}
		x[j] /= Ux[Up[j+1]-1]
	}
	return true
}

// Name function in CSparse: CS_FLIP
func flip(i int) int {
	return -i - 2
}

// Name function in CSparse: CS_UNFLIP
func unflip(i int) int {
	if i < 0 {
		return flip(i)
	}
	return i
}

// Name function in CSparse: CS_MARKED
func marked(w []int, j int) bool {
	return w[j] < 0
}

// Name function in CSparse: CS_MARK
func mark(w []int, j int) int {
	w[j] = flip(w[j])
	return w[j]
}
