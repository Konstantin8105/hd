package hd

import "github.com/disiqueira/gotree"

type ErrorTree struct {
	Name string
	errs []error
}

func (e *ErrorTree) Add(err error) {
	e.errs = append(e.errs, err)
}

func (e ErrorTree) Error() (s string) {
	return e.getTree().Print()
}

func (e ErrorTree) IsError() bool {
	return len(e.errs) > 0
}

func (e ErrorTree) getTree() gotree.Tree {
	name := "+"
	if e.Name != "" {
		name = e.Name
	}
	t := gotree.New(name)
	for _, err := range e.errs {
		if et, ok := err.(ErrorTree); ok {
			t.AddTree(et.getTree())
			continue
		}
		t.Add(err.Error())
	}
	return t
}
