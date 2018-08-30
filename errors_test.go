package hd

import (
	"fmt"
	"testing"
)

func TestErrorTree(t *testing.T) {
	var et ErrorTree
	for i := 0; i < 10; i++ {
		et.Add(fmt.Errorf("Error %d", i))
		if i%3 == 0 {
			var ett ErrorTree
			for j := 0; j < i/3+1; j++ {
				ett.Add(fmt.Errorf("Inside error %d", j))
				if j%2 == 0 {
					var ettt ErrorTree
					ettt.Name = "Some deep deep errors"
					for k := 0; k < 1+j/2; k++ {
						ettt.Add(fmt.Errorf("Deep error %d", k))
					}
					ett.Add(ettt)
				}
			}
			et.Add(ett)
		}
	}
	fmt.Println(et.Error())
}
