package integrals

import (
	"scientificgo.org/special"
)

func boys(n, t float64) float64 {
	return special.HypPFQ([]float64{n + 0.5}, []float64{n + 1.5}, -t) / (2.0*n + 1.0)
}
