package integrals

import "math"

var SQRT_PI = math.Sqrt(math.Pi)

func S_ij(i, j int8, alpha, beta, ax, bx float64) float64 {
	ab_diff := ax - bx
	ab_diff_squared := math.Pow(ab_diff, 2)
	ab_sum := alpha + beta
	ab_product := alpha * beta

	if i == 0 && j == 0 {
		return SQRT_PI * math.Exp((-ab_diff_squared * ab_product / ab_sum)) / math.Sqrt(ab_sum)
	}
	if i == 0 && j == 1 {
		return SQRT_PI * ab_diff * alpha * math.Exp(-ab_diff_squared*ab_product/ab_sum) / math.Pow(ab_sum, 3./2.)
	}
	if i == 1 && j == 0 {
		return -SQRT_PI * ab_diff * beta * math.Exp(-ab_diff_squared*ab_product/ab_sum) / math.Pow(ab_sum, 3./2.)
	}
	if i == 1 && j == 1 {
		return (1. / 2.) * SQRT_PI * (-2*math.Pow(ab_diff, 2)*ab_product + ab_sum) * math.Exp(-ab_diff_squared*ab_product/ab_sum) / math.Pow(ab_sum, 5./2.)
	}
	panic("not implemented")
}
