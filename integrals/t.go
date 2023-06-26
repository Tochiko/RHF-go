package integrals

import "math"

func T_ij(i, j int8, alpha, beta, ax, bx float64) float64 {
	ab_diff := ax - bx
	ab_diff_squared := math.Pow(ab_diff, 2)
	ab_sum := alpha + beta
	ab_product := alpha * beta
	ab_diff_ab_product := 2 * math.Pow(ab_diff, 2) * ab_product
	ab_diff_sq_ab_product_per_sum := -ab_diff_squared * ab_product / ab_sum

	if i == 0 && j == 0 {
		return -SQRT_PI*beta*(1+alpha*(ab_diff_ab_product-ab_sum)/math.Pow(ab_sum, 2))*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Sqrt(ab_sum) + SQRT_PI*beta*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Sqrt(ab_sum)
	}
	if i == 0 && j == 1 {
		return -SQRT_PI*ab_diff*alpha*beta*(3+alpha*(ab_diff_ab_product-3*alpha-3*beta)/math.Pow(ab_sum, 2))*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 3/2) + 3*SQRT_PI*ab_diff*alpha*beta*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 3/2)
	}
	if i == 1 && j == 0 {
		return -SQRT_PI*ab_diff*math.Pow(beta, 2)*(-1+alpha*(-ab_diff_ab_product+3*alpha+3*beta)/math.Pow(ab_sum, 2))*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 3/2) - SQRT_PI*ab_diff*math.Pow(beta, 2)*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 3/2)
	}
	if i == 1 && j == 1 {
		return (3/2)*SQRT_PI*beta*(-ab_diff_ab_product+ab_sum)*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 5/2) - 1/2*SQRT_PI*beta*(-3*ab_diff_ab_product+3*alpha+3*beta+alpha*(3*ab_diff_ab_product*ab_sum-ab_diff_ab_product*(ab_diff_ab_product-3*alpha-3*beta)-3*math.Pow(ab_sum, 2))/math.Pow(ab_sum, 2))*math.Exp(ab_diff_sq_ab_product_per_sum)/math.Pow(ab_sum, 5/2)
	}
	panic("not implemented")
}
