package integrals

import (
	"gonum.org/v1/gonum/mat"
	"math"
	"scientificgo.org/special"
)

var TwoPi = PI * 2

func boys(n, t float64) float64 {
	return special.HypPFQ([]float64{n + 0.5}, []float64{n + 1.5}, -t) / (2.0*n + 1.0)
}

func v_ij(i, j, k, l, m, n int8, alpha, beta float64, Aa, Bb, Cc []float64) float64 {
	A := mat.NewVecDense(3, Aa)
	B := mat.NewVecDense(3, Bb)
	C := mat.NewVecDense(3, Cc)
	var AB *mat.VecDense
	var alphaA *mat.VecDense
	//alphaA.CopyVec(A)
	var betaB *mat.VecDense
	//betaB.CopyVec(B)
	AB.SubVec(A, B)
	r_AB := mat.Dot(AB, AB)
	alphaA.ScaleVec(alpha, A)
	betaB.ScaleVec(beta, B)

	p := alpha + beta
	q := alpha * beta

	var Pp *mat.VecDense
	Pp.AddVec(alphaA, betaB)
	var P *mat.VecDense
	P.ScaleVec(1/p, Pp)
	var PC *mat.VecDense
	PC.SubVec(P, C)
	p_RPC := p * mat.Dot(PC, PC)

	//AB := A - B
	//r_AB := np.dot(AB, AB)
	//P = (alpha*A + beta*B) / p
	//PC = P - C
	//p_RPC = p * np.dot(PC, PC)
	A_x, A_y, A_z := A.AtVec(0), A.AtVec(1), A.AtVec(2)
	B_x, B_y, B_z := B.AtVec(0), B.AtVec(1), B.AtVec(2)
	C_x, C_y, C_z := C.AtVec(0), C.AtVec(1), C.AtVec(2)

	if i == 0 && j == 0 && k == 0 && l == 0 && m == 0 && n == 0 {
		return TwoPi * math.Exp(-q*r_AB/p) * boys(0, p_RPC) / p
	}
	if i == 0 && j == 0 && k == 0 && l == 0 && m == 0 && n == 1 {
		return TwoPi * (-alpha*(-A_z+B_z)*boys(0, p_RPC) + (-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 0 && k == 0 && l == 0 && m == 1 && n == 0 {
		return TwoPi * (-alpha*(-A_y+B_y)*boys(0, p_RPC) + (-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 0 && k == 0 && l == 1 && m == 0 && n == 0 {
		return TwoPi * (-alpha*(-A_x+B_x)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 0 && k == 1 && l == 0 && m == 0 && n == 0 {
		return TwoPi * (-beta*(A_z-B_z)*boys(0, p_RPC) + (-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 1 && k == 0 && l == 0 && m == 0 && n == 0 {
		return TwoPi * (-beta*(A_y-B_y)*boys(0, p_RPC) + (-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 1 && j == 0 && k == 0 && l == 0 && m == 0 && n == 0 {
		return TwoPi * (-beta*(A_x-B_x)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 0 && k == 1 && l == 0 && m == 0 && n == 1 {
		return PI * (-2*alpha*(-A_z+B_z)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC)/p - 2*beta*(A_z-B_z)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC)/p + boys(0, p_RPC) - boys(1, p_RPC) + 2*q*(-A_z+B_z)*(A_z-B_z)*boys(0, p_RPC)/p + 2*math.Pow((-A_z*alpha-B_z*beta+C_z*p), 2)*boys(2, p_RPC)/p) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 0 && k == 1 && l == 0 && m == 1 && n == 0 {
		return TwoPi * (-alpha*(-A_y+B_y)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC) - beta*(A_z-B_z)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC) + q*(-A_y+B_y)*(A_z-B_z)*boys(0, p_RPC) + (-A_y*alpha-B_y*beta+C_y*p)*(-A_z*alpha-B_z*beta+C_z*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 0 && j == 0 && k == 1 && l == 1 && m == 0 && n == 0 {
		return TwoPi * (-alpha*(-A_x+B_x)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC) - beta*(A_z-B_z)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC) + q*(-A_x+B_x)*(A_z-B_z)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*(-A_z*alpha-B_z*beta+C_z*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 0 && j == 1 && k == 0 && l == 0 && m == 0 && n == 1 {
		return TwoPi * (-alpha*(-A_z+B_z)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC) - beta*(A_y-B_y)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC) + q*(A_y-B_y)*(-A_z+B_z)*boys(0, p_RPC) + (-A_y*alpha-B_y*beta+C_y*p)*(-A_z*alpha-B_z*beta+C_z*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 0 && j == 1 && k == 0 && l == 0 && m == 1 && n == 0 {
		return PI * (-2*alpha*(-A_y+B_y)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC)/p - 2*beta*(A_y-B_y)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC)/p + boys(0, p_RPC) - boys(1, p_RPC) + 2*q*(-A_y+B_y)*(A_y-B_y)*boys(0, p_RPC)/p + 2*math.Pow((-A_y*alpha-B_y*beta+C_y*p), 2)*boys(2, p_RPC)/p) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	if i == 0 && j == 1 && k == 0 && l == 1 && m == 0 && n == 0 {
		return TwoPi * (-alpha*(-A_x+B_x)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC) - beta*(A_y-B_y)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC) + q*(-A_x+B_x)*(A_y-B_y)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*(-A_y*alpha-B_y*beta+C_y*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 1 && j == 0 && k == 0 && l == 0 && m == 0 && n == 1 {
		return TwoPi * (-alpha*(-A_z+B_z)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC) - beta*(A_x-B_x)*(-A_z*alpha-B_z*beta+C_z*p)*boys(1, p_RPC) + q*(A_x-B_x)*(-A_z+B_z)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*(-A_z*alpha-B_z*beta+C_z*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 1 && j == 0 && k == 0 && l == 0 && m == 1 && n == 0 {
		return TwoPi * (-alpha*(-A_y+B_y)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC) - beta*(A_x-B_x)*(-A_y*alpha-B_y*beta+C_y*p)*boys(1, p_RPC) + q*(A_x-B_x)*(-A_y+B_y)*boys(0, p_RPC) + (-A_x*alpha-B_x*beta+C_x*p)*(-A_y*alpha-B_y*beta+C_y*p)*boys(2, p_RPC)) * math.Exp(-q*r_AB/p) / math.Pow(p, 3)
	}
	if i == 1 && j == 0 && k == 0 && l == 1 && m == 0 && n == 0 {
		return PI * (-2*alpha*(-A_x+B_x)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC)/p - 2*beta*(A_x-B_x)*(-A_x*alpha-B_x*beta+C_x*p)*boys(1, p_RPC)/p + boys(0, p_RPC) - boys(1, p_RPC) + 2*q*(-A_x+B_x)*(A_x-B_x)*boys(0, p_RPC)/p + 2*math.Pow((-A_x*alpha-B_x*beta+C_x*p), 2)*boys(2, p_RPC)/p) * math.Exp(-q*r_AB/p) / math.Pow(p, 2)
	}
	panic("not implemented")
}
