package integrals

import (
	"gonum.org/v1/gonum/mat"
	"math"
)

var Pi5_2 = math.Pow(PI, 5/2)

func VElecIJ(ii, jj, kk, ll, mm, nn, oo, pp, qq, rr, ss, tt int8, alpha, beta, gamma, delta float64, Aa, Bb, Cc, Dd [3]float64) float64 {
	p := alpha + beta
	q := gamma + delta
	rho := p * q / (p + q)

	A := mat.NewVecDense(3, Aa[:])
	B := mat.NewVecDense(3, Bb[:])
	C := mat.NewVecDense(3, Cc[:])
	D := mat.NewVecDense(3, Dd[:])
	AB := mat.NewVecDense(3, make([]float64, 3))
	AB.SubVec(A, B)
	CD := mat.NewVecDense(3, make([]float64, 3))
	CD.SubVec(C, D)
	r_AB := mat.Dot(AB, AB)
	r_CD := mat.Dot(CD, CD)
	alphaA := mat.NewVecDense(3, make([]float64, 3))
	betaB := mat.NewVecDense(3, make([]float64, 3))
	alphaA.ScaleVec(alpha, A)
	betaB.ScaleVec(beta, B)
	gammaC := mat.NewVecDense(3, make([]float64, 3))
	deltaD := mat.NewVecDense(3, make([]float64, 3))
	gammaC.ScaleVec(gamma, C)
	deltaD.ScaleVec(delta, D)

	Pp := mat.NewVecDense(3, make([]float64, 3))
	Qq := mat.NewVecDense(3, make([]float64, 3))
	Pp.AddVec(alphaA, betaB)
	Qq.AddVec(gammaC, deltaD)

	P := mat.NewVecDense(3, make([]float64, 3))
	Q := mat.NewVecDense(3, make([]float64, 3))
	P.ScaleVec(1/p, Pp)
	Q.ScaleVec(1/q, Qq)

	PQ := mat.NewVecDense(3, make([]float64, 3))
	PQ.SubVec(P, Q)
	r_PQ := mat.Dot(PQ, PQ)

	A_x, A_y, A_z := A.AtVec(0), A.AtVec(1), A.AtVec(2)
	B_x, B_y, B_z := B.AtVec(0), B.AtVec(1), B.AtVec(2)
	C_x, C_y, C_z := C.AtVec(0), C.AtVec(1), C.AtVec(2)
	D_x, D_y, D_z := D.AtVec(0), D.AtVec(1), D.AtVec(2)

	if ii == 0 && jj == 0 && kk == 0 && ll == 0 && mm == 0 && nn == 0 && oo == 0 && pp == 0 && qq == 0 && rr == 0 && ss == 0 && tt == 0 {
		return 2 * Pi5_2 * math.Exp(-alpha*beta*r_AB/p-delta*gamma*r_CD/q) * boys(0, (math.Pow((p*(C_x*gamma+D_x*delta)-q*(A_x*alpha+B_x*beta)), 2)+math.Pow((p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta)), 2)+math.Pow((p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta)), 2))/(p*q*(p+q))) / (p * q * math.Sqrt(p+q))
	}
	if ii == 0 && jj == 0 && kk == 0 && ll == 0 && mm == 0 && nn == 0 && oo == 0 && pp == 0 && qq == 0 && rr == 0 && ss == 0 && tt == 1 {
		return (1 / 2) * (-2*Pi5_2*delta*gamma*(-2*C_z+2*D_z)*math.Exp(-alpha*beta*r_AB/p-delta*gamma*r_CD/q)*boys(0, (math.Pow((p*(C_x*gamma+D_x*delta)-q*(A_x*alpha+B_x*beta)), 2)+math.Pow((p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta)), 2)+math.Pow((p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta)), 2))/(p*q*(p+q)))/(p*math.Pow(q, 2)*math.Sqrt(p+q)) - 4*Pi5_2*delta*(p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta))*math.Exp(-alpha*beta*r_AB/p-delta*gamma*r_CD/q)*boys(1, (math.Pow((p*(C_x*gamma+D_x*delta)-q*(A_x*alpha+B_x*beta)), 2)+math.Pow((p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta)), 2)+math.Pow((p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta)), 2))/(p*q*(p+q)))/(p*math.Pow(q, 2)*math.Pow((p+q), 3/2))) / delta
	}
	if ii == 0 && jj == 0 && kk == 0 && ll == 0 && mm == 0 && nn == 0 && oo == 0 && pp == 0 && qq == 0 && rr == 0 && ss == 1 && tt == 0 {
		return (1 / 2) * (-2*Pi5_2*delta*gamma*(-2*C_y+2*D_y)*math.Exp(-alpha*beta*r_AB/p-delta*gamma*r_CD/q)*boys(0, (math.Pow((p*(C_x*gamma+D_x*delta)-q*(A_x*alpha+B_x*beta)), 2)+math.Pow((p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta)), 2)+math.Pow((p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta)), 2))/(p*q*(p+q)))/(p*math.Pow(q, 2)*math.Sqrt(p+q)) - 4*Pi5_2*delta*(p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta))*math.Exp(-alpha*beta*r_AB/p-delta*gamma*r_CD/q)*boys(1, (math.Pow((p*(C_x*gamma+D_x*delta)-q*(A_x*alpha+B_x*beta)), 2)+math.Pow((p*(C_y*gamma+D_y*delta)-q*(A_y*alpha+B_y*beta)), 2)+math.Pow((p*(C_z*gamma+D_z*delta)-q*(A_z*alpha+B_z*beta)), 2))/(p*q*(p+q)))/(p*math.Pow(q, 2)*math.Pow((p+q), 3/2))) / delta
	}
	panic("not implemented")
}
