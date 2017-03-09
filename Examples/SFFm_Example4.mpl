libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a]):
rs := RegularSystem(R):

p := (x+1)^2*(x^2+x+1)*(x+a):

# Compute the square-free factorization of p
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS'):
printResult_SFFm(result, p, rs, R, false, 4);