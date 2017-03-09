libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a,b]):
rs := RegularSystem(R):

p := (x+a)*(x+b)^2*(x+1)^3:

# Compute the square-free factorization of p
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS'):
printResult_SFFm(result, p, rs, R, false, 2);