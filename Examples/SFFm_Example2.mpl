libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a,b]):
rs := RegularSystem(R):

p := (x+a)*(x+b)^2*(x+1)^3:

# Compute the square-free factorization of p in the sense of Kalkbrener
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS', 'output'='kalkbrener'):
printResult_SFFm(result, p, rs, R, false, 2);

# Compute the square-free factorization of p in the sense of Kalkbrener
printf("\n");
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS', 'output'='lazard'):
printResult_SFFm(result, p, rs, R, true, 2);