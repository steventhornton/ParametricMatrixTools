libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a]):
rs := RegularSystem(R):

p := (x+1)^2*(x+a):

# Compute the square-free factorization of p in the sense of Kalkbrener
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS', 'output'='kalkbrener'):
printResult_SFFm(result, p, rs, R, false, 1);

# Compute the square-free factorization of p in the sense of Kalkbrener
printf("\n");
result := SquareFreeFactorization_monic(p, x, rs, R, 'outputType'='RS', 'output'='lazard'):
printResult_SFFm(result, p, rs, R, true, 1);