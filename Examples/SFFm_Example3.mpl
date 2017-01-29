libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a,b]):
rs := RegularSystem(R):

p := (x+a)*(a*b+x)^2:

# Compute the square-free factorization of p in the sense of Kalkbrener
result := SquarefreeFactorization_monic(p, x, rs, R, false):
printResult_SFFm(result, p, rs, R, false, 3);

# Compute the square-free factorization of p in the sense of Kalkbrener
printf("\n");
result := SquarefreeFactorization_monic(p, x, rs, R, true):
printResult_SFFm(result, p, rs, R, true, 3);