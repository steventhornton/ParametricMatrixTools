libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):


R := PolynomialRing([x,a]):
rs := RegularSystem(R):

p := (x+1)^2*(x^2+x+1)*(x+a):

# Compute the square-free factorization of p in the sense of Kalkbrener
result := SquarefreeFactorization_monic(p, x, rs, R, false):
printResult_SFFm(result, p, rs, R, false, 4);

# Compute the square-free factorization of p in the sense of Kalkbrener
printf("\n");
result := SquarefreeFactorization_monic(p, x, rs, R, true):
printResult_SFFm(result, p, rs, R, true, 4);