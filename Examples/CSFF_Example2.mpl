libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a]):

p := (x+1)^2 * (x+a):

cs := GeneralConstruct([],[],R):

result := ComprehensiveSquareFreeFactorization(p, x, cs, R, 'outputType'='RS'):

printResult_CSFF(result, p, cs, R, 2);