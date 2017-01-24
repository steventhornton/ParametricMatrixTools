libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a, b]):

p := (x-b)^4+a*x^2:

cs := GeneralConstruct([],[],R):

result := ComprehensiveSquareFreeFactorization(p, x, cs, R, 'outputType'='RS'):

printResult_CSFF(result, p, cs, R, 3);