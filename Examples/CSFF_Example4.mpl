libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a, b]):

p := a*x^2+x^3-a*b-a-x:

cs := GeneralConstruct([],[],R):

result := ComprehensiveSquareFreeFactorization(p, x, cs, R, 'outputType'='RS'):

printResult_CSFF(result, p, cs, R, 4);