libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a]):

p1 := 2:
p2 := x + a:

cs := GeneralConstruct([],[], R):

result, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS'):

printResult_CGCD(result, cs_zero, p1, p2, cs, R, 8):