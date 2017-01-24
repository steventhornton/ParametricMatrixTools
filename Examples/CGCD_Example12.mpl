libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a]):

p1 := (a+1)*x + (a+1):
p2 := (a+1):

cs := GeneralConstruct([],[], R):

result, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS'):

printResult_CGCD(result, cs_zero, p1, p2, cs, R, 12):