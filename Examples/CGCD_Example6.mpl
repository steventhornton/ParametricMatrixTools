libname := libname, currentdir():
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):
read("Examples/printResult.mpl"):

R := PolynomialRing([x, a, b, c, d]):

p1 := a*x+b:
p2 := c*x+d:

cs := GeneralConstruct([],[], R):

result, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R, 'outputType'='RS'):

printResult_CGCD(result, cs_zero, p1, p2, cs, R, 6):