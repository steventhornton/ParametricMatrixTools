# ParametricMatrixTools
Many common computations in linear algebra such as matrix rank, Jordan canonical form and Smith normal form, are discontinuous when a matrix contain parameters. The __ParametricMatrixTools__ package provides functions for such computations where an appropriate case discussion is provided in the answer. ParametricMatrixTools is implemented in the  [Maple](http://www.maplesoft.com/products/maple/) computer algebra system, and works along side the [LinearAlgebra](http://www.maplesoft.com/support/help/Maple/view.aspx?path=LinearAlgebra) and [RegularChains](http://regularchains.org/) packages.

The parametric matrices this package is designed to work with are matrices where the entries are multivariate polynomials in the parameters. All of the methods allow constraints to be imposed on the parameter values.

 Some methods (such as Smith form) may assume the entries are parametric univariate polynomials. A parametric univariate polynomial is a multivariate polynomial that is treated as univariate in one variable, where the coefficients of this variable are multivariate polynomials in the remaining variables, known as parameters. Computations involving parametric univariate polynomials do not allow constraints on the main variable, only on the parameters.

## Using the Package
From the `ParametricMatrixTools` directory, run `make mla`. This will create a Maple library file called `ParametricMatrixTools.mla` containing the package. See the `Examples.mw` worksheet for examples of using the functions in the package.

## Methods
The following methods are currently implemented:

- ComprehensiveGcd

Brief descriptions of each of these methods can be found in the headers of each file. An example Maple worksheet has been provided to show how these methods are used.

## ComprehensiveGcd
The `ComprehensiveGcd` function will compute a complete case discussion of the gcd of two parametric univariate polynomials.

### Example
```
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):

# Set up the polynomials with variable order x > a
R := PolynomialRing([x, a]):
p1 := (x+1)*(x+a):
p2 := (x+2)^2:

# Constraints on the parameters (a <> 10)
cs := GeneralConstruct([], [a-10], R):

# Compute the gcd
g, cs_zero := ComprehensiveGcd(p1, p2, x, cs, R):
IsEmpty(cs_zero, R);
    true
Display(g, R);
    [[-a*x-2*a+2*x+4, a-2 <> 0], [x^2+4*x+4, a-2 = 0]]
```


__This package is being actively developed and likely contains bugs :bug:__
