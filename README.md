# ParametricMatrixTools
Many common computations in linear algebra such as matrix rank, Jordan canonical form and Smith normal form, are discontinuous when a matrix contain parameters. The __ParametricMatrixTools__ package provides functions for such computations where an appropriate case discussion is provided in the answer. ParametricMatrixTools is implemented in the  [Maple](http://www.maplesoft.com/products/maple/) computer algebra system, and works along side the [LinearAlgebra](http://www.maplesoft.com/support/help/Maple/view.aspx?path=LinearAlgebra) and [RegularChains](http://regularchains.org/) packages.

The parametric matrices this package is designed to work with are matrices where the entries are multivariate polynomials in the parameters. All of the methods allow constraints to be imposed on the parameter values.

 Some methods (such as Smith form) may assume the entries are parametric univariate polynomials. A parametric univariate polynomial is a multivariate polynomial that is treated as univariate in one variable, where the coefficients of this variable are multivariate polynomials in the remaining variables, known as parameters. Computations involving parametric univariate polynomials do not allow constraints on the main variable, only on the parameters.

## Using the Package
From the `ParametricMatrixTools` directory, run `make mla`. This will create a Maple library file called `ParametricMatrixTools.mla` containing the package. See the `Examples.mw` worksheet for examples of using the functions in the package.

## Methods

- [ComprehensiveSmithForm](Documentation/ComprehensiveSmithForm.md)
- [ComprehensiveFrobeniusForm](Documentation/ComprehensiveFrobeniusForm.md)
- [JordanForm](Documentation/JordanForm.md)
- [ComprehensiveJordanForm](Documentation/ComprehensiveJordanForm.md)
- [ComprehensiveGcd](Documentation/ComprehensiveGcd.md)
- [ListComprehensiveGcd](Documentation/ListComprehensiveGcd.md)
- [ComprehensiveSquareFreeFactorization](Documentation/ComprehensiveSquareFreeFactorization.md)

Detailed descriptions of each of these methods can be found in the [Documentation](Documentation) directory. An Maple worksheet [Examples.mw](Examples.mw) contains numerous examples of each of the above methods in use.

__This package is being actively developed and likely contains bugs :bug:__
