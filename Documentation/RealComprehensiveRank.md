# RealComprehensiveRank
The `RealComprehensiveRank` function will compute a
complete case discussion of the rank of a matrix where the entries are multivariate polynomials whose indeterminants are treated as parameters. Parameters are assumed to be real valued and computation is done over polynomial equality, inequation and inequality (strict and non-strict) constraints on the parameters.

## Calling Sequence
```
RealComprehensiveRank(A, R)
RealComprehensiveRank(A, rsas, R)
RealComprehensiveRank(A, lrsas, R)
RealComprehensiveRank(A, F, R)
RealComprehensiveRank(A, F, H, R)
RealComprehensiveRank(A, F, N, P, H, R)    
```

## Input

| Variable | Description |
| --- | --- |
| `A`  | Square matrix of multivariate polynomials over `R` |
| `rsas` | Regular semi-algebraic system |
| `lrsas` | List of regular semi-algebraic systems |
| `F`  | List of polynomials over `R` representing equality constraints on the parameters |
| `F`  | List of polynomials over `R` representing equality constraints on the parameters |
| `N`  | List of polynomials over `R` representing non-negativity constraints on the parameters |
| `P`  | List of polynomials over `R` representing positivity constraints on the parameters |
| `R`  | Polynomial ring |

## Output

A list with elements of the form `[r, lrsas]` where `r` is the reank of the input matrix for all parameter values in the zero set of `lrsas`. Together, all the lists of regular semi-algebraic systems form a partition of the input constraints on the parameters.

## Example
```
with(RegularChains):
with(ParametricMatrixTools):

# Set up the matrix
R := PolynomialRing([a]):
A := Matrix([[1, 2, 3],
             [4, 5, 6],
             [7, 8, a]]):

# Compute the rank of A
result := RealComprehensiveRank(A, [], [], [a], [], R):

# Print the result
result[1][1], Display(result[1][2], R);
    2, [a - 9 = 0]

result[2][1], Display(result[2][2], R);
    3, [a > 0 and a - 9 <> 0]
```
