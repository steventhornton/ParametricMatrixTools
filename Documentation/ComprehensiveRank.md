# ComprehensiveRank
The `ComprehensiveRank` function will compute a complete case discussion of the rank of a matrix of multivariate polynomials.

## Calling Sequence
```
ComprehensiveRank(A, R, options)
ComprehensiveRank(A, rs, R, options)
ComprehensiveRank(A, cs, R, options)
ComprehensiveRank(A, F, R, options)
ComprehensiveRank(A, F, H, R, options)
```

## Input

| Variable | Description |
| --- | --- |
| `A`  | Square matrix of multivariate polynomials over `R` |
| `rs` | Regular system |
| `cs` | Constructible set |
| `F`  | List of polynomials over `R` representing equality constraints on the parameters |
| `H`  | List of polynomials over `R` representing inequation constraints on the parameters |
| `R`  | Polynomial ring |

## Options

| Option name | Default value | Description |
| --- | --- | --- |
| `'outputType'` | `'CS'` | When set to `'CS'` or `'ConstructibleSet'` the output will contain constructible sets, when set to `'RS'` or `'RegularSystem'`, the output will contain regular systems. |
## Output

| Output structure | Options |
| --- | --- |
| `[r, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` |
| `[r, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` |

Where `r` is the rank of the input matrix for all parameter values in the zero set of `cs` or `rs`. Together, all the constructible sets or regular systems form a partition of the input constructible set, regular system or set of equations and inequations.

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
result := ComprehensiveRank(A, R):

# Print the result
result[1][1], Display(result[1][2], R);
    2, a - 9 = 0

result[2][1], Display(result[2][2], R);
    3, a - 9 <> 0
```
