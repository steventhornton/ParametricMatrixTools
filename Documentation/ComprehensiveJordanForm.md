# ComprehensiveJordanForm
The `ComprehensiveJordanForm` function will compute a complete case discussion of the [Jordan normal form](https://en.wikipedia.org/wiki/Jordan_normal_form) of a matrix of multivariate polynomials.

## Calling Sequence
```
ComprehensiveJordanForm(A, R, options)
ComprehensiveJordanForm(A, rs, R, options)
ComprehensiveJordanForm(A, cs, R, options)
ComprehensiveJordanForm(A, F, R, options)
ComprehensiveJordanForm(A, F, H, R, options)
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
| `'outputMatrices'` | `'J'` | Which matrices to return, value can be either `'J'` or `'Q'` or a list containing any combination of `'J'`, and `'Q'` |
| `'lazy'` | `false` | When `true`, this method will only return one branch of the computation. This option is not currently implemented. |
| `'algorithm'` | `'snf_minors'` | Specify which algorithm to use, currently the `'snf_minors'` algorithm is the only one implemented. |

## Option Compatibility
- Only the Jordan form can be returned when the `snf_minors` algorithm is used. The transformation matrix can not be computed with this algorithm.

## Output

| Output structure | Options |
| --- | --- |
| `[J, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = 'J'` |
| `[J, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = 'J'` |
| `[Q, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = 'Q'` |
| `[Q, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = 'Q'` |
| `[J, Q, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = ['J', 'Q']` |
| `[J, Q, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = ['J', 'Q']` |

Where `J` is the Jordan form of the input matrix for all parameter values in the zero set of `cs` or `rs`. `Q` is the similarity transformation matrix  such that `J = Q^(-1) A Q`. Together, all the constructible sets or regular systems form a partition of the input constructible set, regular system or set of equations and inequations.

## Example
```
with(RegularChains):
with(ParametricMatrixTools):

# Set up the matrix
R := PolynomialRing([a]):
A := Matrix([[1, a, 0, 0, 0],
             [0, 1, a, 0, 0],
             [0, 0, 1, a, 0],
             [0, 0, 0, 1, a],
             [0, 0, 0, 0, 1]]:

# Compute the Jordan form of A
result := ComprehensiveJordanForm(A, R):

# Print the result
result[1][1], Display(result[1][2], R);
    [1, 0, 0, 0, 0]
    [0, 1, 0, 0, 0]
    [0, 0, 1, 0, 0]             {a = 0
    [0, 0, 0, 1, 0]
    [0, 0, 0, 0, 0]

result[2][1], Display(result[2][2], R);
    [1, 1, 0, 0, 0]
    [0, 1, 1, 0, 0]
    [0, 0, 1, 1, 0]             {a <> 0
    [0, 0, 0, 1, 1]
    [0, 0, 0, 0, 1]
```
