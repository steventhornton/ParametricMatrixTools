# ComprehensiveFrobeniusForm
The `ComprehensiveFrobeniusForm` function will compute a complete case discussion of the [Frobenius (rational) normal form](https://en.wikipedia.org/wiki/Frobenius_normal_form) of a matrix of multivariate polynomials.

## Calling Sequence
```
ComprehensiveFrobeniusForm(A, R, options)
ComprehensiveFrobeniusForm(A, rs, R, options)
ComprehensiveFrobeniusForm(A, cs, R, options)
ComprehensiveFrobeniusForm(A, F, R, options)
ComprehensiveFrobeniusForm(A, F, H, R, options)
```

## Input

| Variable | Description |
| - | - |
| `A`  | Square matrix of multivariate polynomials over `R` |
| `rs` | Regular system |
| `cs` | Constructible set |
| `F`  | List of polynomials over `R` representing equality constraints on the parameters |
| `H`  | List of polynomials over `R` representing inequation constraints on the parameters |
| `R`  | Polynomial ring |

## Options

| Option name | Default value | Description |
| - | - | - |
| `'outputType'` | `'CS'` | When set to `'CS'` or `'ConstructibleSet'` the output will contain constructible sets, when set to `'RS'` or `'RegularSystem'`, the output will contain regular systems. |
| `'outputMatrices'` | `'F'` | Which matrices to return, value can be either `'F'` or `'Q'` or a list containing any combination of `'F'`, and `'Q'` |
| `'lazy'` | `false` | When `true`, this method will only return one branch of the computation. This option is not currently implemented. |
| `'algorithm'` | `'snf_minors'` | Specify which algorithm to use, currently the `'snf_minors'` algorithm is the only one implemented. |

## Option Compatibility
- Only the Frobenius form can be returned when the `snf_minors` algorithm is used. The transformation matrix can not be computed with this algorithm.

## Output

| Output structure | Options |
| - | - |
| `[F, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = 'F'` |
| `[F, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = 'F'` |
| `[Q, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = 'Q'` |
| `[Q, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = 'Q'` |
| `[F, Q, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'outputMatrices' = ['F', 'Q']` |
| `[F, Q, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` and `'outputMatrices' = ['F', 'Q']` |

Where `F` is the Frobenius form of the input matrix for all parameter values in the zero set of `cs` or `rs`. `Q` is the similarity transformation matrix  such that `F = Q^(-1) A Q`. Together, all the constructible sets or regular systems form a partition of the input constructible set, regular system or set of equations and inequations.

## Example
```
with(RegularChains):
with(ParametricMatrixTools):

# Set up the matrix
R := PolynomialRing([a]):
A := Matrix([[ 0,  a,  1,  1,  1],
             [ 2, -2,  0, -2, -4],
             [ 0,  0,  1,  1,  3],
             [-6,  0, -3, -1, -3],
             [ 2,  2,  2,  2,  4]]):

# Compute the Frobenius form of A
result := ComprehensiveFrobeniusForm(A, [], R)

# Print the result
result[1][1], Display(result[1][2], R);
    [0, 0, 0, 0, -8a + 16]
    [1, 0, 0, 0,  4a -  8]          {a-4 <> 0
    [0, 1, 0, 0, -4a + 12]          {a-1 <> 0
    [0, 0, 1, 0,  2a -  6]          {23*a^2-61*a+68 <> 0
    [0, 0, 0, 1,        2]

result[2][1], Display(result[2][2], R);
    [0, 0, 0, 0, -8a + 16]
    [1, 0, 0, 0,  4a -  8]
    [0, 1, 0, 0, -4a + 12]          {23*a^2-61*a+68 = 0
    [0, 0, 1, 0,  2a -  6]
    [0, 0, 0, 1,        2]

result[3][1], Display(result[3][2], R);
    [0, 0, 0, 8, 0]
    [1, 0, 0, 0, 0]
    [0, 1, 0, 2, 0]                 {a-4 = 0
    [0, 0, 1, 0, 0]
    [0, 0, 0, 0, 2]

result[4][1], Display(result[4][2], R);
    [0, 0,  4, 0,  0]
    [1, 0, -2, 0,  0]
    [0, 1,  2, 0,  0]                 {a-1 = 0
    [0, 0,  0, 0, -2]
    [0, 0,  0, 1,  0]
```
