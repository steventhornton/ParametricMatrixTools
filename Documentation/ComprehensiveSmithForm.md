# ComprehensiveSmithForm
The `ComprehensiveSmithForm` function will compute a complete case discussion of the [Smith normal form](https://en.wikipedia.org/wiki/Smith_normal_form) of a matrix of monic, parametric univariate polynomials.

## Calling Sequence
```
ComprehensiveSmithForm(A, v, R, options)
ComprehensiveSmithForm(A, v, rs, R, options)
ComprehensiveSmithForm(A, v, cs, R, options)
ComprehensiveSmithForm(A, v, F, R, options)
ComprehensiveSmithForm(A, v, F, H, R, options)
```

## Input
| Variable | Description |
| --- | --- |
|`A` | Square matrix of monic, parametric univariate polynomials in `v` |
|`v` | The variable that the entries of `A` are monic polynomials in |
|`rs`| Regular system |
|`cs`| Constructible set |
|`F` | List of polynomials over `R` representing equality constraints on the parameters |
|`H` | List of polynomials over `R` representing inequation constraints on the parameters |
|`R` | Polynomial ring where `v` is the greatest variable |

## Options
| Option name | Default value | Description |
| --- | --- | --- |
| `'outputType'` | `'CS'` | When set to `'CS'` or `'ConstructibleSet'` the output will contain constructible sets, when set to `'RS'` or `'RegularSystem'`, the output will contain regular systems. |
| `'outputMatrices'` | `'S'` | Which matrices to return, value can be either `'S'`, `'U'`, or `'V'` or a list containing any combination of `'S'`, `'U'`, and `'V'` |
| `'lazy'` | `false` | When `true`, this method will only return one branch of the computation. This option is not currently implemented. |
| `'algorithm'` | `'minors'` | Specify which algorithm to use, currently the `'minors'` algorithm is the only one implemented. |

## Option Compatibility
- Only the Smith form can be returned when the `minors` algorithm is used. The transformation matrices can not be computed with this algorithm.

## Output
A list with elements in one of the following forms:

| Output structure | Options |
| --- | --- |
| `[S, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = 'S'` |
| `[S, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = 'S'` |
| `[U, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = 'U'` |
| `[U, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = 'U'` |
| `[V, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = 'V'` |
| `[V, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = 'V'` |
| `[S, U, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = ['S', 'U']` |
| `[S, U, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = ['S', 'U']` |
| `[S, V, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = ['S', 'V']` |
| `[S, V, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = ['S', 'V']` |
| `[S, U, V, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = ['S', 'U', 'V']` |
| `[S, U, V, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = ['S', 'U', 'V']` |
| `[U, V, rs]` | `'outputType'` is `'RegularSystem'` or ``'RS'`` and `'outputMatrices' = ['U', 'V']` |
| `[U, V, cs]` | `'outputType'` is `'ConstructibleSet'` or ``'CS'`` and `'outputMatrices' = ['U', 'V']` |

Where `S` is the Smith form of the input matrix for all parameter values in the zero set of `cs` or `rs`. `U` and `V` are the left and right transformation matrices respectively such that `S = U.A.V`. Together, all the constructible sets or regular systems form a partition of the input constructible set, regular system or set of equations and inequations.

## Assumptions
1. `v` is the greatest variable of `R`
2. `A` is a square matrix
3. All entries of `A` are monic in `v`
4. `F` and `H` do not contain any polynomials containing `v`
5. If a constructible set or regular system is input, v cannot be a variable in any equations or inequations

## Example
```
with(RegularChains):
with(ParametricMatrixTools):

# Set up the matrix
R := PolynomialRing([x,a]):
A := Matrix([[ 0,  a,  1,  1,  1],
             [ 2, -2,  0, -2, -4],
             [ 0,  0,  1,  1,  3],
             [-6,  0, -3, -1, -3],
             [ 2,  2,  2,  2,  4]]):
Ax := x*IdentityMatrix(5) - A:

# Compute the Smith form of Ax
result := ComprehensiveSmithForm(Ax, x, [], R):

# Print the result
map(factor, result[1][1]), Display(result[1][2], R);
    [1,0,0,0,0]
    [0,1,0,0,0]                              {a-4 <> 0
    [0,0,1,0,0]                              {a-1 <> 0
    [0,0,0,1,0]                              {23a^2-61a+68 <> 0
    [0,0,0,0,-(x-2)(x^2+2)(-x^2+2a-4)]

map(factor, result[2][1]), Display(result[2][2], R);
    [1,0,0,0,0]
    [0,1,0,0,0]
    [0,0,1,0,0]                              {23a^2-61a+68 = 0
    [0,0,0,1,0]
    [0,0,0,0,-(x-2)(x^2+2)(-x^2+2a-4)]

map(factor, result[3][1]), Display(result[3][2], R);
    [1,0,0,0,  0]
    [0,1,0,0,  0]
    [0,0,1,0,  0]                     {a-4 = 0
    [0,0,0,x-2,0]
    [0,0,0,0,  (x-2)*(x+2)*(x^2+2)]

map(factor, result[4][1]), Display(result[4][2], R);
    [1,0,0,0,    0]
    [0,1,0,0,    0]
    [0,0,1,0,    0]                     {a-1 = 0
    [0,0,0,x^2+2,0]
    [0,0,0,0,    (x-2)*(x^2+2)]
```
