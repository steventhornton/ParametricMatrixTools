# ComprehensiveSquareFreeFactorization
The `ComprehensiveSquareFreeFactorization` function will compute a complete case discussion of the [square-free factorization](https://en.wikipedia.org/wiki/Square-free_polynomial) of a parametric univariate polynomial.

## Calling Sequence
```
ComprehensiveSquareFreeFactorization(p, v, R, options)
ComprehensiveSquareFreeFactorization(p, v, rs, R, options)
ComprehensiveSquareFreeFactorization(p, v, cs, R, options)
ComprehensiveSquareFreeFactorization(p, v, F, R, options)
ComprehensiveSquareFreeFactorization(p, v, F, H, R, options)
```

## Input

| Variable | Description |
| --- | --- |
| `p` | A parametric, univariate polynomial in `v` |
| `v`  | The variable that `p1` and `p2` are polynomials in |
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

A list with elements of the form:
| Output structure | Options |
| --- | --- |
| `[lp, rs_i]` | `'outputType'` is `'RegularSystem'` or `'RS'` |
| `[lp, cs_i]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` |

Where `lp` is a list with elements of the form `[q_j, m]` such that `p = product(q_j^m)` and `q_j` are the square-free factoris in the zero set of `cs_i` or `rs_i`.

## Assumptions
1. `degree(p, mvar(R)) > 0`

## Example
```
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):

# Set up the polynomials with variable order x > a
R := PolynomialRing([x, a]):
p := (x+1)^2 * (x+1):

# Constraints on the parameters (a <> 10)
cs := GeneralConstruct([], [a-10], R):

# Compute the square-free factorization
result := ComprehensiveSquareFreeFactorization(p, x, cs, R)
nops(result);
    2
Display(result[1], R);
    [[x+a, 1], [x+1, 2]], a-1 <> 0, a-10 <> 0]
Display(result[2], R);
    [[[x+1, 3]], a-1 = 0]
```
