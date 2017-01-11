# ListComprehensiveGcd
The `ListComprehensiveGcd` function will compute a complete case discussion of the [greatest common divisor](https://en.wikipedia.org/wiki/Polynomial_greatest_common_divisor) of a list of parametric univariate polynomials.

## Calling Sequence
```
ListComprehensiveGcd(lp, v, R, options)
ListComprehensiveGcd(lp, v, rs, R, options)
ListComprehensiveGcd(lp, v, cs, R, options)
ListComprehensiveGcd(lp, v, F, R, options)
ListComprehensiveGcd(lp, v, F, H, R, options)
```

## Input

| Variable | Description |
| --- | --- |
| `lp` | A list of parametric, univariate polynomial in `v` |
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

A sequence `gcdList, cs_zero` where `gcdList` is a list with elements of the form:

| Output structure | Options |
| --- | --- |
| `[g, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` |
| `[g, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` |
| `[g, cof_p1, cof_p2, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` and `'cofactors' = true` |

Where g_i is the gcd of p1 and p2 for all values in the zero set of `cs` or `rs`, and `cof_p1` and `cof_p2` are the cofactors of `p1` and `p2` respectively: `cof_p1 = p1/g`, `cof_p2 = p2/g`. `cs_zero` is the constructible set where both `p1` and `p2` vanish for all values in its zero set. The set of constructible sets returned (including `cs_zero`) form a partition of the input constructible set.

## Assumptions
1. `v` is the greatest variable of `R`

## Example
```
with(RegularChains):
with(ConstructibleSetTools):
with(ParametricMatrixTools):

# Set up the polynomials with variable order x > a
R := PolynomialRing([x, a]):
lp := [x+1, (x+a)*(x-a), x-a]:

# Compute the gcd
g, cs_zero := ListComprehensiveGcd(lp, x, [], R);

# Print the result
IsEmpty(cs_zero, R);
    true
Display(g, R);
    [[1, a+1 <> 0], [x-a, a+1 = 0]]
```
