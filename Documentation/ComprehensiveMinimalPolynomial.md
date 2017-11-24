# ComprehensiveMinimalPolynomial
The `ComprehensiveMinimalPolynomial` function will compute a complete case discussion of the minimal polynomial of a matrix of multivariate polynomials.

## Calling Sequence
```
ComprehensiveMinimalPolynomial(A, v, R, options)
ComprehensiveMinimalPolynomial(A, v, rs, R, options)
ComprehensiveMinimalPolynomial(A, v, cs, R, options)
ComprehensiveMinimalPolynomial(A, v, F, R, options)
ComprehensiveMinimalPolynomial(A, v, F, H, R, options)
```

## Input

| Variable | Description |
| --- | --- |
| `A`  | Square matrix of multivariate polynomials over `R` |
| `v` | Variable to use for the computed minimal polynmoials |
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
| `[p, rs]` | `'outputType'` is `'RegularSystem'` or `'RS'` |
| `[p, cs]` | `'outputType'` is `'ConstructibleSet'` or `'CS'` |

Where `p` is the minimal polynomial of the input matrix for all parameter values in the zero set of `cs` or `rs`. Together, all the constructible sets or regular systems form a partition of the input constructible set, regular system or set of equations and inequations.
