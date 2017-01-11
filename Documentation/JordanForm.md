# JordanForm
The `JordanForm` method is an improved version of the `JordanForm` method from the `LinearAlgebra` package for computing the [Jordan normal form](https://en.wikipedia.org/wiki/Jordan_normal_form). It returns the correct result when the matrix has parameters.

## Calling Sequence
```
JordanForm(A)
```

## Input

| Variable | Description |
| --- | --- |
| `A`  | Square matrix |

## Output

The Jordan form of the input matrix.

## Example
```
with(LinearAlgebra):
with(ParametricMatrixTools):

p := x^5 + x^4 + x^3 + x^2 + x + a:
F := CompanionMatrix(p, x):

# The LinearAlgebra implementation fails with a parameter
LinearAlgebra:-JordanForm(F);
                             ([0 0 0 0 -a])
                             ([1 0 0 0 -1])
    LinearAlgebra:-JordanForm([0 1 0 0 -1])
                             ([0 0 1 0 -1])
                             ([0 0 0 1 -1])

# The ParametricMatrixTools implementation gives the correct result
ParametricMatrixTools:-JordanForm(F);
    [RootOf(_Z^5+_Z^4+_Z^3+_Z^2+_Z+a, index = 1), 0, 0, 0, 0]
    [0, RootOf(_Z^5+_Z^4+_Z^3+_Z^2+_Z+a, index = 2), 0, 0, 0]
    [0, 0, RootOf(_Z^5+_Z^4+_Z^3+_Z^2+_Z+a, index = 3), 0, 0]
    [0, 0, 0, RootOf(_Z^5+_Z^4+_Z^3+_Z^2+_Z+a, index = 4), 0]
    [0, 0, 0, 0, RootOf(_Z^5+_Z^4+_Z^3+_Z^2+_Z+a, index = 5)]
```
