# Changelog

## Unreleased

- Requires Core 0.12: the `@ratmath/core` peer range is now `>=0.12.0` (Core 0.12 returns Integers from `floor`/`ceil`/`trunc`/`round` and rounds half up by default).
- Tests: the rational reconstruction test no longer computes an unused
  5^(m − 2) (about 10⁸ bits, past some engines' BigInt cap) and now checks
  several residues, an unrecoverable one and the bounds; the large Dixon test
  compares A·x with the given right-hand side instead of with itself.
- `@ratmath/core` is now a peer dependency (`>=0.10.0 <0.12.0`): the package
  uses the program's own copy of Core instead of bringing its own, so Core 0.11
  works and there is never a second copy whose `Rational` is not the same class.

### Fixed

- The `det` preview no longer shows 0 for a matrix whose determinant is tiny
  but not zero (the 12×12 Hilbert matrix now previews about 2.8e-78).
- The Krawczyk step keeps its midpoint inside the box (through
  `@ratmath/algebraic`).
- README: install instructions, links, the shape of `expm`'s value, and the
  names `value.rational` and `toNumbers()`.

### Added

- `methods`: Jacobi, Gauss–Seidel, SOR, conjugate gradients, the power
  method, inverse iteration, the QR algorithm and pivoting, run in doubles
  with the exact solution, the exact-rational iterate and the rounding error
  beside every step (roadmap 13.3). Docs page "Iterative and eigenvalue
  methods".
- Matrix functions with certified interval entries: `sinm`, `cosm`, `sinhm`,
  `coshm`, the principal square root `sqrtm`, the principal logarithm `logm`
  and rational powers `powm(A, t)`. A rational square root is returned
  exactly. Docs page "Matrix functions" and workbench functions of the same
  names.
- `eigenvectors`: every entry has `basis`, a full exact basis of the
  eigenspace. Repeated irrational eigenvalues and complex eigenvalues are
  solved over the number field Q(λ) instead of approximated.
- `jordan` and `similarity`: `algebraicForm()` returns J and P as exact
  algebraic numbers when the eigenvalues are not rational.
- `norm(A, p)` and `cond(A, p)` for p = 1, `"inf"`, `"fro"` and 2, as exact
  algebraic numbers.
- `definiteness(A)`, `cholesky(A)` with exact square roots, and
  `symmetricEigen(A)` with an exact orthonormal eigenbasis.
- `integerKernel(A)`, `integerImage(A)` and `lllReduce(A)`.
- `complexField`, `algebraicMatrix(rows)` and `numberFieldOps(K)`: matrices of
  real and complex algebraic numbers through the `field` functions.
- Docs: a code snippet on every page, the new functions in the reference,
  page descriptions and links to the other libraries.
- `CONTRIBUTING.md`, this changelog and a CI workflow.

### Changed

- `Refinement` and `live` now come from `@ratmath/core/refine` (Core 0.10.0),
  one class shared by every RatMath library. The exports here are unchanged.
- Needs `@ratmath/core` 0.10.0 or later.
- An `eigenvectors` entry with an irrational eigenvalue of multiplicity above
  one used to carry approximate vectors with `exact: false`; it now carries
  exact ones.

## 0.1.0

First version.
