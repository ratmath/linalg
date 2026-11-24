# Changelog

## Unreleased

### Fixed

- The `det` preview no longer shows 0 for a matrix whose determinant is tiny
  but not zero (the 12×12 Hilbert matrix now previews about 2.8e-78).
- The Krawczyk step keeps its midpoint inside the box (through
  `@ratmath/algebraic`).
- README: install instructions, links, the shape of `expm`'s value, and the
  names `value.rational` and `toNumbers()`.

### Added

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

- Needs `@ratmath/core` 0.9.0 or later.
- An `eigenvectors` entry with an irrational eigenvalue of multiplicity above
  one used to carry approximate vectors with `exact: false`; it now carries
  exact ones.

## 0.1.0

First version.
