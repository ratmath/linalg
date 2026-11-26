# @ratmath/linalg

Exact rational linear algebra for JavaScript, built to feel instant.

Every operation answers twice: first a Float64 **preview** that is ready in
microseconds (what the screen shows while a slider is moving), then the
**exact** rational answer, computed in small slices that never block the page
and are dropped the moment the question changes. The exact answer is always
the one that stays.

```js
import { solve, det, live } from "@ratmath/linalg";

const job = solve([[2, 1, 1], [1, 3, 2], [1, 0, 0]], [4, 5, 6]);
job.quick;          // { stage: "float", value: { status: "unique", x: [[6], [15.000000000000004], [-23.000000000000004]] }, trust: "good" }
job.exact().x;      // Matrix [[6], [15], [-23]]   (exact rationals)
await job.settle(); // same, but in 8 ms slices that keep the page responsive

det([["1/2", "1/3"], ["1/4", "1/5"]]).exact().toString(); // "1/60"
```

Docs and live demos: https://ratmath.com/linalg/ · Source:
https://github.com/ratmath/linalg · Siblings:
[core](https://ratmath.com/core/), [algebraic](https://ratmath.com/algebraic/),
[analytic](https://ratmath.com/analytic/), [graphics](https://ratmath.com/graphics/)

```sh
npm install @ratmath/linalg @ratmath/core
```

ES modules for Bun, Node 22+ and browsers, with TypeScript declarations.
While the version is 0.x, a minor version may change an API.

Part of [RatMath](https://ratmath.com). Entries are exact rationals from
[`@ratmath/core`](https://github.com/ratmath/core): numbers, bigints, strings
in Core's literal grammar (`"3/4"`, `"1.25"`, `"0.#3"`), or Core `Rational`s.
A decimal JS number is read by its shortest spelling, so `0.1` is `1/10`.

## Fast first, exact last

Each operation returns a `Refinement`:

| member | what it gives |
| --- | --- |
| `quick` | float preview `{ stage: "float", value, trust }`, computed synchronously on first read |
| `exact()` | the exact answer, synchronously |
| `exactStage()` | `{ stage: "exact", value, via }`; `via` names the method that settled it |
| `settle({ signal, sliceMs, onStage, onProgress })` | promise of the exact answer, run in slices; rejects with `AbortError` when `signal` aborts |
| `stages()` | generator of every stage: float, progress ticks, exact |
| `settled` | true once the exact answer is cached |

`trust` on the preview says how far to believe it:

- `"good"`: well conditioned; the digits shown are very likely right.
- `"rough"`: ill conditioned; the shape is right, trailing digits are not.
- `"unsure"`: a pivot or rank decision fell inside rounding noise, so the exact
  pass may change the kind of answer (singular or not, the rank).

A display can use it directly: show `good` previews normally, show `rough` ones
with fewer digits, and mark `unsure` ones as provisional until `exact` lands.

### Sliders: `live`

```js
const update = live(
  (t) => solve([[1, t], [t, 2]], [1, 0]),
  (stage) => draw(stage.value, { exact: stage.exact }),
  { idleMs: 30 },
);
slider.oninput = () => update(slider.value);
```

Each update draws the preview at once. Exact work starts once the input has
been still for `idleMs`, runs in slices, and is cancelled by the next update,
so only the latest question is ever finished. `update.cancel()` stops pending
work.

The `Refinement` shape is meant to be shared across RatMath: a graphics or
analytic computation that exposes `quick`, `stages()`, `exact()` and
`settle()` (build one with `new Refinement(quickFn, exactGeneratorFn)`) can be
driven by the same `live` loop.

## Operations

| function | exact value |
| --- | --- |
| `det(A)` | Core `Rational` |
| `rank(A)` | integer |
| `rref(A, { columns })` | `{ matrix: Matrix, pivots: number[] }` |
| `solve(A, b)` | `{ status: "unique" \| "infinite" \| "none", x: Matrix \| null, basis: Matrix[] }` |
| `inverse(A)` | `Matrix`, or `null` when singular |
| `nullspace(A)` | `Matrix[]` (column vectors) |
| `multiply(A, B)` | `Matrix` (immediate, no refinement needed) |

`solve` accepts several right-hand sides as columns of `b`. When there are
infinitely many solutions, `x` is one of them (free unknowns set to 0) and
every `x + Σ cᵢ·basisᵢ` is also a solution. Float previews have the same
shapes with `number[][]` in place of `Matrix`.

`Matrix` is immutable: `Matrix.from(rows)`, `Matrix.identity(n)`,
`Matrix.zeros(r, c)`, `get(i, j)` (a Core `Rational`), `toArray()`,
`toStrings()`, `toNumbers()`, `add`, `subtract`, `scale`, `multiply`,
`transpose`, `augment`, `column`, `equals`. A flat array is read as a column.

## How the exact answers are found

- **Fraction-free elimination (Bareiss).** Each row is cleared of
  denominators, then eliminated with integer arithmetic only. Every division is
  exact and every intermediate entry is a minor of the input, so numbers grow
  as little as they can. Used for `rref`, `rank`, general `solve`, `inverse`.
- **Multi-modular determinant.** The determinant is computed modulo word-size
  primes (below 2²⁶, so all the O(n³) work is on plain doubles) and rebuilt
  with the Chinese remainder theorem once the primes pass twice the Hadamard
  bound. Deterministic, no guessing.
- **Guess and check.** When the preview of a square system or inverse is
  trustworthy, its floats are turned into the nearest small fractions and
  checked exactly: `A·x = b` plus a modular proof that `A` is nonsingular, or
  `A·X = I`. The check is O(n²) (O(n³) for the inverse) with small numbers,
  so answers with modest denominators, the common case in a calculator, become
  exact almost at once. A failed check falls back to elimination; a guess is
  never returned unchecked.

Rough timings on random matrices with entries like `-7/3` (Bun, one core):

| n | float preview | exact det | exact solve | exact solve, nice answer |
| --- | --- | --- | --- | --- |
| 10 | ≈1 ms | 3 ms | 5 ms | 0.5 ms |
| 50 | ≈2 ms | 21 ms | 35 ms | 5 ms |
| 100 | ≈3 ms | 91 ms | 318 ms | 6 ms |

The low-level generators are exported as `exact` (`bareiss`, `rrefPairs`,
`solvePairs`, `inversePairs`, `determinantModular`, `detMod`, `run`, ...) and
the float routines as `float`.

## Inexact entries: `field`

For entries that are ranges or enclosures rather than rationals,
`field.rref`, `field.det`, `field.inverse` and `field.multiply` take a field
object (`zero, one, add, sub, mul, div, neg, sign`). `field.rationalField`
and `field.intervalField` (Core `RationalInterval`s) are provided. When a pivot
cannot be told apart from 0, the `undecided(row, col)` option is called, so the
caller can ask for more digits and run again.

## Beyond solving

All of these follow the same pattern: `quick` is a float preview, `exact()` /
`settle()` the exact or certified answer.

### Eigenvalues, eigenvectors, singular values

| function | preview | exact |
| --- | --- | --- |
| `charpoly(A)` | float coefficients | `Rational[]` of det(xI − A), constant term first |
| `eigenvalues(A)` | `[{ re, im }]` (QR algorithm) | `{ real: [{ value, multiplicity }], complex, complexCount, charpoly }` |
| `eigenvectors(A)` | `[{ value, vector }]` (inverse iteration) | `[{ value, multiplicity, exact, vectors, approximate(bits) }]` |
| `singularValues(A)` | `number[]` (Jacobi) | `[{ value, multiplicity }]`, descending, zeros included |
| `svd(A)` | `{ U, s, V }` | the same plus certified `values` |

Real eigenvalues and singular values are `AlgebraicReal`s: a square-free
integer polynomial and an interval holding exactly one of its roots.
`value.refine(bits)` shrinks the interval as far as you ask,
`value.interval(bits)` returns it as two Core `Rational`s. When the number is
rational (decided exactly, not guessed), `value.toRational()` is the Core
`Rational` and `value.rational` is the same value as a `[numerator,
denominator]` pair of bigints; both are `null` otherwise. `toNumber()` and
`toString(digits)` print a real number; a complex one has `toNumbers()`,
which returns `{ re, im }`. This is
the fast-then-exact idea applied to irrational answers: the number is exact
(its polynomial and isolating interval), and you pick the digits.

How: the characteristic polynomial is found exactly (Hessenberg reduction
modulo word-size primes, joined by the Chinese remainder theorem), split into
square-free factors, and each factor's real roots are isolated with Sturm
sequences. Complex eigenvalues are `ComplexAlgebraic`s: Durand–Kerner finds
them in floating point, then a Krawczyk test on the real and imaginary parts of
p(x + iy) proves each box holds exactly one root (boxes are checked disjoint
and conjugate pairs matched). `refine(bits)` shrinks a box. If certification
fails the entry is `{ value: null, approx, multiplicity }`.

Eigenvectors are exact for every eigenvalue. For a rational λ, `vectors` is a
basis of the eigenspace. For an irrational or complex λ, `basis` is a list of
`AlgebraicVector`s over Q(λ): each entry is a rational polynomial in λ
(`toStrings()` prints them, `approximate(bits)` and `toNumbers()` evaluate
them), and `algebraic` is the first of them. A simple eigenvalue's vector is
read off a column of the adjugate of λI − A (Faddeev–LeVerrier). A repeated
one is found by row reducing A − λI over the number field Q(λ), which gives
the whole eigenspace whatever its dimension.

### Factorizations and constructions

| function | exact value |
| --- | --- |
| `lu(A)` | `{ perm, L, U }` with PA = LU, any shape or rank |
| `ldl(A)` | `{ L, D, positiveDefinite }` with A = L·D·Lᵀ (symmetric A; `null` if it needs pivoting) |
| `qr(A)` | `{ Q, R, normsSquared }`, square-root-free: Q has orthogonal columns, R is unit upper triangular |
| `gramSchmidt(A)` | orthogonal basis of the column space (exact, unnormalised) |
| `power(A, k)` | Aᵏ for any integer k (negative uses the inverse) |
| `pseudoinverse(A)` | Moore–Penrose A⁺ (rational for rational A) |
| `leastSquares(A, b)` | `{ x, residual }`, x of least norm minimising ‖Ax − b‖ |
| `trace(A)`, `spaces(A)` | immediate: trace; rank, nullity, column-space and row-space bases |
| `similarity(A)` / `jordan(A)` | `{ invariantFactors, minimalPolynomial, rationalCanonical, jordan, J, P, algebraicForm() }` |
| `expm(A, { bits })` | a matrix (array of rows) of Core `RationalInterval`s, each certainly containing the entry of exp(A) and narrower than 2^-bits |
| `sinm`, `cosm`, `sinhm`, `coshm(A, { bits })` | the same for sin, cos, sinh and cosh of a matrix |
| `sqrtm(A, { bits })`, `logm(A, { bits })` | the principal square root and logarithm as certified intervals; single points when the root is rational; a `RangeError` when A has an eigenvalue on the negative real axis or at 0 |
| `powm(A, t, { bits })` | A to a rational power t, as exp(t·log A) |
| `norm(A, p)` | an exact algebraic number: p = 1, `"inf"`, `"fro"` or 2 |
| `cond(A, p)` | the condition number ‖A‖·‖A⁻¹‖, exact; `Infinity` when singular |
| `definiteness(A)` | `{ kind, positive, negative, zero }` for symmetric A |
| `cholesky(A)` | `{ L, D }` with A = L·Lᵀ and exact square roots, or `null` |
| `symmetricEigen(A)` | `{ values, Q }` with A = Q·diag(values)·Qᵀ and QᵀQ = I, exact |

These return at once (no Refinement), on integer matrices:

| function | value |
| --- | --- |
| `hermite(A)` | `{ H, U, pivots }`, U·A = H |
| `smith(A)` | `{ S, U, V, diagonal }`, U·A·V = S with d₁ ∣ d₂ ∣ … |
| `integerKernel(A)` | a basis of the integer solutions of A·x = 0, as columns, LLL-reduced |
| `integerImage(A)` | a basis of the lattice spanned by the columns of A, LLL-reduced |
| `lllReduce(A)` | the rows of A replaced by a short, nearly orthogonal basis of the same lattice |

`similarity` computes the Smith form of xI − A over Q[x]. Its invariant factors
give the minimal polynomial and the rational canonical form (always rational);
their square-free parts give the Jordan block sizes of every eigenvalue,
irrational and complex ones included. When all eigenvalues are rational, `J`
and `P` (from Jordan chains, A = PJP⁻¹ checked exactly) are filled in as
Matrices. For any eigenvalues, `algebraicForm()` returns `{ J, P }` with exact
algebraic entries (√2, i, …), built from Jordan chains over each Q(λ); it is
computed when called. `expm` uses scaling and squaring on interval Taylor
series with a rigorous tail bound and outward rounding, so the true exp(A)
lies in the returned intervals entrywise.

The other matrix functions follow the same pattern. `sinm` and `cosm` sum the
series of A/2ˢ and undo the scaling with the double-angle rules. `sqrtm` takes
a float guess, proves with a Krawczyk test on X² − A = 0 that a small box holds
exactly one root, and proves with a Lyapunov test (a positive definite P with
XP + PXᵀ positive definite) that the root's eigenvalues have positive real
part, which makes it the principal root. `logm` takes proven square roots
until the matrix is near I, sums the series of log(I + X) with a tail bound,
and doubles back. These need nothing from `@ratmath/analytic`: a matrix
function only needs the series and a norm bound, not Taylor models.

```js
import { sqrtm, logm, sinm } from "@ratmath/linalg";

sqrtm([[2, 1], [1, 3]]).quick.value;                   // [[1.3763…, 0.3249…], [0.3249…, 1.7013…]]
sqrtm([[1, 1], [0, 1]]).exact()[0][1].low.toString();   // "1/2", a rational root is exact
logm([[0, -1], [1, 0]], { bits: 100 }).exact()[1][0];   // an interval of width < 2^-100 around π/2
sinm([[0, 1], [0, 0]]).exact();                         // encloses [[0, 1], [0, 0]]
```

```js
import { jordan, norm, cond, cholesky, symmetricEigen, algebraic } from "@ratmath/linalg";

const { J, P } = jordan([[0, 2], [1, 0]]).exact().algebraicForm();
J.map((row) => row.map(algebraic.format));           // [["−√2", "0"], ["0", "√2"]]
algebraic.format(norm([[1, 2], [3, 4]], "fro").exact());   // "√30"
algebraic.format(cond([[1, 2], [3, 4]], 2).exact());       // "(15 + √221)/2"
cholesky([[4, 2], [2, 3]]).exact().L;                // [[2, 0], [1, √2]] as algebraic numbers
symmetricEigen([[2, 1], [1, 2]]).exact();            // values 1, 3; Q with entries ±√2/2
```

The norms: 1 is the largest column sum and `"inf"` the largest row sum (both
rational), `"fro"` is the square root of the sum of squares, and 2 is the
largest singular value. `definiteness` counts the positive eigenvalues of a
symmetric matrix with Descartes' rule of signs on the characteristic
polynomial, which is exact because all its roots are real. An "algebraic
matrix" is an array of rows of exact algebraic numbers; `field.multiply`,
`field.det`, `field.inverse` and `field.rref` work on them with
`complexField`, and `algebraicMatrix(rows)` builds one from rationals and
algebraic numbers.

The QR preview is the usual Householder form with orthonormal Q. The exact
form avoids square roots; the orthonormal Q is Q·diag(1/√normsSquared).

### Nonlinear support

Nonlinear solvers keep dipping into linear algebra: Newton steps solve with the
Jacobian, and certifying a root needs interval linear algebra. Two pieces are
here:

- `certifyRoot({ f, jacobian, guess, bits })` finds a zero of f : Rⁿ → Rⁿ.
  `f` receives Core `Rational`s and returns exact values (or enclosures);
  `jacobian` receives Core `RationalInterval`s (a box) and returns a matrix
  enclosing the Jacobian over it. The preview is float Newton from `guess`.
  The exact stage grows a box around the Newton point until the **Krawczyk
  test** proves exactly one zero inside, then contracts it to width 2^-bits:
  `{ status: "unique" | "none" | "unknown", box: RationalInterval[] }`.
- `solveIntervals(A, b)` encloses every solution of A′x = b′ for A′, b′ ranging
  over interval entries (Core `RationalInterval`s or `[lo, hi]`), which is what
  the calculator's inexact coefficients need. `status: "unique"` also proves
  every A′ in the range nonsingular.

`interval` exports the pieces (interval arithmetic on exact endpoints, outward
rounding, `krawczykStep`, `krawczyk`) for solvers that want to drive their own
loop, for example bisecting boxes the test leaves `"unknown"`.

### Polynomials

`poly` exports the exact polynomial tools used above: `squareFree`,
`isolate`, `countRealRoots`, `polyGcd`, `charpolyModular`,
`simplestBetween`, and `AlgebraicReal`.

## Solving at scale

`solve` tries, in order: guess the answer from the float preview with continued
fractions and check it exactly (a modular proof shows A is nonsingular); for
n ≥ 60, Dixon p-adic lifting with rational reconstruction; otherwise
fraction-free Bareiss. `via` on the exact stage says which won. A 100×100
integer system takes about 130 ms with Dixon against 330 ms for Bareiss.

Rough timings for exact eigenvalues of a random 20×20 integer matrix: 50 ms
for the real ones, about 160 ms with every complex one certified; singular
values about 95 ms.

## Algebraic numbers and other fields

`algebraic` is [@ratmath/algebraic](https://github.com/ratmath/algebraic):
exact real and complex algebraic numbers (`add`, `sub`, `mul`, `div`, `inv`,
`pow`, `sqrt`, `nthRoot`, `compare`, `roots`), polynomial factoring over the
integers, number fields such as Q(√2) and Q(i), and finite fields. Results are
put on their minimal polynomials, so degrees stay small. Eigenvalues and the
polynomial and interval machinery underneath them come from the same package.

`field` runs elimination over any field object `{ zero, one, add, sub, mul,
div, neg, sign }`. Provided: `rationalField`, `intervalField` (zero tests may be
undecided), `primeField(p)` (integers mod a prime, with `from(rational)`) and
`algebraicField` (exact zero tests). `NumberField` and `galoisField(p, k)` from
@ratmath/algebraic are field objects too.

## Workers

`@ratmath/linalg/worker` exports `serveLinalg(self)` for the worker side and
`workerClient(worker)` for the page. The client returns Refinements whose
preview is computed locally (instant) and whose exact stage runs in the worker,
so a long exact computation never blocks the page. Values cross the boundary
as plain data (`src/serialize.js`).

## Teaching forms

`teach` (also `@ratmath/linalg/teach`) has slower, explained versions of the
algorithms. Each returns `{ title, summary, steps, result }`, and each step
carries text plus the matrices and highlights to draw: `gaussJordan`,
`solveSystem`, `inverse`, `determinantByElimination`, `cofactorExpansion`,
`bareiss`, `modularDeterminant`, `guessAndCheck`, `gramSchmidt`, `lu`,
`leastSquares`, `eigen`, `sturmWalk`, `newtonKrawczyk`, `powerBySquaring`,
`jordanBlocks`, `exponentialSeries`, `fourSubspaces`, `dixonLifting`. `RowReducer` is an undoable row-operation
sandbox with `hint()` for the next sensible move. The fast functions above
never pay for this bookkeeping.

## Numerical methods for study: `methods`

The classical methods, run in doubles as a textbook runs them, with what
exact arithmetic says beside every step: `jacobi`, `gaussSeidel`, `sor`
(with the iteration matrix and its spectral radius), `conjugateGradient`
(exact in at most n steps without rounding), `powerMethod`,
`inverseIteration`, `qrAlgorithm` and `pivoting`. Each row gives the double
iterate's `error`, the `methodError` of the same iteration in exact rationals,
and the `rounding` between the two, all computed exactly. See the docs page
"Iterative and eigenvalue methods".

```js
import { methods } from "@ratmath/linalg";
methods.conjugateGradient(hilbert4, [1, 1, 1, 1]).exactSteps;   // 4
```

## Docs site

`docs/` is the static site published at https://ratmath.com/linalg/: an
introduction, a reference, topic pages with live demos, and a free-form exact
workbench. How it is built is in [CONTRIBUTING.md](CONTRIBUTING.md).

## TypeScript

`index.d.ts` declares the public API.

## Not yet

Faster root isolation (Descartes / continued fractions) for very large
characteristic polynomials. Matrix functions of a matrix with eigenvalues on
the negative real axis (non-principal roots and logarithms), symbolic
parameters and sparse storage are not planned yet.

## Develop

```sh
git submodule update --init   # @ratmath/algebraic, until it is on npm
bun install
bun test
```

MIT licensed.
