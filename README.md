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

## Not yet

Eigenvalues and characteristic polynomials, LU/QR/Cholesky as exact
factorizations, p-adic (Dixon) lifting for faster large exact solves, a Web
Worker wrapper, TypeScript declarations.

## Develop

```sh
bun install
bun test
```

MIT licensed.
