// The wider toolkit: spectral data, factorizations, powers, least squares,
// integer normal forms, and certified nonlinear root finding. Same pattern as
// ops.js: a float preview first, the exact (or certified) answer after.

import { Matrix } from "./matrix.js";
import { Refinement } from "./refine.js";
import { toRational, toPair, reduce, pairToNumber, addQ } from "./q.js";
import { run, rrefPairs, inversePairs, solvePairs } from "./exact.js";
import { solveFloat, inverseFloat } from "./float.js";
import { eigenvaluesFloat, eigenvectorFloat, svdFloat, qrFloat } from "./eigen-float.js";
import { eigenExact, singularExact } from "./spectral.js";
import {
  luPairs, ldlPairs, gramSchmidtPairs, powerPairs, multiplyPairs, transposePairs, hermiteInts, smithInts,
} from "./decomp.js";
import { iv, toCore, krawczyk, ivWidth } from "./interval.js";
import { structurePairs, adjugatePolys, algebraicEigenvector } from "./structure.js";
import { expmPairs, expmFloat } from "./expm.js";

const M = (x) => Matrix.from(x);
const square = (A, what) => {
  if (!A.isSquare) throw new RangeError(`${what} needs a square matrix, got ${A.rows}×${A.cols}`);
};
const rowsOf = (A) => A.toNumbers();

// ---------------------------------------------------------------------------
// Spectral

/** Coefficients of det(xI − A), constant term first (monic). */
export function charpoly(input) {
  const A = M(input);
  square(A, "charpoly");
  return new Refinement(
    () => {
      const ev = eigenvaluesFloat(A.float(), A.rows);
      if (!ev) return { value: null, trust: "unsure" };
      // expand Π (x − λ) with complex arithmetic
      let re = [1];
      let im = [0];
      for (const { re: a, im: b } of ev) {
        const nr = new Array(re.length + 1).fill(0);
        const ni = new Array(re.length + 1).fill(0);
        for (let k = 0; k < re.length; k += 1) {
          nr[k + 1] += re[k];
          ni[k + 1] += im[k];
          nr[k] -= a * re[k] - b * im[k];
          ni[k] -= a * im[k] + b * re[k];
        }
        re = nr;
        im = ni;
      }
      return { value: re, trust: "rough" };
    },
    function* () {
      const { charpoly: cp } = yield* eigenExact(A.pairs());
      return { value: cp.map(toRational), via: "hessenberg" };
    },
  );
}

/**
 * Eigenvalues. Preview: [{ re, im }] from the QR algorithm.
 * Exact: { real: [{ value: AlgebraicReal, multiplicity }],
 *          complex: [{ value: ComplexAlgebraic, multiplicity }], complexCertified,
 *          complexCount, charpoly }
 * Real eigenvalues are certified: each is a root of a known square-free integer
 * polynomial in an interval that `value.refine(bits)` shrinks as far as asked;
 * `value.rational` is the exact Rational when there is one. Complex
 * eigenvalues are ComplexAlgebraic boxes proven (Krawczyk) to hold exactly one
 * root each; if that ever fails, complexCertified is false and they are
 * given as { value: null, approx: { re, im } }.
 */
export function eigenvalues(input) {
  const A = M(input);
  square(A, "eigenvalues");
  const job = new Refinement(
    () => {
      const ev = eigenvaluesFloat(A.float(), A.rows);
      return { value: ev, trust: ev ? "rough" : "unsure" };
    },
    function* () {
      const e = yield* eigenExact(A.pairs());
      const approx = job.quick.value ?? [];
      return {
        value: {
          real: e.real,
          // certified boxes when possible, otherwise float approximations
          complex: e.complexCertified ? e.complex : approx.filter((z) => z.im !== 0).map((z) => ({ value: null, approx: z, multiplicity: 1 })),
          complexCertified: e.complexCertified,
          complexCount: e.complexCount,
          charpoly: e.charpoly.map(toRational),
        },
        via: "sturm",
      };
    },
  );
  return job;
}

/**
 * Eigenvectors for the real eigenvalues.
 * Preview: [{ value: number, vector: number[] }].
 * Exact: [{ value, multiplicity, exact, vectors, algebraic, approximate(bits) }]
 *   rational eigenvalue: vectors = exact basis of the eigenspace
 *   irrational eigenvalue: algebraic = an AlgebraicVector whose entries are
 *     polynomials in λ (an exact element of Q(λ)ⁿ, from the adjugate of
 *     λI − A); approximate(bits) evaluates it at a rational point near λ
 *   complex eigenvalue (when certified): algebraic as above, vectors = []
 * When the adjugate vanishes at λ (an irrational eigenvalue with several
 * independent eigenvectors) exact is false and approximate() falls back to
 * one exact inverse-iteration step.
 */
export function eigenvectors(input) {
  const A = M(input);
  square(A, "eigenvectors");
  const n = A.rows;
  const job = new Refinement(
    () => {
      const ev = eigenvaluesFloat(A.float(), n);
      if (!ev) return { value: null, trust: "unsure" };
      const real = ev.filter((z) => z.im === 0);
      return { value: real.map(({ re }) => ({ value: re, vector: eigenvectorFloat(A.float(), n, re) })), trust: "rough" };
    },
    function* () {
      const e = yield* eigenExact(A.pairs());
      const out = [];
      let adj = null;
      for (const { value, multiplicity } of e.real) {
        const lam = value.rational;
        if (lam) {
          const shifted = A.pairs().map((row, i) => row.map((x, j) => (i === j ? addQ(x, [-lam[0], lam[1]]) : x)));
          const { basis } = yield* solvePairs(shifted, Array.from({ length: n }, () => [[0n, 1n]]));
          const vectors = basis.map((v) => Matrix.fromPairs(v.map((p) => [p])));
          out.push({ value, multiplicity, exact: true, vectors, approximate: () => vectors[0] });
        } else {
          adj ??= adjugatePolys(A.pairs());
          const algebraic = algebraicEigenvector(adj, value);
          const approximate = (bits = 53) => (algebraic ? algebraic.approximate(bits) : approxVector(A, value, bits));
          out.push({ value, multiplicity, exact: Boolean(algebraic), algebraic, vectors: [approximate(53)], approximate });
        }
      }
      if (e.complexCertified) {
        for (const { value, multiplicity } of e.complex) {
          adj ??= adjugatePolys(A.pairs());
          const algebraic = algebraicEigenvector(adj, value);
          out.push({ value, multiplicity, exact: Boolean(algebraic), algebraic, vectors: [], approximate: null });
        }
      }
      return { value: out, via: "nullspace" };
    },
  );
  return job;
}

function approxVector(A, value, bits) {
  const n = A.rows;
  value.refine(bits);
  const lt = reduce(value.lo[0] * value.hi[1] + value.hi[0] * value.lo[1], 2n * value.lo[1] * value.hi[1]);
  const start = eigenvectorFloat(A.float(), n, pairToNumber(...lt)).map((x) => [toPair(x)]);
  for (let nudge = 0n; nudge < 5n; nudge += 1n) {
    const shift = nudge ? reduce(lt[0] * 2n ** BigInt(bits + 8) + nudge * lt[1], lt[1] * 2n ** BigInt(bits + 8)) : lt;
    const shifted = A.pairs().map((row, i) => row.map((x, j) => (i === j ? addQ(x, [-shift[0], shift[1]]) : x)));
    const { status, x } = run(solvePairs(shifted, start));
    if (status !== "unique") continue;
    // normalise so the largest entry is 1
    const nums = x.map(([p]) => Math.abs(pairToNumber(...p)));
    const k = nums.indexOf(Math.max(...nums));
    const [a, b] = x[k][0];
    return Matrix.fromPairs(x.map(([[p, q]]) => [reduce(p * b, q * a)]));
  }
  return null;
}

/**
 * Singular values, descending. Preview: numbers. Exact: [{ value: AlgebraicReal,
 * multiplicity }], each σ certified as a root of f(x²) with f a factor of the
 * characteristic polynomial of AᵀA (or AAᵀ, whichever is smaller).
 */
export function singularValues(input) {
  const A = M(input);
  return new Refinement(
    () => ({ value: svdFloat(A.float(), A.rows, A.cols).s, trust: "good" }),
    function* () {
      return { value: yield* singularExact(A.pairs()), via: "sturm" };
    },
  );
}

/**
 * Thin SVD A = U·diag(s)·Vᵀ. Preview: { U, s, V } in floats. Exact stage keeps
 * U and V as floats (they are irrational in general) and adds certified
 * singular values: { U, s, V, values }.
 */
export function svd(input) {
  const A = M(input);
  const job = new Refinement(
    () => ({ value: svdFloat(A.float(), A.rows, A.cols), trust: "good" }),
    function* () {
      const values = yield* singularExact(A.pairs());
      return { value: { ...job.quick.value, values }, via: "sturm" };
    },
  );
  return job;
}

// ---------------------------------------------------------------------------
// Factorizations

function luFloat(A) {
  const m = A.rows;
  const n = A.cols;
  const U = rowsOf(A);
  const L = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? 1 : 0)));
  const perm = Array.from({ length: m }, (_, i) => i);
  let r = 0;
  for (let col = 0; col < n && r < m; col += 1) {
    let p = r;
    for (let i = r + 1; i < m; i += 1) if (Math.abs(U[i][col]) > Math.abs(U[p][col])) p = i;
    if (U[p][col] === 0) continue;
    [U[p], U[r]] = [U[r], U[p]];
    [perm[p], perm[r]] = [perm[r], perm[p]];
    for (let j = 0; j < r; j += 1) [L[p][j], L[r][j]] = [L[r][j], L[p][j]];
    for (let i = r + 1; i < m; i += 1) {
      const f = U[i][col] / U[r][col];
      L[i][r] = f;
      for (let j = col; j < n; j += 1) U[i][j] -= f * U[r][j];
    }
    r += 1;
  }
  return { perm, L, U };
}

/** PA = LU. Value: { perm, L, U } (row i of PA is row perm[i] of A). The preview pivots for size, the exact form on the first nonzero entry. */
export function lu(input) {
  const A = M(input);
  return new Refinement(
    () => ({ value: luFloat(A), trust: "good" }),
    function* () {
      const { perm, L, U } = yield* luPairs(A.pairs());
      return { value: { perm, L: Matrix.fromPairs(L), U: Matrix.fromPairs(U) }, via: "doolittle" };
    },
  );
}

/** Symmetric A = L·D·Lᵀ (no square roots, so exact). Value { L, D: Rational[] }, or null if pivoting is needed. positiveDefinite tells whether all D > 0. */
export function ldl(input) {
  const A = M(input);
  square(A, "ldl");
  if (!A.equals(A.transpose())) throw new RangeError("ldl needs a symmetric matrix");
  return new Refinement(
    () => {
      const n = A.rows;
      const a = rowsOf(A);
      const L = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
      const D = new Array(n).fill(0);
      for (let j = 0; j < n; j += 1) {
        let d = a[j][j];
        for (let k = 0; k < j; k += 1) d -= L[j][k] ** 2 * D[k];
        D[j] = d;
        for (let i = j + 1; i < n; i += 1) {
          let s = a[i][j];
          for (let k = 0; k < j; k += 1) s -= L[i][k] * L[j][k] * D[k];
          L[i][j] = d ? s / d : 0;
        }
      }
      return { value: { L, D, positiveDefinite: D.every((d) => d > 0) }, trust: "rough" };
    },
    function* () {
      const r = yield* ldlPairs(A.pairs());
      if (!r) return { value: null, via: "ldl" };
      return {
        value: { L: Matrix.fromPairs(r.L), D: r.D.map(toRational), positiveDefinite: r.D.every(([n]) => n > 0n) },
        via: "ldl",
      };
    },
  );
}

/**
 * QR. Preview: Householder, Q with orthonormal columns, R upper with
 * non-negative diagonal. Exact: square-root-free Gram–Schmidt,
 * { Q (orthogonal columns), R (unit upper), normsSquared }, where the
 * orthonormal Q is Q·diag(1/√normsSquared).
 */
export function qr(input) {
  const A = M(input);
  return new Refinement(
    () => ({ value: A.rows >= A.cols ? qrFloat(A.float(), A.rows, A.cols) : null, trust: "good" }),
    function* () {
      const { Q, R, normsSquared } = yield* gramSchmidtPairs(A.pairs());
      return { value: { Q: Matrix.fromPairs(Q), R: Matrix.fromPairs(R), normsSquared: normsSquared.map(toRational) }, via: "gram-schmidt" };
    },
  );
}

/** Orthogonal basis (exact, not normalised) of the span of the columns, zero columns dropped. */
export function gramSchmidt(input) {
  const A = M(input);
  const job = qr(A);
  return new Refinement(
    () => ({ value: job.quick.value?.Q ?? null, trust: job.quick.trust }),
    function* () {
      let last;
      for (const st of job.stages()) {
        if (st.stage === "progress") yield st;
        last = st;
      }
      const { Q, normsSquared } = last.value;
      const keep = normsSquared.map((n, j) => (n.numerator !== 0n ? j : -1)).filter((j) => j >= 0);
      return { value: keep.map((j) => Q.column(j)), via: "gram-schmidt" };
    },
  );
}

/** A^k for any integer k (negative uses the inverse; null if singular). */
export function power(input, k) {
  const A = M(input);
  square(A, "power");
  const e = BigInt(k);
  return new Refinement(
    () => {
      let base = rowsOf(A);
      if (e < 0n) {
        base = inverseFloat(A.float(), A.rows).value;
        if (!base) return { value: null, trust: "unsure" };
      }
      const n = A.rows;
      const mul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((s, x, t) => s + x * Y[t][j], 0)));
      let result = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
      let p = e < 0n ? -e : e;
      while (p > 0n) {
        if (p & 1n) result = mul(result, base);
        p >>= 1n;
        if (p) base = mul(base, base);
      }
      return { value: result, trust: "rough" };
    },
    function* () {
      let base = A.pairs();
      if (e < 0n) {
        base = yield* inversePairs(base);
        if (!base) return { value: null, via: "inverse" };
      }
      const r = yield* powerPairs(base, e < 0n ? -e : e);
      return { value: Matrix.fromPairs(r), via: "squaring" };
    },
  );
}

/** Trace, exact and immediate. */
export function trace(input) {
  const A = M(input);
  square(A, "trace");
  let s = [0n, 1n];
  for (let i = 0; i < A.rows; i += 1) s = addQ(s, A.pair(i, i));
  return toRational(s);
}

// ---------------------------------------------------------------------------
// Least squares and the pseudoinverse

function pinvPairs(A, rref) {
  // full-rank factorisation A = C·F: F = nonzero rows of rref(A), C = pivot columns of A
  const { rows, pivots } = rref;
  const r = pivots.length;
  if (r === 0) return Array.from({ length: A[0]?.length ?? 0 }, () => A.map(() => [0n, 1n]));
  const F = rows.slice(0, r);
  const C = A.map((row) => pivots.map((j) => row[j]));
  const Ft = transposePairs(F);
  const Ct = transposePairs(C);
  const FFt = run(inversePairs(multiplyPairs(F, Ft)));
  const CtC = run(inversePairs(multiplyPairs(Ct, C)));
  return multiplyPairs(multiplyPairs(Ft, FFt), multiplyPairs(CtC, Ct));
}

/** Moore–Penrose pseudoinverse A⁺, exact (it is rational for rational A). */
export function pseudoinverse(input) {
  const A = M(input);
  return new Refinement(
    () => {
      const { U, s, V } = svdFloat(A.float(), A.rows, A.cols);
      const tol = Math.max(A.rows, A.cols) * (s[0] || 0) * 2 ** -52;
      const P = Array.from({ length: A.cols }, (_, i) => Array.from({ length: A.rows }, (_, j) =>
        s.reduce((acc, sk, k) => (sk > tol ? acc + (V[i][k] * U[j][k]) / sk : acc), 0)));
      return { value: P, trust: "good" };
    },
    function* () {
      const pairs = A.pairs();
      const red = yield* rrefPairs(pairs, A.cols);
      return { value: Matrix.fromPairs(pinvPairs(pairs, red)), via: "rank-factorization" };
    },
  );
}

/**
 * Least squares: the x minimising ‖A·x − b‖ (the one of least norm when
 * there are many). Exact: x = A⁺·b; when A has full column rank this is the
 * solution of the normal equations AᵀA·x = Aᵀb.
 * Value: { x, residual } with residual = b − A·x.
 */
export function leastSquares(aInput, bInput) {
  const A = M(aInput);
  let B = M(bInput);
  if (B.rows !== A.rows && B.rows === 1 && B.cols === A.rows) B = B.transpose();
  const job = new Refinement(
    () => {
      const P = pseudoinverse(A).quick.value;
      const b = B.toNumbers();
      const x = P.map((row) => b[0].map((_, c) => row.reduce((s, p, j) => s + p * b[j][c], 0)));
      const a = A.toNumbers();
      const residual = b.map((row, i) => row.map((v, c) => v - a[i].reduce((s, aij, j) => s + aij * x[j][c], 0)));
      return { value: { x, residual }, trust: "good" };
    },
    function* () {
      const pairs = A.pairs();
      const red = yield* rrefPairs(pairs, A.cols);
      let x;
      if (red.pivots.length === A.cols) {
        const At = transposePairs(pairs);
        const r = yield* solvePairs(multiplyPairs(At, pairs), multiplyPairs(At, B.pairs()));
        x = Matrix.fromPairs(r.x);
      } else {
        x = Matrix.fromPairs(multiplyPairs(pinvPairs(pairs, red), B.pairs()));
      }
      return { value: { x, residual: B.subtract(A.multiply(x)) }, via: red.pivots.length === A.cols ? "normal-equations" : "pseudoinverse" };
    },
  );
  return job;
}

/** Bases for the column space (pivot columns of A) and row space (nonzero rows of rref). Exact, immediate. */
export function spaces(input) {
  const A = M(input);
  const { rows, pivots } = run(rrefPairs(A.pairs(), A.cols));
  return {
    rank: pivots.length,
    columnSpace: pivots.map((j) => A.column(j)),
    rowSpace: rows.slice(0, pivots.length).map((r) => Matrix.fromPairs([r])),
    nullity: A.cols - pivots.length,
  };
}

// ---------------------------------------------------------------------------
// Integer normal forms (exact, immediate)

const intMatrix = (rows) => Matrix.fromPairs(rows.map((r) => r.map((v) => [v, 1n])));

/** Hermite normal form: { H, U, pivots } with U·A = H, U unimodular. Integer entries only. */
export function hermite(input) {
  const { H, U, pivots } = hermiteInts(M(input).pairs());
  return { H: intMatrix(H), U: intMatrix(U), pivots };
}

/** Smith normal form: { S, U, V, diagonal } with U·A·V = S, d₁ | d₂ | …. Integer entries only. */
export function smith(input) {
  const { S, U, V, diagonal } = smithInts(M(input).pairs());
  return { S: intMatrix(S), U: intMatrix(U), V: intMatrix(V), diagonal };
}

// ---------------------------------------------------------------------------
// Nonlinear: Newton preview, Krawczyk certificate

/**
 * Find and certify a zero of f : Rⁿ → Rⁿ near `guess`.
 *   f(x)          x: Core Rationals → exact values (Rationals) or enclosures
 *   jacobian(X)   X: Core RationalIntervals → matrix enclosing the Jacobian over X
 * Preview: float Newton iterates from guess. Exact stage: a box around the
 * Newton point is grown until the Krawczyk test certifies a unique zero in
 * it, then contracted to width 2^-bits.
 * Value: { status: "unique" | "none" | "unknown", box: RationalInterval[] | null }.
 */
export function certifyRoot({ f, jacobian, guess, bits = 53, steps = 30 }) {
  const n = guess.length;
  const toR = (x) => toRational(toPair(x));
  const job = new Refinement(
    () => {
      let x = guess.map(Number);
      for (let k = 0; k < steps; k += 1) {
        const fx = f(x.map(toR)).map((v) => pairToNumber(...midPair(iv(v))));
        const J = jacobian(x.map((v) => toCore(iv(v)))).map((row) => row.map((v) => pairToNumber(...midPair(iv(v)))));
        const { value } = solveFloat(Float64Array.from(J.flat()), n, n, fx.map((v) => [v]));
        if (value.status !== "unique") break;
        const dx = value.x.map((r) => r[0]);
        x = x.map((v, i) => v - dx[i]);
        if (Math.max(...dx.map(Math.abs)) <= 1e-15 * Math.max(1, ...x.map(Math.abs))) break;
      }
      return { value: x, trust: "unsure" };
    },
    function* () {
      const center = job.quick.value.map((v) => toPair(v));
      let radius = Math.max(1e-12, 1e-9 * Math.max(...job.quick.value.map(Math.abs)));
      let last = { status: "unknown", box: null };
      for (let attempt = 0; attempt < 12; attempt += 1) {
        const r = toPair(radius);
        const X = center.map((c) => [addQ(c, [-r[0], r[1]]), addQ(c, r)]);
        last = yield* krawczyk(f, jacobian, X, { bits });
        if (last.status === "unique") break;
        if (last.status === "none" && attempt > 4) break;
        radius *= 16;
      }
      return { value: { status: last.status, box: last.box && last.box.map(toCore) }, via: "krawczyk" };
    },
  );
  return job;
}

const midPair = (x) => reduce(x[0][0] * x[1][1] + x[1][0] * x[0][1], 2n * x[0][1] * x[1][1]);

/**
 * Enclose the solutions of A·x = b where entries of A and b are intervals
 * (Core RationalIntervals or [lo, hi]) or exact values. Value:
 * { status: "unique" | "unknown", box } with box ⊇ { x : A'x = b', A' ∈ A, b' ∈ b }
 * when status is "unique" (every A' in A is then nonsingular).
 */
export function solveIntervals(aInput, bInput, { bits = 53 } = {}) {
  const A = aInput.map((row) => row.map(iv));
  const b = bInput.map((x) => iv(Array.isArray(x) && x.length === 1 ? x[0] : x));
  const n = A.length;
  const f = (x) => {
    const xp = x.map(toPair);
    return A.map((row, i) => {
      let lo = [0n, 1n];
      let hi = [0n, 1n];
      row.forEach(([l, h], j) => {
        const c = [mulPair(l, xp[j]), mulPair(h, xp[j])];
        const [cl, ch] = lessPair(c[0], c[1]) ? c : [c[1], c[0]];
        lo = addQ(lo, cl);
        hi = addQ(hi, ch);
      });
      return toCore([addQ(lo, [-b[i][1][0], b[i][1][1]]), addQ(hi, [-b[i][0][0], b[i][0][1]])]);
    });
  };
  const jac = () => A.map((row) => row.map(toCore));
  const mids = A.map((row) => row.map((x) => pairToNumber(...midPair(x))));
  return new Refinement(
    () => {
      const { value, trust } = solveFloat(Float64Array.from(mids.flat()), n, n, b.map((x) => [pairToNumber(...midPair(x))]));
      return { value: value.status === "unique" ? value.x.map((r) => r[0]) : null, trust };
    },
    function* () {
      const guess = solveFloat(Float64Array.from(mids.flat()), n, n, b.map((x) => [pairToNumber(...midPair(x))])).value;
      if (guess.status !== "unique") return { value: { status: "unknown", box: null }, via: "krawczyk" };
      const center = guess.x.map((r) => toPair(r[0]));
      // start from the spread the interval data allows, and grow
      let radius = Math.max(1e-12, ...A.flat().map((x) => pairToNumber(...ivWidth(x))), ...b.map((x) => pairToNumber(...ivWidth(x))));
      radius *= 1 + Math.max(...center.map((c) => Math.abs(pairToNumber(...c))));
      let last = { status: "unknown", box: null };
      for (let attempt = 0; attempt < 10; attempt += 1) {
        const r = toPair(radius);
        const X = center.map((c) => [addQ(c, [-r[0], r[1]]), addQ(c, r)]);
        last = yield* krawczyk(f, jac, X, { bits });
        if (last.status === "unique") break;
        radius *= 8;
      }
      return { value: { status: last.status === "unique" ? "unique" : "unknown", box: last.box && last.box.map(toCore) }, via: "krawczyk" };
    },
  );
}

const mulPair = ([a, b], [c, d]) => reduce(a * c, b * d);
const lessPair = ([a, b], [c, d]) => a * d < c * b;


// ---------------------------------------------------------------------------
// Similarity structure and the matrix exponential

/**
 * Everything about A up to similarity, exactly:
 * { invariantFactors, minimalPolynomial (Rational coefficients, constant first),
 *   rationalCanonical (Matrix: companion blocks of the invariant factors),
 *   jordan: [{ value, multiplicity, blocks: sizes }] for every eigenvalue,
 *   J, P (Matrices with A = P·J·P⁻¹ when every eigenvalue is rational, else null) }
 * Preview: float eigenvalues only.
 */
export function similarity(input) {
  const A = M(input);
  square(A, "similarity");
  return new Refinement(
    () => ({ value: { eigenvalues: eigenvaluesFloat(A.float(), A.rows) }, trust: "rough" }),
    function* () {
      const e = yield* eigenExact(A.pairs());
      const s = yield* structurePairs(A.pairs(), e);
      const polys = (p) => p.map(toRational);
      return {
        value: {
          invariantFactors: s.invariant.filter((d) => d.length > 1).map(polys),
          minimalPolynomial: polys(s.minimal),
          rationalCanonical: Matrix.fromPairs(s.rational),
          jordan: s.jordan,
          J: s.J && Matrix.fromPairs(s.J),
          P: s.P && Matrix.fromPairs(s.P),
        },
        via: "smith-over-Q[x]",
      };
    },
  );
}

/** Jordan form: same as similarity(A), kept as its own name. */
export const jordan = similarity;

/**
 * exp(A). Preview: floats. Exact stage: a matrix of Core RationalIntervals,
 * each certainly containing the true entry and narrower than 2^-bits.
 */
export function expm(input, { bits = 53 } = {}) {
  const A = M(input);
  square(A, "expm");
  return new Refinement(
    () => ({ value: expmFloat(A.toNumbers()), trust: "good" }),
    function* () {
      const E = yield* expmPairs(A.pairs(), bits);
      return { value: E.map((row) => row.map(toCore)), via: "taylor-interval" };
    },
  );
}
