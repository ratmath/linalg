// Exact results that leave the rationals: eigenspaces over Q(λ), Jordan
// forms with algebraic eigenvalues, norms and condition numbers, Cholesky
// with exact square roots, orthonormal eigenvectors of symmetric matrices,
// and lattice tools on integer matrices.
//
// A number here is an exact algebraic number from @ratmath/algebraic
// (AlgebraicReal or ComplexAlgebraic): a polynomial with integer
// coefficients plus a region that holds exactly one of its roots. "Algebraic
// matrix" means an array of rows of such numbers.

import * as ALG from "@ratmath/algebraic";
import { NumberField, lll } from "@ratmath/algebraic";
import { Matrix } from "./matrix.js";
import { Refinement } from "@ratmath/core/refine";
import { AlgebraicVector } from "./structure.js";
import { eigenExact, singularExact } from "./spectral.js";
import { hermiteInts } from "./decomp.js";
import { svdFloat } from "./eigen-float.js";
import { rref as fieldRref } from "./field.js";
import { toRational } from "./q.js";
import { inverse } from "./ops.js";
import { charpolyModular, fromPairs } from "./poly.js";

const M = (x) => (x instanceof Matrix ? x : Matrix.from(x));
const Q = (p) => ALG.fromRational(p);
const ZERO = [0n, 1n];

/** The operations of a NumberField as the plain record field.js works with. */
export function numberFieldOps(K) {
  return {
    zero: K.zero,
    one: K.one,
    add: (a, b) => a.add(b),
    sub: (a, b) => a.sub(b),
    mul: (a, b) => a.mul(b),
    div: (a, b) => a.div(b),
    neg: (a) => a.neg(),
    sign: (a) => (a.isZero() ? 0 : 1),
  };
}

/** A basis of the null space from a reduced row echelon form. */
function nullBasis(rows, pivots, n, F) {
  const basis = [];
  const isPivot = new Set(pivots);
  for (let free = 0; free < n; free += 1) {
    if (isPivot.has(free)) continue;
    const v = Array.from({ length: n }, () => F.zero);
    v[free] = F.one;
    pivots.forEach((col, r) => {
      v[col] = F.neg(rows[r][free]);
    });
    basis.push(v);
  }
  return basis;
}

/** Q(λ) for an eigenvalue, with λ as its generator. */
function fieldOf(value) {
  const m = ALG.minimalPolynomial(value);
  if (m.length === 2) return null; // rational
  return new NumberField(m, { name: "λ", root: value, check: false });
}

/** A − λI with entries in K. */
function shifted(pairs, K) {
  return pairs.map((row, i) => row.map((x, j) => (i === j ? K.of(x).sub(K.gen) : K.of(x))));
}

/**
 * A basis of the eigenspace of an irrational or complex eigenvalue, as
 * vectors with entries in Q(λ) (polynomials in λ). Works for any geometric
 * multiplicity: the matrix A − λI is row reduced over the field Q(λ).
 */
export function eigenspace(pairs, value) {
  const K = fieldOf(value);
  if (!K) throw new RangeError("eigenspace: the eigenvalue is rational; use the rational null space");
  const F = numberFieldOps(K);
  const n = pairs.length;
  const { rows, pivots } = fieldRref(shifted(pairs, K), F);
  return nullBasis(rows, pivots, n, F).map((v) => new AlgebraicVector(value, v.map((e) => padded(e.c, K.n))));
}

const padded = (c, n) => Array.from({ length: Math.max(1, Math.min(n, c.length)) }, (_, i) => c[i] ?? ZERO);

// ---- Jordan form over algebraic eigenvalues ------------------------------------------

const matVec = (A, v, F) => A.map((row) => row.reduce((s, x, j) => F.add(s, F.mul(x, v[j])), F.zero));
const matMul = (A, B, F) => A.map((row) => B[0].map((_, j) => row.reduce((s, x, k) => F.add(s, F.mul(x, B[k][j])), F.zero)));
const rankOf = (vectors, F) => (vectors.length ? fieldRref(vectors, F).pivots.length : 0);

/** Jordan chains of N = A − λI over a field: [{ size, chain: [v1, …, vk] }] with N·v1 = 0, N·v(i+1) = v(i). */
function chainsOver(N, multiplicity, F) {
  const n = N.length;
  const kernels = [[]];
  let Nk = N;
  for (let k = 1; ; k += 1) {
    const { rows, pivots } = fieldRref(Nk, F);
    const basis = nullBasis(rows, pivots, n, F);
    kernels.push(basis);
    if (basis.length >= multiplicity || k > n) break;
    Nk = matMul(Nk, N, F);
  }
  const chains = [];
  for (let k = kernels.length - 1; k >= 1; k -= 1) {
    const level = chains.flatMap((c) => (c.size > k ? [c.chain[k - 1]] : []));
    const span = [...kernels[k - 1], ...level];
    let r = rankOf(span, F);
    for (const v of kernels[k]) {
      const r2 = rankOf([...span, v], F);
      if (r2 === r) continue;
      span.push(v);
      r = r2;
      const chain = [v];
      for (let i = 1; i < k; i += 1) chain.unshift(matVec(N, chain[0], F));
      chains.push({ size: k, chain });
    }
  }
  return chains.sort((a, b) => b.size - a.size);
}

/**
 * J and P with A·P = P·J when some eigenvalue is irrational or complex.
 * Entries are exact algebraic numbers. `eigen` is the result of eigenExact.
 * Returns null when a complex eigenvalue could not be certified.
 */
export function jordanAlgebraic(pairs, eigen) {
  if (!eigen.complexCertified) return null;
  const n = pairs.length;
  const blocks = [];
  const columns = [];
  const rationalOps = { zero: ZERO, one: [1n, 1n] };
  for (const { value, multiplicity } of [...eigen.real, ...eigen.complex]) {
    const K = fieldOf(value) ?? new NumberField([0n, 1n], { check: false });
    const F = numberFieldOps(K);
    const lam = K.n === 1 ? K.of(value.rational ?? ALG.rationalOf(value)) : K.gen;
    const N = pairs.map((row, i) => row.map((x, j) => (i === j ? K.of(x).sub(lam) : K.of(x))));
    const out = (e) => (K.n === 1 ? Q(e.c[0] ?? rationalOps.zero) : e.toAlgebraic());
    for (const { size, chain } of chainsOver(N, multiplicity, F)) {
      blocks.push({ value: ALG.simplify(value), size });
      for (const v of chain) columns.push(v.map(out));
    }
  }
  if (columns.length !== n) return null;
  const zero = Q(ZERO);
  const one = Q([1n, 1n]);
  const J = Array.from({ length: n }, () => Array.from({ length: n }, () => zero));
  let at = 0;
  for (const { value, size } of blocks) {
    for (let i = 0; i < size; i += 1) {
      J[at + i][at + i] = value;
      if (i + 1 < size) J[at + i][at + i + 1] = one;
    }
    at += size;
  }
  const P = Array.from({ length: n }, (_, i) => columns.map((c) => c[i]));
  return { J, P };
}

// ---- norms and condition numbers ------------------------------------------------------

const absQ = ([a, b]) => [a < 0n ? -a : a, b];
const addQ = ([a, b], [c, d]) => [a * d + c * b, b * d];
const ltQ = ([a, b], [c, d]) => a * d < c * b;

function exactNorm(A, p, singular) {
  const pairs = A.pairs();
  if (p === 1 || p === "inf") {
    const lines = p === 1 ? pairs[0].map((_, j) => pairs.map((row) => row[j])) : pairs;
    let best = ZERO;
    for (const line of lines) {
      const s = line.reduce((t, x) => addQ(t, absQ(x)), ZERO);
      if (ltQ(best, s)) best = s;
    }
    return Q(best);
  }
  if (p === "fro") {
    const s = pairs.flat().reduce((t, [a, b]) => addQ(t, [a * a, b * b]), ZERO);
    return ALG.sqrt(Q(s));
  }
  if (p === 2) return singular[0]?.value ?? Q(ZERO);
  throw new RangeError('norm: p must be 1, 2, "inf" or "fro"');
}

function floatNorm(A, p) {
  const a = A.toNumbers();
  if (p === 1) return Math.max(...a[0].map((_, j) => a.reduce((s, row) => s + Math.abs(row[j]), 0)));
  if (p === "inf") return Math.max(...a.map((row) => row.reduce((s, x) => s + Math.abs(x), 0)));
  if (p === "fro") return Math.hypot(...a.flat());
  if (p === 2) return svdFloat(A.float(), A.rows, A.cols).s[0] ?? 0;
  throw new RangeError('norm: p must be 1, 2, "inf" or "fro"');
}

/**
 * A matrix norm, exactly: p = 1 (largest column sum) and "inf" (largest row
 * sum) are rationals, "fro" is the square root of a rational, and 2 is the
 * largest singular value. The value is an algebraic number in every case
 * (`.rational` gives the pair when it is rational).
 */
export function norm(input, p = 2) {
  const A = M(input);
  return new Refinement(
    () => ({ value: floatNorm(A, p), trust: "good" }),
    function* () {
      const singular = p === 2 ? yield* singularExact(A.pairs()) : null;
      return { value: exactNorm(A, p, singular), via: p === 2 ? "singular-values" : "exact" };
    },
  );
}

/**
 * The condition number ‖A‖·‖A⁻¹‖ in the same norms, exactly. For p = 2 it
 * is σmax/σmin. A singular matrix gives Infinity.
 */
export function cond(input, p = 2) {
  const A = M(input);
  if (!A.isSquare) throw new RangeError(`cond needs a square matrix, got ${A.rows}×${A.cols}`);
  return new Refinement(
    () => {
      const s = svdFloat(A.float(), A.rows, A.cols).s;
      if (p === 2) return { value: s[s.length - 1] > 0 ? s[0] / s[s.length - 1] : Infinity, trust: "rough" };
      return { value: null, trust: "unsure" };
    },
    function* () {
      if (p === 2) {
        const singular = yield* singularExact(A.pairs());
        const last = singular[singular.length - 1].value;
        if (ALG.isZero(last)) return { value: Infinity, via: "singular-values" };
        return { value: ALG.div(singular[0].value, last), via: "singular-values" };
      }
      const inv = yield* inverseOf(A);
      if (!inv) return { value: Infinity, via: "inverse" };
      return { value: ALG.mul(exactNorm(A, p), exactNorm(inv, p)), via: "inverse" };
    },
  );
}

function* inverseOf(A) {
  for (const s of inverse(A).stages()) {
    if (s.stage === "progress") yield s;
    else if (s.exact) return s.value;
  }
  return null;
}

// ---- symmetric matrices -----------------------------------------------------------------

function symmetric(A, name) {
  if (!A.isSquare || !A.equals(A.transpose())) throw new RangeError(`${name} needs a symmetric matrix`);
}

/** Sign changes in a coefficient list (zeros skipped). */
function variations(coeffs) {
  let count = 0;
  let last = 0;
  for (const c of coeffs) {
    const s = c > 0n ? 1 : c < 0n ? -1 : 0;
    if (s === 0) continue;
    if (last !== 0 && s !== last) count += 1;
    last = s;
  }
  return count;
}

/**
 * Definiteness of a symmetric matrix from the signs of its eigenvalues:
 * { kind, positive, negative, zero } with kind "positive-definite",
 * "positive-semidefinite", "negative-definite", "negative-semidefinite",
 * "indefinite" or "zero". Every eigenvalue of a symmetric matrix is real, so
 * Descartes' rule of signs on the characteristic polynomial counts the
 * positive ones exactly; no eigenvalue is computed.
 */
export function definiteness(input) {
  const A = M(input);
  symmetric(A, "definiteness");
  const n = A.rows;
  return new Refinement(
    () => ({ value: null, trust: "unsure" }),
    function* () {
      const cp = fromPairs(yield* charpolyModular(A.pairs())); // integer coefficients, constant first
      let zero = 0;
      while (zero < cp.length - 1 && cp[zero] === 0n) zero += 1;
      const rest = cp.slice(zero);
      const positive = variations(rest);
      const negative = n - zero - positive;
      const kind = n === zero ? "zero"
        : negative === 0 ? (zero ? "positive-semidefinite" : "positive-definite")
        : positive === 0 ? (zero ? "negative-semidefinite" : "negative-definite")
        : "indefinite";
      return { value: { kind, positive, negative, zero }, via: "descartes" };
    },
  );
}

/**
 * Cholesky factorization A = L·Lᵀ of a symmetric positive definite matrix,
 * with exact square roots. Value: { L } as an algebraic matrix (lower
 * triangular, positive diagonal), or null when A is not positive definite.
 * Preview: L in floats.
 */
export function cholesky(input) {
  const A = M(input);
  symmetric(A, "cholesky");
  const n = A.rows;
  return new Refinement(
    () => {
      const a = A.toNumbers();
      const L = Array.from({ length: n }, () => new Array(n).fill(0));
      for (let j = 0; j < n; j += 1) {
        let d = a[j][j];
        for (let k = 0; k < j; k += 1) d -= L[j][k] ** 2;
        if (!(d > 0)) return { value: null, trust: "unsure" };
        L[j][j] = Math.sqrt(d);
        for (let i = j + 1; i < n; i += 1) {
          let s = a[i][j];
          for (let k = 0; k < j; k += 1) s -= L[i][k] * L[j][k];
          L[i][j] = s / L[j][j];
        }
      }
      return { value: { L }, trust: "good" };
    },
    function* () {
      // rational LDLᵀ first: A = L1·D·L1ᵀ, then L = L1·√D
      const a = A.pairs();
      const sub = ([p, q], [r, s]) => [p * s - r * q, q * s];
      const mul = ([p, q], [r, s]) => [p * r, q * s];
      const div = ([p, q], [r, s]) => (s * r < 0n ? [-p * s, -q * r] : [p * s, q * r]);
      const L1 = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? [1n, 1n] : ZERO)));
      const D = [];
      for (let j = 0; j < n; j += 1) {
        let d = a[j][j];
        for (let k = 0; k < j; k += 1) d = sub(d, mul(mul(L1[j][k], L1[j][k]), D[k]));
        if (d[0] * d[1] <= 0n) return { value: null, via: "ldl" };
        D[j] = d;
        for (let i = j + 1; i < n; i += 1) {
          let s = a[i][j];
          for (let k = 0; k < j; k += 1) s = sub(s, mul(mul(L1[i][k], L1[j][k]), D[k]));
          L1[i][j] = div(s, d);
        }
        yield { step: "column", done: j + 1, of: n };
      }
      const roots = D.map((d) => ALG.sqrt(Q(toPairReduced(d))));
      const L = L1.map((row, i) => row.map((x, j) => (j > i || x[0] === 0n ? Q(ZERO) : ALG.mul(Q(toPairReduced(x)), roots[j]))));
      return { value: { L, D: D.map((d) => toRational(toPairReduced(d))) }, via: "ldl" };
    },
  );
}

function toPairReduced([a, b]) {
  if (b < 0n) [a, b] = [-a, -b];
  const g = gcd(a < 0n ? -a : a, b);
  return [a / g, b / g];
}
const gcd = (a, b) => {
  while (b) [a, b] = [b, a % b];
  return a || 1n;
};

/**
 * Orthogonal diagonalization of a symmetric matrix: A = Q·diag(values)·Qᵀ
 * with QᵀQ = I, exactly. Value: { values, Q } where values are the
 * eigenvalues (ascending, repeated by multiplicity) and Q is an algebraic
 * matrix whose columns are orthonormal eigenvectors. Each eigenspace is
 * found over Q(λ), made orthogonal there by Gram–Schmidt, and normalised
 * with an exact square root.
 */
export function symmetricEigen(input) {
  const A = M(input);
  symmetric(A, "symmetricEigen");
  const n = A.rows;
  return new Refinement(
    () => ({ value: null, trust: "unsure" }),
    function* () {
      const e = yield* eigenExact(A.pairs());
      const values = [];
      const columns = [];
      for (const { value, multiplicity } of e.real) {
        const K = fieldOf(value) ?? new NumberField([0n, 1n], { check: false });
        const F = numberFieldOps(K);
        const lam = K.n === 1 ? K.of(value.rational) : K.gen;
        const N = A.pairs().map((row, i) => row.map((x, j) => (i === j ? K.of(x).sub(lam) : K.of(x))));
        const { rows, pivots } = fieldRref(N, F);
        const basis = nullBasis(rows, pivots, n, F);
        if (basis.length !== multiplicity) throw new Error("symmetricEigen: eigenspace has the wrong dimension");
        const dot = (u, v) => u.reduce((s, x, i) => s.add(x.mul(v[i])), K.zero);
        const ortho = [];
        for (let v of basis) {
          for (const u of ortho) {
            const c = dot(v, u).div(dot(u, u));
            v = v.map((x, i) => x.sub(c.mul(u[i])));
          }
          ortho.push(v);
        }
        const out = (x) => (K.n === 1 ? Q(x.c[0] ?? ZERO) : x.toAlgebraic());
        for (const v of ortho) {
          const length = ALG.sqrt(out(dot(v, v)));
          columns.push(v.map((x) => (x.isZero() ? Q(ZERO) : ALG.div(out(x), length))));
          values.push(ALG.simplify(value));
        }
        yield { step: "eigenspace", done: columns.length, of: n };
      }
      const Qm = Array.from({ length: n }, (_, i) => columns.map((c) => c[i]));
      return { value: { values, Q: Qm }, via: "number-field" };
    },
  );
}

// ---- integer matrices -------------------------------------------------------------------

const intMatrix = (rows) => Matrix.fromPairs(rows.map((row) => row.map((x) => [x, 1n])));

/**
 * A basis of the integer solutions of A·x = 0, as the columns of a Matrix
 * (LLL-reduced, so the vectors are short). Every integer solution is an
 * integer combination of them. Integer entries only.
 */
export function integerKernel(input) {
  const A = M(input);
  const n = A.cols;
  // Hermite form of the transpose: U·Aᵀ = H; rows of U with a zero row of H span the kernel over Z.
  const { H, U } = hermiteInts(A.transpose().pairs());
  const kernel = U.filter((_, i) => H[i].every((x) => x === 0n));
  if (!kernel.length) return Matrix.zeros(n, 0);
  const short = lll(kernel);
  return intMatrix(short).transpose();
}

/**
 * A basis of the lattice spanned by the columns of A (all integer
 * combinations of them), as the columns of a Matrix, LLL-reduced.
 */
export function integerImage(input) {
  const A = M(input);
  const { H } = hermiteInts(A.transpose().pairs());
  const rows = H.filter((row) => row.some((x) => x !== 0n));
  if (!rows.length) return Matrix.zeros(A.rows, 0);
  return intMatrix(lll(rows)).transpose();
}

/** LLL-reduce a lattice basis given as the rows of an integer matrix: short, nearly orthogonal rows spanning the same lattice. */
export function lllReduce(input) {
  const A = M(input);
  const rows = A.pairs().map((row) => row.map(([a, b]) => {
    if (b !== 1n) throw new RangeError("lllReduce needs integer entries");
    return a;
  }));
  return intMatrix(lll(rows));
}

// ---- matrices of algebraic numbers ---------------------------------------------------------

/** Exact algebraic numbers, real or complex, as a field for the functions in `field`. */
export const complexField = {
  get zero() { return Q(ZERO); },
  get one() { return Q([1n, 1n]); },
  add: (a, b) => ALG.add(a, b),
  sub: (a, b) => ALG.sub(a, b),
  mul: (a, b) => ALG.mul(a, b),
  div: (a, b) => ALG.div(a, b),
  neg: (a) => ALG.neg(a),
  sign: (a) => (ALG.isZero(a) ? 0 : 1),
};

/** Rows of exact numbers from rationals ("3/4", 2, 5n), Core Rationals, or algebraic numbers already built. */
export function algebraicMatrix(rows) {
  return rows.map((row) => row.map((x) => (x instanceof ALG.AlgebraicReal || x instanceof ALG.ComplexAlgebraic ? x : ALG.fromRational(x))));
}
