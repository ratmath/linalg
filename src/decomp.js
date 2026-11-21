// Exact factorizations and other exact constructions. Each public function in
// ops2.js wraps these in a Refinement with a float preview; here is only the
// exact arithmetic, on rows of [n, d] pairs.

import { addQ, subQ, mulQ, divQ, reduce, abs } from "./q.js";

const ZERO = [0n, 1n];
const ONE = [1n, 1n];
const isZero = ([n]) => n === 0n;
const neg = ([n, d]) => [-n, d];

/**
 * PA = LU with L unit lower triangular (m×m) and U upper (m×n), pivoting on
 * the first nonzero entry. Works for any shape and rank.
 * Returns { perm, L, U } where row i of PA is row perm[i] of A.
 */
export function* luPairs(A) {
  const m = A.length;
  const n = A[0]?.length ?? 0;
  const U = A.map((r) => r.map((x) => [...x]));
  const L = Array.from({ length: m }, (_, i) => Array.from({ length: m }, (_, j) => (i === j ? ONE : ZERO)));
  const perm = Array.from({ length: m }, (_, i) => i);
  let r = 0;
  for (let col = 0; col < n && r < m; col += 1) {
    let p = r;
    while (p < m && isZero(U[p][col])) p += 1;
    if (p === m) continue;
    if (p !== r) {
      [U[p], U[r]] = [U[r], U[p]];
      [perm[p], perm[r]] = [perm[r], perm[p]];
      for (let j = 0; j < r; j += 1) [L[p][j], L[r][j]] = [L[r][j], L[p][j]];
    }
    for (let i = r + 1; i < m; i += 1) {
      if (isZero(U[i][col])) continue;
      const f = divQ(U[i][col], U[r][col]);
      L[i][r] = f;
      for (let j = col; j < n; j += 1) U[i][j] = subQ(U[i][j], mulQ(f, U[r][j]));
    }
    r += 1;
    yield { step: "pivot", done: r, of: Math.min(m, n) };
  }
  return { perm, L, U };
}

/**
 * Symmetric A = L·D·Lᵀ, L unit lower, D diagonal, without pivoting.
 * Returns null when a zero pivot blocks it (the matrix then needs pivoting).
 */
export function* ldlPairs(A) {
  const n = A.length;
  const L = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : ZERO)));
  const D = new Array(n).fill(ZERO);
  for (let j = 0; j < n; j += 1) {
    let d = A[j][j];
    for (let k = 0; k < j; k += 1) d = subQ(d, mulQ(mulQ(L[j][k], L[j][k]), D[k]));
    D[j] = d;
    for (let i = j + 1; i < n; i += 1) {
      let s = A[i][j];
      for (let k = 0; k < j; k += 1) s = subQ(s, mulQ(mulQ(L[i][k], L[j][k]), D[k]));
      if (isZero(d)) {
        if (!isZero(s)) return null;
        L[i][j] = ZERO;
      } else L[i][j] = divQ(s, d);
    }
    yield { step: "column", done: j + 1, of: n };
  }
  return { L, D };
}

/**
 * Gram–Schmidt without square roots: A = Q·R with Q's columns mutually
 * orthogonal (not normalised) and R unit upper triangular. Dependent columns
 * give zero columns in Q. normsSquared[j] = ‖q_j‖², so the orthonormal form is
 * Q·diag(1/√n_j) and R scaled by √n_j row-wise.
 */
export function* gramSchmidtPairs(A) {
  const m = A.length;
  const n = A[0]?.length ?? 0;
  const col = (M, j) => M.map((row) => row[j]);
  const dot = (u, v) => u.reduce((s, x, i) => addQ(s, mulQ(x, v[i])), ZERO);
  const Qc = [];
  const N = [];
  const R = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : ZERO)));
  for (let j = 0; j < n; j += 1) {
    let v = col(A, j);
    const a = v;
    for (let k = 0; k < j; k += 1) {
      if (isZero(N[k])) continue;
      const c = divQ(dot(Qc[k], a), N[k]);
      R[k][j] = c;
      if (!isZero(c)) v = v.map((x, i) => subQ(x, mulQ(c, Qc[k][i])));
    }
    Qc.push(v);
    N.push(dot(v, v));
    yield { step: "column", done: j + 1, of: n };
  }
  const Q = Array.from({ length: m }, (_, i) => Qc.map((c) => c[i]));
  return { Q, R, normsSquared: N };
}

export function multiplyPairs(A, B) {
  const m = A.length;
  const k = B.length;
  const n = B[0]?.length ?? 0;
  const out = [];
  for (let i = 0; i < m; i += 1) {
    const row = [];
    for (let j = 0; j < n; j += 1) {
      let num = 0n;
      let den = 1n;
      for (let t = 0; t < k; t += 1) {
        const [a, b] = A[i][t];
        if (a === 0n) continue;
        const [c, d] = B[t][j];
        if (c === 0n) continue;
        const dd = b * d;
        if (dd === den) num += a * c;
        else {
          num = num * dd + a * c * den;
          den *= dd;
        }
      }
      row.push(reduce(num, den));
    }
    out.push(row);
  }
  return out;
}

export const transposePairs = (A) => (A.length ? A[0].map((_, j) => A.map((r) => r[j])) : []);
export const identityPairs = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : ZERO)));

/** A^k for k ≥ 0 by repeated squaring (yields once per squaring). */
export function* powerPairs(A, k) {
  let result = identityPairs(A.length);
  let base = A;
  let e = BigInt(k);
  const steps = e.toString(2).length;
  let done = 0;
  while (e > 0n) {
    if (e & 1n) result = multiplyPairs(result, base);
    e >>= 1n;
    if (e > 0n) base = multiplyPairs(base, base);
    done += 1;
    yield { step: "square", done, of: steps };
  }
  return result;
}

// ---------------------------------------------------------------------------
// Integer normal forms. Entries must be integers (denominator 1).

function requireIntegers(A) {
  return A.map((row) => row.map(([n, d]) => {
    if (d !== 1n) throw new RangeError("Hermite and Smith forms need integer entries");
    return n;
  }));
}

const idInt = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1n : 0n)));

/** Extended gcd: [g, x, y] with a·x + b·y = g ≥ 0. */
function egcd(a, b) {
  let [r0, r1, s0, s1, t0, t1] = [a, b, 1n, 0n, 0n, 1n];
  while (r1 !== 0n) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
    [t0, t1] = [t1, t0 - q * t1];
  }
  return r0 < 0n ? [-r0, -s0, -t0] : [r0, s0, t0];
}

/** Floor division for BigInt. */
const fdiv = (a, b) => {
  const q = a / b;
  return (a % b !== 0n && (a < 0n) !== (b < 0n)) ? q - 1n : q;
};

/**
 * Row-style Hermite normal form: U·A = H with U unimodular, H in row echelon
 * form, pivots positive, entries above each pivot reduced into [0, pivot).
 */
export function hermiteInts(A0) {
  const A = requireIntegers(A0).map((r) => [...r]);
  const m = A.length;
  const n = A[0]?.length ?? 0;
  const U = idInt(m);
  const rowOp = (i, j, a, b, c, d) => {
    // [row_i; row_j] ← [[a, b], [c, d]] · [row_i; row_j]
    for (const M of [A, U]) {
      const ri = M[i];
      const rj = M[j];
      for (let k = 0; k < ri.length; k += 1) {
        const x = ri[k];
        const y = rj[k];
        ri[k] = a * x + b * y;
        rj[k] = c * x + d * y;
      }
    }
  };
  const pivots = [];
  let r = 0;
  for (let col = 0; col < n && r < m; col += 1) {
    for (let i = r + 1; i < m; i += 1) {
      if (A[i][col] === 0n) continue;
      const [g, x, y] = egcd(A[r][col], A[i][col]);
      const a = A[r][col] / g;
      const b = A[i][col] / g;
      rowOp(r, i, x, y, -b, a);
    }
    if (A[r][col] === 0n) continue;
    if (A[r][col] < 0n) for (const M of [A, U]) M[r] = M[r].map((v) => -v);
    for (let i = 0; i < r; i += 1) {
      const q = fdiv(A[i][col], A[r][col]);
      if (q !== 0n) for (const M of [A, U]) M[i] = M[i].map((v, k) => v - q * M[r][k]);
    }
    pivots.push(col);
    r += 1;
  }
  return { H: A, U, pivots };
}

/**
 * Smith normal form: U·A·V = S, S diagonal with d1 | d2 | …, all ≥ 0,
 * U and V unimodular.
 */
export function smithInts(A0) {
  const A = requireIntegers(A0).map((r) => [...r]);
  const m = A.length;
  const n = A[0]?.length ?? 0;
  const U = idInt(m);
  const V = idInt(n);
  const rows = (i, j, a, b, c, d) => {
    for (const M of [A, U]) {
      const ri = M[i];
      const rj = M[j];
      for (let k = 0; k < ri.length; k += 1) {
        const x = ri[k];
        const y = rj[k];
        ri[k] = a * x + b * y;
        rj[k] = c * x + d * y;
      }
    }
  };
  const cols = (i, j, a, b, c, d) => {
    for (const M of [A, V]) {
      for (const row of M) {
        const x = row[i];
        const y = row[j];
        row[i] = a * x + b * y;
        row[j] = c * x + d * y;
      }
    }
  };
  const k = Math.min(m, n);
  for (let t = 0; t < k; t += 1) {
    // bring a nonzero entry of smallest size to (t, t)
    for (;;) {
      let bi = -1;
      let bj = -1;
      for (let i = t; i < m; i += 1) for (let j = t; j < n; j += 1) {
        if (A[i][j] !== 0n && (bi < 0 || abs(A[i][j]) < abs(A[bi][bj]))) { bi = i; bj = j; }
      }
      if (bi < 0) return finish();
      if (bi !== t) rows(t, bi, 0n, 1n, 1n, 0n);
      if (bj !== t) cols(t, bj, 0n, 1n, 1n, 0n);
      let clean = true;
      for (let i = t + 1; i < m; i += 1) {
        if (A[i][t] === 0n) continue;
        const [g, x, y] = egcd(A[t][t], A[i][t]);
        const a = A[t][t] / g;
        const b = A[i][t] / g;
        rows(t, i, x, y, -b, a);
      }
      for (let j = t + 1; j < n; j += 1) {
        if (A[t][j] === 0n) continue;
        const [g, x, y] = egcd(A[t][t], A[t][j]);
        const a = A[t][t] / g;
        const b = A[t][j] / g;
        cols(t, j, x, y, -b, a);
      }
      for (let i = t + 1; i < m; i += 1) if (A[i][t] !== 0n) clean = false;
      if (!clean) continue;
      // divisibility: the pivot must divide every remaining entry
      let bad = -1;
      for (let i = t + 1; i < m && bad < 0; i += 1) for (let j = t + 1; j < n; j += 1) {
        if (A[i][j] % A[t][t] !== 0n) { bad = i; break; }
      }
      if (bad < 0) break;
      rows(t, bad, 1n, 1n, 0n, 1n); // row_t += row_bad, then repeat
    }
  }
  return finish();

  function finish() {
    for (let t = 0; t < k; t += 1) if (A[t][t] < 0n) for (const M of [A, U]) M[t] = M[t].map((v) => -v);
    return { S: A, U, V, diagonal: Array.from({ length: k }, (_, t) => A[t][t]) };
  }
}

export { ZERO, ONE, neg, isZero };
