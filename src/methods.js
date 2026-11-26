// The classical numerical methods of linear algebra, for study.
//
// The rest of this library returns exact answers and avoids these methods on
// purpose. Here they are run as a textbook would run them, in doubles, and
// beside each step stands what exact arithmetic says:
//
//   error       how far the double iterate is from the exact answer
//   methodError how far the same iteration carried out in exact rationals is
//               from the exact answer: the error of the METHOD alone
//   rounding    how far the double iterate is from the exact iterate: the
//               part of the error that rounding added
//
// Every double is a rational, so these distances are computed exactly and
// only then turned into numbers. Exact iterates grow long quickly, so they are
// carried for `exactIterations` steps (and while they stay below a size
// limit); after that the rows give the error only.
//
//   jacobi, gaussSeidel, sor     stationary iterations, with the iteration
//                                matrix and its spectral radius
//   conjugateGradient            exact in at most n steps without rounding
//   powerMethod, inverseIteration  the dominant (or nearest) eigenvalue
//   qrAlgorithm                  all eigenvalues, by repeated QR in doubles
//   pivoting                     elimination with and without row exchanges

import { Matrix } from "./matrix.js";
import { toPair, reduce, addQ, subQ, mulQ, divQ, pairToNumber, abs, bitLength } from "./q.js";
import { eigenvaluesFloat, qrFloat } from "./eigen-float.js";
import { eigenvalues } from "./ops2.js";

const Z = [0n, 1n];
const ONE = [1n, 1n];
const isZ = (p) => p[0] === 0n;
const num = (p) => pairToNumber(p[0], p[1]);
const absQ = (p) => (p[0] < 0n ? [-p[0], p[1]] : p);
const cmpQ = (a, b) => {
  const l = a[0] * b[1];
  const r = b[0] * a[1];
  return l < r ? -1 : l > r ? 1 : 0;
};
const maxQ = (xs) => xs.reduce((m, v) => (cmpQ(v, m) > 0 ? v : m), Z);
const sizeOf = (v) => v.reduce((s, p) => s + bitLength(abs(p[0])) + bitLength(p[1]), 0);

/** The exact rational value of a double. */
export function exactDouble(x) {
  if (!Number.isFinite(x)) throw new RangeError(`${x} has no exact value`);
  if (x === 0) return Z;
  const view = new DataView(new ArrayBuffer(8));
  view.setFloat64(0, x);
  const bits = view.getBigUint64(0);
  const negative = bits >> 63n === 1n;
  const e = Number((bits >> 52n) & 0x7ffn);
  const frac = bits & ((1n << 52n) - 1n);
  const m = e === 0 ? frac : frac | (1n << 52n);
  const p = e === 0 ? -1074 : e - 1075;
  const n = negative ? -m : m;
  return p >= 0 ? [n << BigInt(p), 1n] : reduce(n, 1n << BigInt(-p));
}

const rowsOf = (A) => {
  const M = Matrix.from(A);
  return Array.from({ length: M.rows }, (_, i) => Array.from({ length: M.cols }, (_, j) => M.pair(i, j)));
};
const vecOf = (b, n, what) => {
  const v = Array.isArray(b) ? b.map((x) => (Array.isArray(x) && typeof x[0] !== "bigint" ? toPair(x[0]) : toPair(x))) : null;
  if (!v || v.length !== n) throw new RangeError(`${what} must have ${n} entries`);
  return v;
};
const square = (A, what) => {
  if (!A.length || A.some((r) => r.length !== A.length)) throw new RangeError(`${what} needs a square matrix`);
  return A.length;
};
/** ‖u − v‖∞ exactly, as a pair. */
const distance = (u, v) => maxQ(u.map((p, i) => absQ(subQ(p, v[i]))));
const distanceToFloat = (exactV, floats) => distance(exactV, floats.map(exactDouble));

/** Gaussian elimination with exact rationals; null when singular. */
function solveExact(A, b) {
  const n = A.length;
  const M = A.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    const p = M.findIndex((r, i) => i >= c && !isZ(r[c]));
    if (p < 0) return null;
    [M[c], M[p]] = [M[p], M[c]];
    for (let i = 0; i < n; i++) {
      if (i === c || isZ(M[i][c])) continue;
      const f = divQ(M[i][c], M[c][c]);
      M[i] = M[i].map((v, k) => subQ(v, mulQ(f, M[c][k])));
    }
  }
  return M.map((r, i) => divQ(r[n], r[i]));
}
const matVec = (A, x) => A.map((r) => r.reduce((s, a, j) => (isZ(a) || isZ(x[j]) ? s : addQ(s, mulQ(a, x[j]))), Z));
const fMatVec = (A, x) => A.map((r) => r.reduce((s, a, j) => s + a * x[j], 0));
const fNormInf = (v) => v.reduce((m, x) => Math.max(m, Math.abs(x)), 0);

// ---- stationary iterations ------------------------------------------------------------

function stationary(kind, Ain, bin, opts = {}) {
  const A = rowsOf(Ain);
  const n = square(A, kind);
  const b = vecOf(bin, n, "b");
  const omega = kind === "sor" ? toPair(opts.omega ?? "1") : ONE;
  const omegaF = num(omega);
  const iterations = opts.iterations ?? 25;
  const exactIterations = Math.min(opts.exactIterations ?? 20, iterations);
  const maxBits = opts.maxBits ?? 200000;
  const tolerance = opts.tolerance ?? 0;
  for (let i = 0; i < n; i++) if (isZ(A[i][i])) throw new RangeError(`${kind} divides by the diagonal, and entry (${i + 1}, ${i + 1}) is 0`);
  const solution = solveExact(A, b);
  if (!solution) throw new RangeError("the matrix is singular, so there is no single solution to approach");
  const Af = A.map((r) => r.map(num));
  const bf = b.map(num);
  const x0 = opts.x0 ? vecOf(opts.x0, n, "x0") : A.map(() => Z);
  const step = (x, sum, sub, mul, div, zero, w, combine) => {
    const next = [...x];
    for (let i = 0; i < n; i++) {
      const src = kind === "jacobi" ? x : next;
      let s = b[i] === undefined ? zero : sum.b(i);
      for (let j = 0; j < n; j++) if (j !== i) s = sub(s, mul(sum.a(i, j), src[j]));
      const gs = div(s, sum.a(i, i));
      next[i] = combine(x[i], gs, w);
    }
    return next;
  };
  const exactOps = { b: (i) => b[i], a: (i, j) => A[i][j] };
  const floatOps = { b: (i) => bf[i], a: (i, j) => Af[i][j] };
  const exactCombine = (old, gs, w) => (kind === "sor" ? addQ(mulQ(subQ(ONE, w), old), mulQ(w, gs)) : gs);
  const floatCombine = (old, gs, w) => (kind === "sor" ? (1 - w) * old + w * gs : gs);
  let xe = x0;
  let xf = x0.map(num);
  const rows = [];
  let exactOn = true;
  for (let k = 0; k <= iterations; k++) {
    if (k > 0) {
      xf = step(xf, floatOps, (a, c) => a - c, (a, c) => a * c, (a, c) => a / c, 0, omegaF, floatCombine);
      if (exactOn) xe = step(xe, exactOps, subQ, mulQ, divQ, Z, omega, exactCombine);
    }
    if (exactOn && (k > exactIterations || sizeOf(xe) > maxBits)) exactOn = false;
    const error = num(distanceToFloat(solution, xf));
    const r = fNormInf(bf.map((v, i) => v - fMatVec(Af, xf)[i]));
    rows.push({
      n: k,
      x: [...xf],
      error,
      methodError: exactOn ? num(distance(solution, xe)) : null,
      rounding: exactOn ? num(distanceToFloat(xe, xf)) : null,
      residual: r,
    });
    if (k > 0 && error <= tolerance) break;
  }
  // the iteration matrix T with x_{k+1} = T x_k + c: T = (D + ωL)⁻¹((1 − ω)D − ωU)
  const T = iterationMatrix(kind, A, omega);
  const rho = spectralRadius(T);
  const diagonal = A.every((r, i) => {
    let off = Z;
    r.forEach((a, j) => (j !== i ? (off = addQ(off, absQ(a))) : null));
    return cmpQ(absQ(r[i]), off) > 0;
  });
  return {
    method: kind === "jacobi" ? "Jacobi" : kind === "gauss-seidel" ? "Gauss–Seidel" : `SOR (ω = ${pairStr(omega)})`,
    solution: solution.map(pairStr),
    rows,
    iterationMatrix: T.map((r) => r.map(pairStr)),
    spectralRadius: rho,
    converges: rho < 1,
    strictlyDiagonallyDominant: diagonal,
    /** the observed factor by which the error shrinks per step, which tends to the spectral radius */
    ratio: observedRatio(rows),
  };
}
const pairStr = ([n, d]) => (d === 1n ? `${n}` : `${n}/${d}`);
function observedRatio(rows) {
  const e = rows.map((r) => r.methodError ?? r.error).filter((v) => v > 0);
  if (e.length < 4) return null;
  // over two steps, so that an error that alternates in size (as Jacobi's can) is averaged
  return Math.sqrt(e[e.length - 1] / e[e.length - 3]);
}
function iterationMatrix(kind, A, omega) {
  const n = A.length;
  // columns of T: T e_j, found by one sweep from x = e_j with b = 0
  const cols = [];
  for (let j = 0; j < n; j++) {
    const x = A.map((_, i) => (i === j ? ONE : Z));
    const next = [...x];
    for (let i = 0; i < n; i++) {
      const src = kind === "jacobi" ? x : next;
      let s = Z;
      for (let k = 0; k < n; k++) if (k !== i) s = subQ(s, mulQ(A[i][k], src[k]));
      const gs = divQ(s, A[i][i]);
      next[i] = kind === "sor" ? addQ(mulQ(subQ(ONE, omega), x[i]), mulQ(omega, gs)) : gs;
    }
    cols.push(next);
  }
  return A.map((_, i) => cols.map((c) => c[i]));
}
function spectralRadius(T) {
  const n = T.length;
  const flat = new Float64Array(n * n);
  T.forEach((r, i) => r.forEach((v, j) => (flat[i * n + j] = num(v))));
  const ev = eigenvaluesFloat(flat, n) ?? [];
  return ev.reduce((m, e) => Math.max(m, Math.hypot(e.re, e.im)), 0);
}

/** Jacobi's method: every component updated from the previous iterate. */
export const jacobi = (A, b, opts) => stationary("jacobi", A, b, opts);
/** Gauss–Seidel: each component updated from the newest values. */
export const gaussSeidel = (A, b, opts) => stationary("gauss-seidel", A, b, opts);
/** Successive over-relaxation with factor ω ({ omega: "5/4" }); ω = 1 is Gauss–Seidel. */
export const sor = (A, b, opts) => stationary("sor", A, b, opts);

// ---- conjugate gradients -----------------------------------------------------------------

/**
 * Conjugate gradients for a symmetric positive definite matrix. In exact
 * arithmetic the residuals are orthogonal and the method ends at the exact
 * solution in at most n steps; in doubles it does not, which the rows show.
 */
export function conjugateGradient(Ain, bin, opts = {}) {
  const A = rowsOf(Ain);
  const n = square(A, "conjugateGradient");
  const b = vecOf(bin, n, "b");
  if (A.some((r, i) => r.some((v, j) => cmpQ(v, A[j][i]) !== 0))) throw new RangeError("conjugate gradients needs a symmetric matrix");
  const solution = solveExact(A, b);
  if (!solution) throw new RangeError("the matrix is singular");
  const iterations = opts.iterations ?? n + 3;
  const Af = A.map((r) => r.map(num));
  const dot = (u, v) => u.reduce((s, p, i) => addQ(s, mulQ(p, v[i])), Z);
  const fdot = (u, v) => u.reduce((s, p, i) => s + p * v[i], 0);
  // exact
  let xe = A.map(() => Z);
  let re = [...b];
  let pe = [...re];
  let rre = dot(re, re);
  let exactDone = false;
  // float
  let xf = A.map(() => 0);
  let rf = b.map(num);
  let pf = [...rf];
  let rrf = fdot(rf, rf);
  const rows = [{ n: 0, x: [...xf], error: num(distanceToFloat(solution, xf)), methodError: num(distance(solution, xe)), rounding: 0, residual: Math.sqrt(rrf) }];
  for (let k = 1; k <= iterations; k++) {
    if (!exactDone) {
      if (isZ(rre)) exactDone = true;
      else {
        const Ap = matVec(A, pe);
        const pAp = dot(pe, Ap);
        if (cmpQ(pAp, Z) <= 0) throw new RangeError("the matrix is not positive definite: pᵀAp ≤ 0 for a search direction");
        const alpha = divQ(rre, pAp);
        xe = xe.map((v, i) => addQ(v, mulQ(alpha, pe[i])));
        re = re.map((v, i) => subQ(v, mulQ(alpha, Ap[i])));
        const rr2 = dot(re, re);
        pe = re.map((v, i) => addQ(v, mulQ(divQ(rr2, rre), pe[i])));
        rre = rr2;
      }
    }
    if (rrf > 0) {
      const Ap = fMatVec(Af, pf);
      const alpha = rrf / fdot(pf, Ap);
      xf = xf.map((v, i) => v + alpha * pf[i]);
      rf = rf.map((v, i) => v - alpha * Ap[i]);
      const rr2 = fdot(rf, rf);
      pf = rf.map((v, i) => v + (rr2 / rrf) * pf[i]);
      rrf = rr2;
    }
    rows.push({ n: k, x: [...xf], error: num(distanceToFloat(solution, xf)), methodError: num(distance(solution, xe)), rounding: num(distanceToFloat(xe, xf)), residual: Math.sqrt(rrf) });
  }
  const exactSteps = rows.findIndex((r) => r.methodError === 0);
  return { method: "conjugate gradients", solution: solution.map(pairStr), rows, exactSteps };
}

// ---- eigenvalues ------------------------------------------------------------------------------

function trueEigenvalues(A) {
  const e = eigenvalues(Matrix.fromPairs(A)).exact();
  const out = [];
  for (const r of e.real) for (let k = 0; k < r.multiplicity; k++) out.push({ re: r.value.toNumber(), im: 0 });
  for (const c of e.complex) {
    const { re, im } = c.value.toNumbers();
    for (let k = 0; k < c.multiplicity; k++) out.push({ re, im });
  }
  return out;
}

/**
 * The power method: x ← Ax, scaled so its largest entry is 1. The Rayleigh
 * quotient xᵀAx/xᵀx estimates the dominant eigenvalue, and the error falls
 * like |λ₂/λ₁|ᵏ (like its square for a symmetric matrix).
 */
export function powerMethod(Ain, opts = {}) {
  const A = rowsOf(Ain);
  const n = square(A, "powerMethod");
  return eigenIteration("power", A, n, opts, null);
}

/**
 * Inverse iteration with a shift: the power method on (A − σI)⁻¹, which finds
 * the eigenvalue nearest σ. { shift: "3/2" }.
 */
export function inverseIteration(Ain, opts = {}) {
  const A = rowsOf(Ain);
  const n = square(A, "inverseIteration");
  const sigma = toPair(opts.shift ?? 0);
  return eigenIteration("inverse", A, n, opts, sigma);
}

function eigenIteration(kind, A, n, opts, sigma) {
  const iterations = opts.iterations ?? 25;
  const exactIterations = Math.min(opts.exactIterations ?? 12, iterations);
  const maxBits = opts.maxBits ?? 100000;
  const truth = trueEigenvalues(A);
  const sigmaF = sigma ? num(sigma) : 0;
  const key = (e) => (kind === "power" ? Math.hypot(e.re, e.im) : -Math.hypot(e.re - sigmaF, e.im));
  const sorted = [...truth].sort((a, b) => key(b) - key(a));
  const target = sorted[0];
  const second = sorted.find((e) => Math.abs(e.re - target.re) + Math.abs(e.im - target.im) > 1e-12);
  const dominant = target.im === 0 && !(second && Math.abs(key(second) - key(target)) < 1e-12);
  const shifted = sigma ? A.map((r, i) => r.map((v, j) => (i === j ? subQ(v, sigma) : v))) : A;
  const Sf = shifted.map((r) => r.map(num));
  if (kind === "inverse" && !solveExact(shifted, A.map(() => ONE)) && !solveExact(shifted, A.map((_, i) => (i === 0 ? ONE : Z)))) {
    throw new RangeError("the shift is an eigenvalue, so A − σI cannot be inverted");
  }
  const applyExact = (x) => (kind === "power" ? matVec(A, x) : solveExact(shifted, x));
  const applyFloat = (x) => (kind === "power" ? fMatVec(Sf, x) : floatSolve(Sf, x));
  const scaleExact = (v) => {
    const m = v.reduce((best, p) => (cmpQ(absQ(p), absQ(best)) > 0 ? p : best), Z);
    return isZ(m) ? v : v.map((p) => divQ(p, m));
  };
  const scaleFloat = (v) => {
    const m = v.reduce((best, p) => (Math.abs(p) > Math.abs(best) ? p : best), 0);
    return m === 0 ? v : v.map((p) => p / m);
  };
  const rayleighE = (x) => {
    const Ax = matVec(A, x);
    const xx = x.reduce((s, p) => addQ(s, mulQ(p, p)), Z);
    return divQ(x.reduce((s, p, i) => addQ(s, mulQ(p, Ax[i])), Z), xx);
  };
  const Af = A.map((r) => r.map(num));
  const rayleighF = (x) => {
    const Ax = fMatVec(Af, x);
    return x.reduce((s, p, i) => s + p * Ax[i], 0) / x.reduce((s, p) => s + p * p, 0);
  };
  let xe = opts.x0 ? vecOf(opts.x0, n, "x0") : A.map(() => ONE);
  let xf = xe.map(num);
  let exactOn = true;
  const rows = [];
  for (let k = 0; k <= iterations; k++) {
    if (k > 0) {
      xf = scaleFloat(applyFloat(xf));
      if (exactOn) xe = scaleExact(applyExact(xe));
    }
    if (exactOn && (k > exactIterations || sizeOf(xe) > maxBits)) exactOn = false;
    const lf = rayleighF(xf);
    const le = exactOn ? rayleighE(xe) : null;
    rows.push({
      n: k,
      x: [...xf],
      estimate: lf,
      exactEstimate: le && pairStr(le),
      error: dominant ? Math.abs(lf - target.re) : null,
      methodError: dominant && le ? Math.abs(num(le) - target.re) : null,
    });
  }
  const ratio = dominant && second ? (kind === "power" ? Math.hypot(second.re, second.im) / Math.abs(target.re) : Math.abs(target.re - sigmaF) / Math.hypot(second.re - sigmaF, second.im)) : null;
  return {
    method: kind === "power" ? "power method" : `inverse iteration (shift ${pairStr(sigma)})`,
    eigenvalues: truth,
    target: dominant ? target.re : null,
    note: dominant ? null : "no single real eigenvalue is strictly largest (or nearest the shift), so the iteration need not settle",
    convergenceFactor: ratio,
    rows,
  };
}
function floatSolve(M, b) {
  const n = M.length;
  const A = M.map((r, i) => [...r, b[i]]);
  for (let c = 0; c < n; c++) {
    let p = c;
    for (let i = c + 1; i < n; i++) if (Math.abs(A[i][c]) > Math.abs(A[p][c])) p = i;
    [A[c], A[p]] = [A[p], A[c]];
    for (let i = c + 1; i < n; i++) {
      const f = A[i][c] / A[c][c];
      for (let k = c; k <= n; k++) A[i][k] -= f * A[c][k];
    }
  }
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    let s = A[i][n];
    for (let k = i + 1; k < n; k++) s -= A[i][k] * x[k];
    x[i] = s / A[i][i];
  }
  return x;
}

/**
 * The unshifted QR algorithm in doubles: A₀ = A, Aₖ = QₖRₖ, Aₖ₊₁ = RₖQₖ. Each
 * Aₖ is similar to A; the entries below the diagonal shrink (like
 * |λᵢ₊₁/λᵢ|ᵏ) and the diagonal tends to the eigenvalues, which are given
 * exactly for comparison. Complex pairs stay as 2×2 blocks.
 */
export function qrAlgorithm(Ain, opts = {}) {
  const A = rowsOf(Ain);
  const n = square(A, "qrAlgorithm");
  const iterations = opts.iterations ?? 30;
  const truth = trueEigenvalues(A);
  let M = new Float64Array(n * n);
  A.forEach((r, i) => r.forEach((v, j) => (M[i * n + j] = num(v))));
  const rows = [];
  for (let k = 0; k <= iterations; k++) {
    if (k > 0) {
      const { Q, R } = qrFloat(M, n, n);
      const next = new Float64Array(n * n);
      for (let i = 0; i < n; i++)
        for (let j = 0; j < n; j++) {
          let s = 0;
          for (let t = 0; t < n; t++) s += R[i][t] * Q[t][j];
          next[i * n + j] = s;
        }
      M = next;
    }
    let below = 0;
    for (let i = 1; i < n; i++) below = Math.max(below, Math.abs(M[i * n + i - 1]));
    rows.push({ n: k, diagonal: Array.from({ length: n }, (_, i) => M[i * n + i]), subdiagonal: below });
  }
  return { method: "QR algorithm (unshifted)", eigenvalues: truth, rows, final: Array.from({ length: n }, (_, i) => Array.from(M.subarray(i * n, i * n + n))) };
}

// ---- pivoting -----------------------------------------------------------------------------------

/**
 * Gaussian elimination in doubles, without row exchanges and with partial
 * pivoting, beside the exact solution: { exact, none, partial }, each with
 * { x, error, growth } where growth is the largest entry met during
 * elimination over the largest entry of A. A tiny pivot makes the
 * multipliers, the growth and the error large.
 */
export function pivoting(Ain, bin) {
  const A = rowsOf(Ain);
  const n = square(A, "pivoting");
  const b = vecOf(bin, n, "b");
  const solution = solveExact(A, b);
  if (!solution) throw new RangeError("the matrix is singular");
  const run = (partial) => {
    const M = A.map((r, i) => [...r.map(num), num(b[i])]);
    const start = Math.max(...A.flat().map((p) => Math.abs(num(p))));
    let big = start;
    for (let c = 0; c < n; c++) {
      let p = c;
      if (partial) for (let i = c + 1; i < n; i++) if (Math.abs(M[i][c]) > Math.abs(M[p][c])) p = i;
      if (M[p][c] === 0) for (let i = c + 1; i < n && M[p][c] === 0; i++) p = i;
      [M[c], M[p]] = [M[p], M[c]];
      for (let i = c + 1; i < n; i++) {
        const f = M[i][c] / M[c][c];
        for (let k = c; k <= n; k++) {
          M[i][k] -= f * M[c][k];
          if (k < n) big = Math.max(big, Math.abs(M[i][k]));
        }
      }
    }
    const x = new Array(n).fill(0);
    for (let i = n - 1; i >= 0; i--) {
      let s = M[i][n];
      for (let k = i + 1; k < n; k++) s -= M[i][k] * x[k];
      x[i] = s / M[i][i];
    }
    const finite = x.every(Number.isFinite);
    return { x, error: finite ? num(distanceToFloat(solution, x)) : Infinity, growth: big / start };
  };
  return { exact: solution.map(pairStr), none: run(false), partial: run(true) };
}
