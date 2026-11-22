// Similarity structure: invariant factors (Smith form of xI − A over Q[x]),
// minimal polynomial, rational canonical form, Jordan structure for every
// eigenvalue, the Jordan form with its transform when the eigenvalues are
// rational, and exact eigenvectors over Q(λ) for irrational eigenvalues.

import { Matrix } from "./matrix.js";
import { addQ, subQ, mulQ, divQ, reduce, toRational, pairToNumber } from "./q.js";
import {
  trimQ, degQ, isZeroPoly, addP, subP, mulP, divmodP, modP, monicP, gcdP, toIntPoly, xP, constP, formatP, evalP,
} from "./qpoly.js";
import { polyGcd, divide, degree, squareFree, isolate } from "./poly.js";
import { certifyComplexRoots, hornerBox, ComplexAlgebraic } from "./complex.js";
import { ivMul, ivAdd } from "./interval.js";
import { run, rrefPairs, nullBasis } from "./exact.js";
import { multiplyPairs } from "./decomp.js";

const Z = [0n, 1n];
const ONE = [1n, 1n];

/** Smith form of xI − A over Q[x]: monic invariant factors d₁ | d₂ | … | dₙ. */
export function* invariantFactorsPairs(A) {
  const n = A.length;
  const M = A.map((row, i) => row.map((a, j) => (i === j ? trimQ([[-a[0], a[1]], ONE]) : trimQ([[-a[0], a[1]]]))));
  const swapRows = (i, j) => ([M[i], M[j]] = [M[j], M[i]]);
  const swapCols = (i, j) => M.forEach((row) => ([row[i], row[j]] = [row[j], row[i]]));
  for (let t = 0; t < n; t += 1) {
    for (;;) {
      let bi = -1;
      let bj = -1;
      for (let i = t; i < n; i += 1) for (let j = t; j < n; j += 1) {
        if (!isZeroPoly(M[i][j]) && (bi < 0 || degQ(M[i][j]) < degQ(M[bi][bj]))) { bi = i; bj = j; }
      }
      if (bi < 0) break;
      if (bi !== t) swapRows(t, bi);
      if (bj !== t) swapCols(t, bj);
      const p = M[t][t];
      let dirty = false;
      for (let i = t + 1; i < n; i += 1) {
        if (isZeroPoly(M[i][t])) continue;
        const [q] = divmodP(M[i][t], p);
        for (let j = t; j < n; j += 1) M[i][j] = subP(M[i][j], mulP(q, M[t][j]));
        if (!isZeroPoly(M[i][t])) dirty = true;
      }
      for (let j = t + 1; j < n; j += 1) {
        if (isZeroPoly(M[t][j])) continue;
        const [q] = divmodP(M[t][j], p);
        for (let i = t; i < n; i += 1) M[i][j] = subP(M[i][j], mulP(q, M[i][t]));
        if (!isZeroPoly(M[t][j])) dirty = true;
      }
      if (dirty) continue;
      let bad = -1;
      for (let i = t + 1; i < n && bad < 0; i += 1) for (let j = t + 1; j < n; j += 1) {
        if (!isZeroPoly(modP(M[i][j], p))) { bad = i; break; }
      }
      if (bad < 0) break;
      for (let j = t; j < n; j += 1) M[t][j] = addP(M[t][j], M[bad][j]);
    }
    yield { step: "smith", done: t + 1, of: n };
  }
  return M.map((row, i) => monicP(row[i]));
}

/** Companion matrix of a monic polynomial (constant first). */
export function companionPairs(p) {
  const k = p.length - 1;
  return Array.from({ length: k }, (_, i) => Array.from({ length: k }, (_, j) => {
    if (j === k - 1) return [-p[i][0], p[i][1]];
    return i === j + 1 ? ONE : Z;
  }));
}

function blockDiag(blocks) {
  const n = blocks.reduce((s, b) => s + b.length, 0);
  const out = Array.from({ length: n }, () => new Array(n).fill(Z));
  let o = 0;
  for (const b of blocks) {
    b.forEach((row, i) => row.forEach((x, j) => (out[o + i][o + j] = x)));
    o += b.length;
  }
  return out;
}

/** Real-interval Horner of an integer polynomial over [lo, hi]. */
function hornerReal(p, X) {
  let acc = [Z, Z];
  for (let i = p.length - 1; i >= 0; i -= 1) acc = ivAdd(ivMul(acc, X), [[p[i], 1n], [p[i], 1n]]);
  return acc;
}
const excludesZero = ([lo, hi]) => lo[0] > 0n || hi[0] < 0n;

/** Is the algebraic number `value` (AlgebraicReal or ComplexAlgebraic) a root of the integer polynomial h? Exact. */
export function isRootOf(value, h) {
  const g = polyGcd(h, value.poly);
  if (degree(g) < 1) return false;
  if (degree(g) === degree(value.poly)) return true;
  const q = divide(value.poly, g);
  for (let bits = 16; bits < 1 << 14; bits *= 2) {
    value.refine(bits);
    if (value instanceof ComplexAlgebraic) {
      const box = [value.re, value.im];
      const [gr, gi] = hornerBox(g, box);
      if (excludesZero(gr) || excludesZero(gi)) return false;
      const [qr, qi] = hornerBox(q, box);
      if (excludesZero(qr) || excludesZero(qi)) return true;
    } else {
      if (value.isRational) return evalP(g.map((c) => [c, 1n]), value.rational)[0] === 0n;
      const X = [value.lo, value.hi];
      if (excludesZero(hornerReal(g, X))) return false;
      if (excludesZero(hornerReal(q, X))) return true;
    }
  }
  throw new Error("isRootOf: could not separate (should not happen for square-free input)");
}

/**
 * Jordan structure from invariant factors: for each eigenvalue, the sizes of
 * its Jordan blocks. `eigen` is the list [{ value, multiplicity }] (real and
 * complex together).
 */
export function jordanStructure(invariant, eigen) {
  const parts = invariant
    .filter((d) => degQ(d) > 0)
    .flatMap((d) => squareFree(toIntPoly(d)).map(({ factor, multiplicity }) => ({ factor, size: multiplicity })));
  return eigen.map(({ value, multiplicity }) => {
    const sizes = parts.filter((p) => isRootOf(value, p.factor)).map((p) => p.size).sort((a, b) => b - a);
    if (sizes.reduce((s, k) => s + k, 0) !== multiplicity) throw new Error("Jordan structure does not add up");
    return { value, multiplicity, blocks: sizes };
  });
}

const rankOf = (cols, n) => (cols.length ? run(rrefPairs(Array.from({ length: n }, (_, i) => cols.map((c) => c[i])), cols.length)).pivots.length : 0);
const matVec = (M, v) => M.map((row) => row.reduce((s, x, j) => addQ(s, mulQ(x, v[j])), Z));

/**
 * Jordan chains for a rational eigenvalue λ: returns [{ size, chain }] with
 * chain = [N^(k−1)v, …, Nv, v] (N = A − λI), the columns of P for one block.
 */
export function jordanChains(A, lambda, multiplicity) {
  const n = A.length;
  const N = A.map((row, i) => row.map((x, j) => (i === j ? subQ(x, lambda) : x)));
  const kernels = [[]];
  let Nk = N;
  for (let k = 1; ; k += 1) {
    const { rows, pivots } = run(rrefPairs(Nk, n));
    const basis = nullBasis(rows, pivots, n);
    kernels.push(basis);
    if (basis.length >= multiplicity || k > n) break;
    Nk = multiplyPairs(Nk, N);
  }
  const s = kernels.length - 1;
  const chains = [];
  for (let k = s; k >= 1; k -= 1) {
    const level = chains.flatMap((c) => (c.size > k ? [c.chain[k - 1]] : []));
    const span = [...kernels[k - 1], ...level];
    let r = rankOf(span, n);
    for (const v of kernels[k]) {
      const r2 = rankOf([...span, v], n);
      if (r2 === r) continue;
      span.push(v);
      r = r2;
      const chain = [v];
      for (let i = 1; i < k; i += 1) chain.unshift(matVec(N, chain[0]));
      chains.push({ size: k, chain });
    }
  }
  return chains.sort((a, b) => b.size - a.size);
}

/**
 * Faddeev–LeVerrier: matrices B_k with adj(xI − A) = Σ_{k=0}^{n−1} B_k x^(n−1−k).
 * The columns of adj(λI − A) are eigenvectors for a simple-enough λ.
 */
export function adjugatePolys(A) {
  const n = A.length;
  const Bs = [];
  let B = A.map((row, i) => row.map((_, j) => (i === j ? ONE : Z)));
  Bs.push(B);
  for (let k = 1; k < n; k += 1) {
    const AB = multiplyPairs(A, B);
    let tr = Z;
    for (let i = 0; i < n; i += 1) tr = addQ(tr, AB[i][i]);
    const c = reduce(-tr[0], tr[1] * BigInt(k)); // coefficient c_{n−k}
    B = AB.map((row, i) => row.map((x, j) => (i === j ? addQ(x, c) : x)));
    Bs.push(B);
  }
  // entry (i, j) as a polynomial in x, constant first
  return Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => trimQ(Bs.map((Bk) => Bk[i][j]).reverse())));
}

/**
 * An exact eigenvector whose entries are polynomials in λ (degree below
 * that of λ's polynomial), i.e. an element of Q(λ)ⁿ.
 */
export class AlgebraicVector {
  constructor(lambda, entries) {
    this.lambda = lambda;
    this.entries = entries; // pairs polys, constant first
  }

  /** Entries as strings in λ, e.g. "λ² − 3λ + 1". */
  toStrings(v = "λ") {
    return this.entries.map((p) => formatP(p, v));
  }

  /** A rational vector within about 2^-bits (entries evaluated at a rational point of λ's interval). */
  approximate(bits = 53) {
    const l = this.lambda;
    if (l instanceof ComplexAlgebraic) throw new Error("approximate: complex eigenvalue, use toNumbers()");
    l.refine(bits + 16);
    const mid = reduce(l.lo[0] * l.hi[1] + l.hi[0] * l.lo[1], 2n * l.lo[1] * l.hi[1]);
    return Matrix.fromPairs(this.entries.map((p) => [evalP(p, mid)]));
  }

  /** Floats; for complex λ, entries are { re, im }. */
  toNumbers() {
    const l = this.lambda;
    if (l instanceof ComplexAlgebraic) {
      const z = l.toNumbers();
      return this.entries.map((p) => {
        let re = 0;
        let im = 0;
        for (let i = p.length - 1; i >= 0; i -= 1) [re, im] = [re * z.re - im * z.im + pairToNumber(...p[i]), re * z.im + im * z.re];
        return { re, im };
      });
    }
    const x = l.toNumber();
    return this.entries.map((p) => p.reduceRight((s, c) => s * x + pairToNumber(...c), 0));
  }
}

/** Exact eigenvector over Q(λ), or null when adj(λI − A) vanishes (geometric multiplicity > 1). */
export function algebraicEigenvector(adj, value) {
  const g = value.poly.map((c) => [c, 1n]);
  const n = adj.length;
  for (let j = 0; j < n; j += 1) {
    const col = adj.map((row) => modP(row[j], g));
    // λ makes the column vanish iff λ is a root of gcd(entries, g)
    let G = g;
    for (const e of col) G = gcdP(G, e);
    if (degQ(G) >= 1 && isRootOf(value, toIntPoly(G))) continue;
    return new AlgebraicVector(value, col);
  }
  return null;
}

/** Everything about similarity structure, exactly. */
export function* structurePairs(A, eigen) {
  const invariant = yield* invariantFactorsPairs(A);
  const nontrivial = invariant.filter((d) => degQ(d) > 0);
  const minimal = nontrivial[nontrivial.length - 1] ?? [ONE];
  const rational = blockDiag(nontrivial.map(companionPairs));
  const all = [...eigen.real, ...eigen.complex.filter((c) => c.value)];
  const jordan = jordanStructure(invariant, all);
  let J = null;
  let P = null;
  const allRational = eigen.complexCount === 0 && eigen.real.every((e) => e.value.isRational);
  if (allRational) {
    const blocks = [];
    const cols = [];
    for (const { value, multiplicity } of eigen.real) {
      const lam = value.rational;
      for (const { size, chain } of jordanChains(A, lam, multiplicity)) {
        blocks.push(Array.from({ length: size }, (_, i) => Array.from({ length: size }, (_, j) => (i === j ? lam : j === i + 1 ? ONE : Z))));
        cols.push(...chain);
      }
    }
    J = blockDiag(blocks);
    P = Array.from({ length: A.length }, (_, i) => cols.map((c) => c[i]));
  }
  return { invariant, minimal, rational, jordan, J, P };
}

export { formatP, toRational, xP, constP, divQ };
