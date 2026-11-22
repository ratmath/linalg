// Exact spectral data: characteristic polynomial, eigenvalues as certified
// real algebraic numbers, singular values likewise.

import { charpolyModular, fromPairs, squareFree, isolate, degree, countRealRoots, AlgebraicReal, signAt } from "./poly.js";
import { sqrtFloor } from "./dixon.js";
import { reduce } from "./q.js";
import { certifyComplexRoots } from "./complex.js";
import { multiplyPairs, transposePairs } from "./decomp.js";

/**
 * Eigen-structure over the reals: the characteristic polynomial (pairs,
 * constant first), and for each square-free factor its real roots (ascending)
 * and how many complex (non-real) roots it has.
 */
export function* eigenExact(pairs) {
  const cp = yield* charpolyModular(pairs);
  const factors = [];
  for (const { factor, multiplicity } of squareFree(fromPairs(cp))) {
    const roots = isolate(factor);
    factors.push({ factor, multiplicity, roots, complexCount: degree(factor) - roots.length });
    yield { step: "factor", done: factors.length, of: factors.length };
  }
  const real = [];
  const complex = [];
  let complexCertified = true;
  for (const f of factors) {
    for (const r of f.roots) real.push({ value: r, multiplicity: f.multiplicity });
    const zs = certifyComplexRoots(f.factor, f.complexCount);
    if (zs === null) complexCertified = false;
    else for (const z of zs) complex.push({ value: z, multiplicity: f.multiplicity });
  }
  real.sort((a, b) => cmpReal(a.value, b.value));
  return { charpoly: cp, factors, real, complex, complexCertified, complexCount: factors.reduce((s, f) => s + f.complexCount * f.multiplicity, 0) };
}

function cmpReal(a, b) {
  // isolating intervals of distinct roots can overlap only across factors; refine until apart
  for (let bits = 8; bits < 4096; bits *= 2) {
    const [al, ah] = [a.refine(bits).lo, a.hi];
    const [bl, bh] = [b.refine(bits).lo, b.hi];
    if (lt(ah, bl) || (eq(ah, bl) && !(eq(al, ah) && eq(bl, bh)))) return -1;
    if (lt(bh, al) || (eq(bh, al) && !(eq(al, ah) && eq(bl, bh)))) return 1;
    if (eq(al, ah) && eq(bl, bh) && eq(al, bl)) return 0;
  }
  return 0;
}
const lt = ([a, b], [c, d]) => a * d < c * b;
const eq = ([a, b], [c, d]) => a * d === c * b;

/** Singular values (descending) as AlgebraicReals with multiplicities; zeros included up to min(m, n). */
export function* singularExact(pairs) {
  const m = pairs.length;
  const n = pairs[0]?.length ?? 0;
  const At = transposePairs(pairs);
  const G = m <= n ? multiplyPairs(pairs, At) : multiplyPairs(At, pairs); // k×k Gram matrix
  const cp = yield* charpolyModular(G);
  // isolate the eigenvalues y of the Gram matrix (degree k, cheap), then map
  // each interval through √ into an isolating interval of a root of f(x²)
  const ys = [];
  for (let { factor, multiplicity } of squareFree(fromPairs(cp))) {
    for (const r of isolate(factor)) ys.push({ y: r, factor, multiplicity });
  }
  ys.sort((a, b) => cmpReal(a.y, b.y));
  // separate neighbouring intervals so the gaps between them are honest
  for (let bits = 16; ; bits *= 2) {
    ys.forEach((e) => e.y.refine(bits));
    if (ys.every((e, i) => i === 0 || lt(ys[i - 1].y.hi, e.y.lo) || (ys[i - 1].y.isRational && e.y.isRational))) break;
    if (bits > 1 << 12) break;
  }
  const out = [];
  ys.forEach(({ y, factor, multiplicity }, i) => {
    if (y.isRational) {
      const [a, b] = y.rational;
      if (a === 0n) return out.push({ value: AlgebraicReal.rational([0n, 1n]), multiplicity });
      const ra = sqrtFloor(a * b);
      if (ra * ra === a * b) return out.push({ value: AlgebraicReal.rational(reduce(ra, b)), multiplicity });
    }
    const h = [];
    factor.forEach((c, k) => {
      h[2 * k] = c;
      if (k < factor.length - 1) h[2 * k + 1] = 0n;
    });
    const lo2 = y.isRational ? y.rational : y.lo; // need lo² ≤ lo2, lo² ≥ previous hi
    const hi2 = y.isRational ? y.rational : y.hi;
    const prev = i > 0 ? (ys[i - 1].y.isRational ? ys[i - 1].y.rational : ys[i - 1].y.hi) : [0n, 1n];
    const next = i + 1 < ys.length ? (ys[i + 1].y.isRational ? ys[i + 1].y.rational : ys[i + 1].y.lo) : null;
    for (let k = 8n; ; k *= 2n) {
      const S = 4n ** k;
      const lo = reduce(sqrtFloor((lo2[0] * S) / lo2[1]), 2n ** k);
      let hn = sqrtFloor((hi2[0] * S + hi2[1] - 1n) / hi2[1]);
      if (hn * hn * hi2[1] < hi2[0] * S) hn += 1n;
      const hi = reduce(hn, 2n ** k);
      const sq = ([a, b]) => [a * a, b * b];
      const okLo = !lt(sq(lo), prev) && (lo[0] > 0n || prev[0] === 0n);
      const okHi = next === null || !lt(next, sq(hi));
      const strict = !eq(sq(lo), lo2) || !y.isRational;
      if (okLo && okHi && strict && signAt(h, lo) !== 0 && signAt(h, hi) !== 0) {
        out.push({ value: new AlgebraicReal(h, lo, hi), multiplicity });
        return;
      }
      if (k > 4096n) throw new Error("singular value isolation failed");
    }
  });
  out.sort((a, b) => -cmpReal(a.value, b.value));
  return out;
}

export { countRealRoots };
