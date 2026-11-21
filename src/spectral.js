// Exact spectral data: characteristic polynomial, eigenvalues as certified
// real algebraic numbers, singular values likewise.

import { charpolyModular, fromPairs, squareFree, isolate, degree, countRealRoots } from "./poly.js";
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
  for (const f of factors) for (const r of f.roots) real.push({ value: r, multiplicity: f.multiplicity });
  real.sort((a, b) => cmpReal(a.value, b.value));
  return { charpoly: cp, factors, real, complexCount: factors.reduce((s, f) => s + f.complexCount * f.multiplicity, 0) };
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
  const out = [];
  for (let { factor, multiplicity } of squareFree(fromPairs(cp))) {
    if (factor[0] === 0n) {
      // y = 0 is a root: a zero singular value; the rest of the factor goes on
      out.push({ value: isolate([0n, 1n])[0], multiplicity });
      factor = factor.slice(1);
      if (factor.length < 2) continue;
    }
    // roots y ≥ 0 of factor(y); σ = √y is a positive root of factor(x²)
    const h = [];
    factor.forEach((c, i) => {
      h[2 * i] = c;
      if (i < factor.length - 1) h[2 * i + 1] = 0n;
    });
    for (const r of isolate(h)) if (r.lo[0] >= 0n && (r.lo[0] > 0n || r.hi[0] > 0n)) out.push({ value: r, multiplicity });
  }
  out.sort((a, b) => -cmpReal(a.value, b.value));
  return out;
}

export { countRealRoots };
