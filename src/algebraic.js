// Exact arithmetic on real algebraic numbers (AlgebraicReal): each result is
// again a square-free integer polynomial with an isolating interval.
//
// a ± b, a·b come from resultants (the polynomial whose roots are all sums or
// products of a root of p with a root of q); the right root is then picked by
// refining interval enclosures until exactly one candidate fits. 1/a reverses
// the coefficients; √a substitutes x². Polynomials are not factored, so
// degrees multiply with each operation: fine for a few steps (√2 + √3 has
// degree 4), and guarded by MAX_DEGREE.

import { AlgebraicReal, squareFree, isolate, primitive, degree, signAt } from "./poly.js";
import { reduce, addQ, subQ, mulQ, divQ } from "./q.js";
import { run, determinantBareiss } from "./exact.js";
import { sqrtFloor } from "./dixon.js";

export const MAX_DEGREE = 48;

const Q = (n, d = 1n) => reduce(BigInt(n), BigInt(d));
const cmp = ([a, b], [c, d]) => {
  const x = a * d - c * b;
  return x > 0n ? 1 : x < 0n ? -1 : 0;
};

/** A rational pair [n, d] as an AlgebraicReal. */
export const fromRational = (pair) => AlgebraicReal.rational(reduce(pair[0], pair[1]));

/** The polynomial to compute with: linear when rational; no stray root at 0. */
function usable(a) {
  const r = a.rational;
  if (r) return [-r[0], r[1]];
  let p = a.poly;
  while (p.length > 1 && p[0] === 0n) p = p.slice(1); // a ≠ 0, so the root 0 is someone else's
  return p;
}

// ---------------------------------------------------------------------------
// Resultants by evaluation and interpolation

/** det of the Sylvester matrix of p and q (integer coefficients, constant first). */
function resultant(p, q) {
  const m = p.length - 1;
  const n = q.length - 1;
  const size = m + n;
  const rows = [];
  const hi = (c) => [...c].reverse(); // highest first
  for (let i = 0; i < n; i += 1) rows.push(Array.from({ length: size }, (_, j) => [j - i >= 0 && j - i <= m ? hi(p)[j - i] : 0n, 1n]));
  for (let i = 0; i < m; i += 1) rows.push(Array.from({ length: size }, (_, j) => [j - i >= 0 && j - i <= n ? hi(q)[j - i] : 0n, 1n]));
  return run(determinantBareiss(rows))[0];
}

const binom = (n, k) => {
  let r = 1n;
  for (let i = 0n; i < BigInt(k); i += 1n) r = (r * (BigInt(n) - i)) / (i + 1n);
  return r;
};

/** Coefficients in y of q(t − y). */
function shifted(q, t) {
  const out = Array(q.length).fill(0n);
  q.forEach((c, k) => {
    for (let j = 0; j <= k; j += 1) out[j] += c * binom(k, j) * t ** BigInt(k - j) * (j % 2 ? -1n : 1n);
  });
  return out;
}

/** Coefficients in y of yⁿ·q(t / y). */
function scaled(q, t) {
  const n = q.length - 1;
  return q.map((_, j) => q[n - j] * t ** BigInt(n - j));
}

/** Polynomial of degree ≤ D through (t, f(t)), t = 0..D, as primitive integer coefficients. */
function interpolate(D, f) {
  const xs = Array.from({ length: D + 1 }, (_, t) => BigInt(t));
  const coef = xs.map((t) => Q(f(t)));
  // Newton divided differences
  for (let j = 1; j <= D; j += 1) for (let i = D; i >= j; i -= 1) coef[i] = divQ(subQ(coef[i], coef[i - 1]), Q(xs[i] - xs[i - j]));
  let poly = [coef[D]];
  for (let i = D - 1; i >= 0; i -= 1) {
    const next = Array(poly.length + 1).fill(Q(0));
    poly.forEach((c, k) => {
      next[k + 1] = addQ(next[k + 1], c);
      next[k] = subQ(next[k], mulQ(c, Q(xs[i])));
    });
    next[0] = addQ(next[0], coef[i]);
    poly = next;
  }
  let L = 1n;
  const g = (a, b) => (b ? g(b, a % b) : a);
  for (const [, d] of poly) L = (L * d) / g(L, d);
  const ints = poly.map(([n, d]) => n * (L / d));
  while (ints.length > 1 && ints[ints.length - 1] === 0n) ints.pop();
  return primitive(ints);
}

// ---------------------------------------------------------------------------
// Picking the right root

/** Interval [lo, hi] of a, width below 2^-bits, as pairs. */
function box(a, bits) {
  const r = a.rational;
  if (r) return [r, r];
  a.refine(bits);
  return [a.lo, a.hi];
}

function pick(poly, enclose) {
  if (degree(poly) > MAX_DEGREE) throw new RangeError(`the result would need a polynomial of degree ${degree(poly)}; the limit is ${MAX_DEGREE}`);
  const candidates = squareFree(poly).flatMap(({ factor }) => isolate(factor));
  for (let bits = 8; bits <= 4096; bits *= 2) {
    const [lo, hi] = enclose(bits);
    const hits = candidates.filter((c) => {
      const [cl, ch] = box(c, bits);
      return cmp(ch, lo) >= 0 && cmp(cl, hi) <= 0;
    });
    if (hits.length === 1) return hits[0];
    if (hits.length === 0) throw new Error("algebraic: lost the root (enclosure missed every candidate)");
  }
  throw new Error("algebraic: could not separate the roots");
}

// ---------------------------------------------------------------------------
// Operations

export function neg(a) {
  const r = a.rational;
  if (r) return fromRational([-r[0], r[1]]);
  const p = a.poly.map((c, i) => (i % 2 ? -c : c));
  return new AlgebraicReal(primitive(p), [-a.hi[0], a.hi[1]], [-a.lo[0], a.lo[1]]);
}

export function add(a, b) {
  const ra = a.rational;
  const rb = b.rational;
  if (ra && rb) return fromRational(addQ(ra, rb));
  const p = usable(a);
  const q = usable(b);
  const D = (p.length - 1) * (q.length - 1);
  const R = interpolate(D, (t) => resultant(p, shifted(q, t)));
  return pick(R, (bits) => {
    const [al, ah] = box(a, bits);
    const [bl, bh] = box(b, bits);
    return [addQ(al, bl), addQ(ah, bh)];
  });
}

export const sub = (a, b) => add(a, neg(b));

export function mul(a, b) {
  const ra = a.rational;
  const rb = b.rational;
  if (ra && rb) return fromRational(mulQ(ra, rb));
  if ((ra && ra[0] === 0n) || (rb && rb[0] === 0n)) return fromRational([0n, 1n]);
  const p = usable(a);
  const q = usable(b);
  const D = (p.length - 1) * (q.length - 1);
  const R = interpolate(D, (t) => resultant(p, scaled(q, t)));
  return pick(R, (bits) => {
    const [al, ah] = box(a, bits);
    const [bl, bh] = box(b, bits);
    const prods = [mulQ(al, bl), mulQ(al, bh), mulQ(ah, bl), mulQ(ah, bh)];
    return [prods.reduce((x, y) => (cmp(x, y) <= 0 ? x : y)), prods.reduce((x, y) => (cmp(x, y) >= 0 ? x : y))];
  });
}

/** Sign of a: -1, 0 or 1, exactly. */
export function sign(a) {
  const r = a.rational;
  if (r) return r[0] > 0n ? 1 : r[0] < 0n ? -1 : 0;
  // an irrational a is not 0; refine until 0 leaves the interval
  for (let bits = 8; ; bits *= 2) {
    const [lo, hi] = box(a, bits);
    if (lo[0] > 0n) return 1;
    if (hi[0] < 0n) return -1;
  }
}

export function inv(a) {
  const r = a.rational;
  if (r) {
    if (r[0] === 0n) throw new RangeError("division by zero");
    return fromRational([r[1], r[0]]);
  }
  sign(a); // moves the interval off 0
  const p = usable(a);
  const rev = primitive([...p].reverse());
  return new AlgebraicReal(rev[rev.length - 1] < 0n ? rev.map((c) => -c) : rev, ...[reduce(a.hi[1], a.hi[0]), reduce(a.lo[1], a.lo[0])].sort(cmp));
}

export const div = (a, b) => mul(a, inv(b));

export function pow(a, k) {
  if (k < 0) return pow(inv(a), -k);
  let out = fromRational([1n, 1n]);
  let base = a;
  while (k > 0) {
    if (k & 1) out = mul(out, base);
    k >>= 1;
    if (k) base = mul(base, base);
  }
  return out;
}

/** -1, 0, 1 as a compares to b. */
export const compare = (a, b) => sign(sub(a, b));

/** Rational bounds on √x for x ≥ 0, about `bits` bits. */
function sqrtBounds([n, d], bits) {
  const s = 2n ** BigInt(bits);
  const lo = sqrtFloor((n * d * s * s));
  return [reduce(lo, d * s), reduce(lo + 1n, d * s)];
}

export function sqrt(a) {
  const sg = sign(a);
  if (sg < 0) throw new RangeError("√ of a negative number is not real (turn on complex numbers)");
  if (sg === 0) return a;
  const r = a.rational;
  if (r) {
    const sn = sqrtFloor(r[0]);
    const sd = sqrtFloor(r[1]);
    if (sn * sn === r[0] && sd * sd === r[1]) return fromRational([sn, sd]);
  }
  const p = usable(a);
  const p2 = [];
  p.forEach((c, i) => {
    p2[2 * i] = c;
    if (i < p.length - 1) p2[2 * i + 1] = 0n;
  });
  return pick(primitive(p2), (bits) => {
    const [lo, hi] = box(a, bits + 4);
    return [sqrtBounds(lo, bits)[0], sqrtBounds(hi, bits)[1]];
  });
}

/** The real roots of an integer or rational polynomial (constant first), ascending, with multiplicities. */
export function realRoots(coeffPairs) {
  let L = 1n;
  const g = (x, y) => (y ? g(y, x % y) : x);
  for (const [, d] of coeffPairs) L = (L * d) / g(L, d);
  const ints = coeffPairs.map(([n, d]) => n * (L / d));
  while (ints.length > 1 && ints[ints.length - 1] === 0n) ints.pop();
  if (ints.length < 2) return [];
  const out = [];
  for (const { factor, multiplicity } of squareFree(ints)) for (const value of isolate(factor)) out.push({ value, multiplicity });
  return out.sort((x, y) => compare(x.value, y.value));
}

export { signAt };
