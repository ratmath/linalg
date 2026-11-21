// Intervals with exact rational endpoints, interval matrices, and the
// Krawczyk operator: the linear-algebra core of certified nonlinear solving.
//
// Given f : Rⁿ → Rⁿ and an enclosure J(X) of its Jacobian over a box X,
//   K(X) = m − Y·f(m) + (I − Y·J(X))·(X − m),   m = mid X,  Y ≈ J(m)⁻¹.
// If K(X) lies inside X, f has exactly one zero in X (and it is in K(X)).
// If K(X) misses X, f has no zero in X. Y only needs to be a good guess, so it
// comes from the float inverse; the certificate is exact interval arithmetic.

import { RationalInterval, Rational } from "@ratmath/core";
import { toPair, reduce, addQ, subQ, mulQ, pairToNumber, toRational } from "./q.js";
import { inverseFloat } from "./float.js";

const lt = ([a, b], [c, d]) => a * d < c * b;
const le = ([a, b], [c, d]) => a * d <= c * b;
const min = (x, y) => (lt(x, y) ? x : y);
const max = (x, y) => (lt(x, y) ? y : x);

/** Read an interval: Core RationalInterval, [lo, hi], or a single exact value. */
export function iv(x) {
  if (x && typeof x === "object" && x.low && x.high) return [toPair(x.low), toPair(x.high)];
  if (Array.isArray(x) && x.length === 2 && typeof x[0] !== "bigint") return [toPair(x[0]), toPair(x[1])];
  const p = toPair(x);
  return [p, p];
}
export const toCore = ([lo, hi]) => new RationalInterval(toRational(lo), toRational(hi));

export const ivAdd = ([a, b], [c, d]) => [addQ(a, c), addQ(b, d)];
export const ivSub = ([a, b], [c, d]) => [subQ(a, d), subQ(b, c)];
export function ivMul([a, b], [c, d]) {
  const p = [mulQ(a, c), mulQ(a, d), mulQ(b, c), mulQ(b, d)];
  return [p.reduce(min), p.reduce(max)];
}
export const ivMid = ([a, b]) => reduce(a[0] * b[1] + b[0] * a[1], 2n * a[1] * b[1]);
export const ivWidth = ([a, b]) => subQ(b, a);
export const ivInside = ([a, b], [c, d]) => lt(c, a) && lt(b, d); // [a,b] ⊂ (c,d)
export const ivIntersect = ([a, b], [c, d]) => {
  const lo = max(a, c);
  const hi = min(b, d);
  return le(lo, hi) ? [lo, hi] : null;
};

/** Round outward to the grid 2^-bits (keeps endpoint sizes bounded). */
export function outward([lo, hi], bits) {
  const s = 2n ** BigInt(bits);
  const fl = (n, d) => (n >= 0n ? n / d : -((-n + d - 1n) / d));
  const l = fl(lo[0] * s, lo[1]);
  const h = -fl(-hi[0] * s, hi[1]);
  return [reduce(l, s), reduce(h, s)];
}

/** Round a rational to the grid 2^-bits (nearest below). */
function roundPair(p, bits) {
  return outward([p, p], bits)[0];
}

/**
 * One Krawczyk test-and-contract on box X (array of intervals).
 * f(point) gets Core Rationals and returns exact values or enclosures;
 * jacobian(box) gets Core RationalIntervals and returns an interval matrix.
 * Returns { status: "unique" | "none" | "unknown", box }.
 */
export function krawczykStep(f, jacobian, X, bits = 64) {
  const n = X.length;
  const m = X.map((x) => roundPair(ivMid(x), bits));
  const J = jacobian(X.map(toCore)).map((row) => row.map(iv));
  const Jmid = new Float64Array(n * n);
  J.forEach((row, i) => row.forEach((x, j) => (Jmid[i * n + j] = pairToNumber(...ivMid(x)))));
  const inv = inverseFloat(Jmid, n).value;
  if (!inv) return { status: "unknown", box: X };
  const Y = inv.map((row) => row.map((y) => toPair(y)));
  const fm = f(m.map(toRational)).map(iv);
  const K = [];
  for (let i = 0; i < n; i += 1) {
    let k = [m[i], m[i]];
    for (let j = 0; j < n; j += 1) k = ivSub(k, ivMul([Y[i][j], Y[i][j]], fm[j]));
    for (let j = 0; j < n; j += 1) {
      // (I − Y·J)_ij
      let c = i === j ? [[1n, 1n], [1n, 1n]] : [[0n, 1n], [0n, 1n]];
      for (let t = 0; t < n; t += 1) c = ivSub(c, ivMul([Y[i][t], Y[i][t]], J[t][j]));
      k = ivAdd(k, ivMul(c, ivSub(X[j], [m[j], m[j]])));
    }
    K.push(outward(k, bits));
  }
  if (K.every((k, i) => ivInside(k, X[i]))) return { status: "unique", box: K };
  const cut = K.map((k, i) => ivIntersect(k, X[i]));
  if (cut.some((c) => c === null)) return { status: "none", box: null };
  return { status: "unknown", box: cut };
}

const widthBits = (box) => Math.min(...box.map((x) => -Math.log2(pairToNumber(...ivWidth(x)) || 2 ** -1e4)));

/**
 * Iterate Krawczyk on X until the box is certified and narrower than 2^-bits,
 * shown to hold no zero, or stuck. A generator: one yield per step.
 */
export function* krawczyk(f, jacobian, X, { bits = 53, maxSteps = 200 } = {}) {
  let status = "unknown";
  let box = X;
  for (let step = 0; step < maxSteps; step += 1) {
    const r = krawczykStep(f, jacobian, box, bits + 16);
    if (r.status === "none") return { status: "none", box: null };
    if (r.status === "unique") status = "unique";
    const before = widthBits(box);
    box = status === "unique" ? r.box.map((k, i) => ivIntersect(k, box[i]) ?? k) : r.box;
    yield { step: "krawczyk", done: step + 1, of: maxSteps };
    if (status === "unique" && widthBits(box) >= bits) break;
    if (status !== "unique" && widthBits(box) <= before + 1e-9) break; // not contracting, not certified
    if (status === "unique" && widthBits(box) <= before + 1e-9) break; // as tight as the arithmetic allows
  }
  return { status, box };
}

export { Rational };
