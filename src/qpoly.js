// Polynomials with rational coefficients ([n, d] pairs, constant first), for
// the places that need field division: Smith form over Q[x], minimal and
// invariant factor polynomials, eigenvectors over Q(λ).

import { addQ, subQ, mulQ, divQ, reduce, lcm } from "./q.js";

const Z = [0n, 1n];
export const isZeroQ = ([n]) => n === 0n;

export function trimQ(p) {
  let n = p.length;
  while (n > 0 && isZeroQ(p[n - 1])) n -= 1;
  return p.slice(0, n);
}
export const degQ = (p) => trimQ(p).length - 1;
export const isZeroPoly = (p) => degQ(p) < 0;

export function addP(a, b) {
  const n = Math.max(a.length, b.length);
  return trimQ(Array.from({ length: n }, (_, i) => addQ(a[i] ?? Z, b[i] ?? Z)));
}
export function subP(a, b) {
  const n = Math.max(a.length, b.length);
  return trimQ(Array.from({ length: n }, (_, i) => subQ(a[i] ?? Z, b[i] ?? Z)));
}
export function mulP(a, b) {
  if (!a.length || !b.length) return [];
  const o = new Array(a.length + b.length - 1).fill(Z);
  a.forEach((x, i) => {
    if (isZeroQ(x)) return;
    b.forEach((y, j) => (o[i + j] = addQ(o[i + j], mulQ(x, y))));
  });
  return trimQ(o);
}
export function divmodP(a, b) {
  a = trimQ(a);
  b = trimQ(b);
  if (!b.length) throw new RangeError("division by the zero polynomial");
  const db = b.length - 1;
  const q = new Array(Math.max(0, a.length - db)).fill(Z);
  let r = [...a];
  while (r.length - 1 >= db && r.length) {
    const k = r.length - 1 - db;
    const c = divQ(r[r.length - 1], b[db]);
    q[k] = c;
    for (let i = 0; i <= db; i += 1) r[k + i] = subQ(r[k + i], mulQ(c, b[i]));
    r = trimQ(r);
  }
  return [trimQ(q), r];
}
export const modP = (a, b) => divmodP(a, b)[1];
export function monicP(p) {
  p = trimQ(p);
  if (!p.length) return p;
  const lc = p[p.length - 1];
  return p.map((c) => divQ(c, lc));
}
export function gcdP(a, b) {
  a = trimQ(a);
  b = trimQ(b);
  while (b.length) [a, b] = [b, modP(a, b)];
  return monicP(a);
}
export const scaleP = (p, s) => trimQ(p.map((c) => mulQ(c, s)));
export const constP = (c) => trimQ([c]);
export const xP = () => [Z, [1n, 1n]];

/** Integer primitive polynomial (BigInt coefficients) with the same roots. */
export function toIntPoly(p) {
  p = trimQ(p);
  let L = 1n;
  for (const [, d] of p) L = lcm(L, d);
  let ints = p.map(([n, d]) => n * (L / d));
  let g = 0n;
  for (const c of ints) {
    let a = c < 0n ? -c : c;
    let b = g;
    while (b) [a, b] = [b, a % b];
    g = a;
  }
  if (ints[ints.length - 1] < 0n) g = -g;
  return g ? ints.map((c) => c / g) : ints;
}
export const fromIntPoly = (p) => p.map((c) => [c, 1n]);

/** p evaluated at the exact rational x. */
export function evalP(p, x) {
  let acc = Z;
  for (let i = p.length - 1; i >= 0; i -= 1) acc = addQ(mulQ(acc, x), p[i]);
  return acc;
}

/** Human form, highest degree first, in the variable `v`. */
export function formatP(p, v = "x") {
  p = trimQ(p);
  if (!p.length) return "0";
  const parts = [];
  for (let i = p.length - 1; i >= 0; i -= 1) {
    const [n, d] = p[i];
    if (n === 0n) continue;
    const neg = n < 0n;
    const a = neg ? -n : n;
    const coef = d === 1n ? `${a}` : `${a}/${d}`;
    const mono = i === 0 ? coef : `${a === 1n && d === 1n ? "" : coef}${v}${i > 1 ? toSup(i) : ""}`;
    parts.push({ neg, mono });
  }
  return parts.map((t, k) => (k === 0 ? (t.neg ? "−" : "") + t.mono : `${t.neg ? " − " : " + "}${t.mono}`)).join("");
}
const SUP = "⁰¹²³⁴⁵⁶⁷⁸⁹";
const toSup = (k) => String(k).split("").map((c) => SUP[c]).join("");

export { reduce };
