// The matrix exponential, certified. exp(A) is irrational in general, so the
// exact stage returns an interval for every entry, guaranteed to contain the
// true value and narrower than asked.
//
// Scaling and squaring: B = A / 2^s with ‖B‖∞ ≤ 1/2, Taylor polynomial of
// degree K for exp(B) plus a rigorous bound on the tail added to every entry,
// then s interval squarings. Endpoints are rounded outward to a dyadic grid so
// their size stays bounded.

import { reduce, addQ, mulQ, abs, pairToNumber } from "./q.js";
import { ivAdd, ivMul, outward } from "./interval.js";

const Z = [0n, 1n];

function normInfPairs(A) {
  // an upper bound for max_i Σ_j |a_ij| as a pair
  let best = Z;
  for (const row of A) {
    let s = Z;
    for (const [n, d] of row) s = addQ(s, [abs(n), d]);
    if (s[0] * best[1] > best[0] * s[1]) best = s;
  }
  return best;
}

function ivMatMul(X, Y, bits) {
  const n = X.length;
  return X.map((row, i) => Y[0].map((_, j) => {
    let acc = [Z, Z];
    for (let k = 0; k < n; k += 1) acc = ivAdd(acc, ivMul(X[i][k], Y[k][j]));
    return outward(acc, bits);
  }));
}

/** Generator: interval matrix enclosing exp(A), every entry narrower than 2^-bits when it can be. */
export function* expmPairs(A, bits = 53) {
  const n = A.length;
  const norm = normInfPairs(A);
  const normF = pairToNumber(...norm);
  let s = 0;
  while (normF / 2 ** s > 0.5) s += 1;
  // the answer can be as large as e^‖A‖; carry enough extra bits for that and for s squarings
  let work = bits + 8 + s + Math.ceil(normF * 1.45) + Math.ceil(Math.log2(n + 1)) * (s + 1);
  for (let attempt = 0; attempt < 4; attempt += 1) {
    const scale = 2n ** BigInt(s);
    const B = A.map((row) => row.map(([a, d]) => reduce(a, d * scale)));
    const Bnorm = reduce(norm[0], norm[1] * scale); // ≤ 1/2
    // Taylor: Σ_{k≤K} B^k / k!, tail ≤ ‖B‖^(K+1)/(K+1)! · 2 (since ‖B‖ ≤ 1/2)
    let K = 1;
    let tail = 2 * pairToNumber(...Bnorm) ** 2 / 2;
    while (tail > 2 ** -(work + 2) && K < 400) {
      K += 1;
      tail *= pairToNumber(...Bnorm) / (K + 1);
    }
    const I = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? [[1n, 1n], [1n, 1n]] : [Z, Z])));
    const Bi = B.map((row) => row.map((x) => [x, x]));
    let term = I;
    let sum = I;
    for (let k = 1; k <= K; k += 1) {
      term = ivMatMul(term, Bi, work).map((row) => row.map(([lo, hi]) => outward([reduce(lo[0], lo[1] * BigInt(k)), reduce(hi[0], hi[1] * BigInt(k))], work)));
      sum = sum.map((row, i) => row.map((x, j) => ivAdd(x, term[i][j])));
      if (k % 4 === 0) yield { step: "taylor", done: k, of: K + s };
    }
    // rigorous tail: 2·‖B‖^(K+1)/(K+1)!, computed with exact rationals then rounded up
    let t = [2n, 1n];
    for (let k = 1; k <= K + 1; k += 1) t = mulQ(t, reduce(Bnorm[0], Bnorm[1] * BigInt(k)));
    const T = outward([[-t[0], t[1]], t], work);
    let E = sum.map((row) => row.map((x) => outward(ivAdd(x, T), work)));
    for (let k = 0; k < s; k += 1) {
      E = ivMatMul(E, E, work);
      yield { step: "square", done: K + k + 1, of: K + s };
    }
    const widest = Math.max(...E.flat().map(([lo, hi]) => pairToNumber(hi[0] * lo[1] - lo[0] * hi[1], hi[1] * lo[1])));
    if (widest <= 2 ** -bits) return E;
    work += Math.ceil(Math.log2(widest / 2 ** -bits)) + 8;
  }
  throw new Error("expm: could not reach the requested width");
}

/** Float preview: scaling and squaring with an 18-term Taylor polynomial. */
export function expmFloat(rows) {
  const n = rows.length;
  const norm = Math.max(0, ...rows.map((r) => r.reduce((s, x) => s + Math.abs(x), 0)));
  let s = 0;
  while (norm / 2 ** s > 0.5) s += 1;
  const B = rows.map((r) => r.map((x) => x / 2 ** s));
  const mul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((acc, x, k) => acc + x * Y[k][j], 0)));
  let term = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  let sum = term.map((r) => [...r]);
  for (let k = 1; k <= 18; k += 1) {
    term = mul(term, B).map((r) => r.map((x) => x / k));
    sum = sum.map((r, i) => r.map((x, j) => x + term[i][j]));
  }
  for (let k = 0; k < s; k += 1) sum = mul(sum, sum);
  return sum;
}
