// p-adic (Dixon) lifting for exact solutions of nonsingular systems.
//
// Invert A once modulo a word-size prime p, then peel off the solution one
// base-p digit at a time: x_i = A⁻¹·r (mod p), r ← (r − A·x_i)/p. After k
// digits X = Σ x_i p^i satisfies A·X ≡ b (mod p^k). Rational reconstruction
// turns X into fractions once p^k is large enough, and an exact check A·x = b
// decides when to stop. All O(n³) work is modulo p on doubles; the BigInt work
// per digit is O(n²) with small multipliers.

import { prime, integerRows, invMod, hadamardBits } from "./exact.js";
import { reduce, abs, bitLength } from "./q.js";

/** A⁻¹ mod p as rows of numbers, or null if A is singular mod p. */
function inverseMod(rows, p) {
  const n = rows.length;
  const P = BigInt(p);
  const pinv = 1 / p;
  const mod = (v) => {
    let r = v - Math.floor(v * pinv) * p;
    if (r < 0) r += p;
    else if (r >= p) r -= p;
    return r;
  };
  const a = rows.map((row, i) => {
    const r = new Float64Array(2 * n);
    row.forEach((v, j) => {
      const x = v === 0n ? 0 : Number(v % P);
      r[j] = x < 0 ? x + p : x;
    });
    r[n + i] = 1;
    return r;
  });
  for (let col = 0; col < n; col += 1) {
    let k = col;
    while (k < n && a[k][col] === 0) k += 1;
    if (k === n) return null;
    [a[k], a[col]] = [a[col], a[k]];
    const inv = invMod(a[col][col], p);
    const prow = a[col];
    for (let j = 0; j < 2 * n; j += 1) prow[j] = mod(prow[j] * inv);
    for (let i = 0; i < n; i += 1) {
      if (i === col || a[i][col] === 0) continue;
      const f = p - a[i][col];
      const row = a[i];
      for (let j = 0; j < 2 * n; j += 1) if (prow[j]) row[j] = mod(row[j] + f * prow[j]);
    }
  }
  return a.map((r) => Array.from(r.subarray(n)));
}

/** Rational reconstruction: a/b ≡ u (mod m) with |a|, b ≤ √(m/2), or null. */
export function ratRecon(u, m) {
  const bound = sqrtFloor(m / 2n);
  let [r0, r1] = [m, ((u % m) + m) % m];
  let [t0, t1] = [0n, 1n];
  while (r1 > bound) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [t0, t1] = [t1, t0 - q * t1];
  }
  if (t1 === 0n || abs(t1) > bound) return null;
  return reduce(t1 < 0n ? -r1 : r1, abs(t1));
}

export function sqrtFloor(n) {
  if (n < 2n) return n;
  let x = 1n << BigInt(Math.ceil(bitLength(n) / 2));
  for (;;) {
    const y = (x + n / x) >> 1n;
    if (y >= x) return x;
    x = y;
  }
}

/**
 * Solve A·x = B exactly for square nonsingular A (pairs). Returns rows of
 * pairs, or null if A looks singular (the caller then uses elimination).
 */
export function* dixonSolve(Apairs, Bpairs) {
  const n = Apairs.length;
  const k = Bpairs[0]?.length ?? 0;
  // clear denominators row by row on [A | B] together
  const { rows } = integerRows(Apairs.map((r, i) => [...r, ...Bpairs[i]]));
  const A = rows.map((r) => r.slice(0, n));
  const B = rows.map((r) => r.slice(n));
  let p = 0;
  let C = null;
  for (let i = 40; i < 46 && !C; i += 1) {
    p = prime(i);
    C = inverseMod(A, p);
  }
  if (!C) return null;
  const P = BigInt(p);
  const pinv = 1 / p;
  // stop for sure once p^k exceeds 2·N·D (Cramer: D ≤ H(A), N ≤ H([A|b]))
  let bMax = 1n;
  for (const r of B) for (const v of r) if (abs(v) > bMax) bMax = abs(v);
  const capBits = 2 * (hadamardBits(A) + bitLength(bMax) + Math.ceil(Math.log2(n + 1)) + 1) + 4;
  const out = Array.from({ length: n }, () => new Array(k));
  for (let c = 0; c < k; c += 1) {
    let r = B.map((row) => row[c]);
    const X = new Array(n).fill(0n);
    let pk = 1n;
    let bits = 0;
    let nextCheck = 8;
    for (let step = 1; ; step += 1) {
      const rm = r.map((v) => {
        const x = v === 0n ? 0 : Number(v % P);
        return x < 0 ? x + p : x;
      });
      const xi = C.map((row) => {
        let s = 0;
        for (let j = 0; j < n; j += 1) {
          if (!rm[j]) continue;
          // s + row·rm < 2^53; reduce with the reciprocal trick
          const v = s + row[j] * rm[j];
          s = v - Math.floor(v * pinv) * p;
          if (s < 0) s += p;
          else if (s >= p) s -= p;
        }
        return s;
      });
      for (let i = 0; i < n; i += 1) X[i] += pk * BigInt(xi[i]);
      const xb = xi.map(BigInt);
      r = r.map((v, i) => {
        let s = v;
        const Ai = A[i];
        for (let j = 0; j < n; j += 1) if (xb[j]) s -= Ai[j] * xb[j];
        return s / P;
      });
      pk *= P;
      bits += Math.log2(p);
      if (step % 4 === 0) yield { step: "lift", done: Math.min(Math.floor(bits), capBits), of: capBits };
      if (step >= nextCheck || bits >= capBits) {
        nextCheck = Math.ceil(step * 1.25) + 1;
        const x = reconstruct(X, pk);
        if (x && checks(A, B, c, x)) {
          x.forEach((v, i) => (out[i][c] = v));
          break;
        }
        if (bits >= capBits) return null;
      }
    }
  }
  return out;
}

function reconstruct(X, m) {
  // common-denominator shortcut: most entries share the denominator found first
  let D = 1n;
  const half = m / 2n;
  const bound = sqrtFloor(half);
  const out = [];
  for (const u of X) {
    let y = (u * D) % m;
    if (y > half) y -= m;
    if (abs(y) <= bound) {
      out.push(reduce(y, D));
      continue;
    }
    const q = ratRecon(u, m);
    if (!q) return null;
    out.push(q);
    D = D * q[1] / gcd(D, q[1]);
  }
  return out;
}

const gcd = (a, b) => {
  while (b) [a, b] = [b, a % b];
  return a < 0n ? -a : a;
};

function checks(A, B, c, x) {
  let L = 1n;
  for (const [, d] of x) L = L * d / gcd(L, d);
  const xi = x.map(([a, d]) => a * (L / d));
  return A.every((row, i) => row.reduce((s, v, j) => s + v * xi[j], 0n) === B[i][c] * L);
}
