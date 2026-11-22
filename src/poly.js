// Exact polynomials and certified real roots come from @ratmath/algebraic
// (AlgebraicReal, Sturm isolation, square-free parts); this module re-exports
// them and adds the characteristic polynomial, which is linear algebra.

export * from "@ratmath/algebraic/poly";
import { abs, reduce, lcm } from "./q.js";

// ---------------------------------------------------------------------------
// Characteristic polynomial over Q by reduction to Hessenberg form
// (Cohen, A Course in Computational Algebraic Number Theory, 2.2.9).
// Exact, O(n³) field operations. Returns coefficient pairs, constant first,
// of det(xI − A) (monic).

import { addQ, subQ, mulQ, divQ } from "./q.js";

const Z = [0n, 1n];
const isZ = ([n]) => n === 0n;

export function* charpolyPairs(pairs) {
  const n = pairs.length;
  const H = pairs.map((r) => r.map((x) => [...x]));
  // reduce to upper Hessenberg by elementary similarity transforms
  for (let m = 1; m < n - 1; m += 1) {
    let i = m + 1;
    while (i <= n - 1 && isZ(H[i][m - 1])) i += 1;
    if (isZ(H[m][m - 1]) && i <= n - 1) {
      // swap rows i, m and columns i, m
      [H[i], H[m]] = [H[m], H[i]];
      for (const row of H) [row[i], row[m]] = [row[m], row[i]];
    }
    if (isZ(H[m][m - 1])) continue;
    const t = H[m][m - 1];
    for (let r = m + 1; r < n; r += 1) {
      if (isZ(H[r][m - 1])) continue;
      const u = divQ(H[r][m - 1], t);
      // row_r -= u·row_m ; col_m += u·col_r
      for (let j = 0; j < n; j += 1) H[r][j] = subQ(H[r][j], mulQ(u, H[m][j]));
      for (let j = 0; j < n; j += 1) H[j][m] = addQ(H[j][m], mulQ(u, H[j][r]));
    }
    yield { step: "hessenberg", done: m, of: n };
  }
  // recurrence on leading principal blocks
  const P = [[[1n, 1n]]];
  for (let m = 1; m <= n; m += 1) {
    // (x − h_mm)·P_{m−1}
    const prev = P[m - 1];
    const next = new Array(m + 1).fill(Z);
    for (let k = 0; k < prev.length; k += 1) {
      next[k + 1] = addQ(next[k + 1], prev[k]);
      next[k] = subQ(next[k], mulQ(H[m - 1][m - 1], prev[k]));
    }
    let t = [1n, 1n];
    for (let i = 1; i <= m - 1; i += 1) {
      t = mulQ(t, H[m - i][m - i - 1]);
      if (isZ(t)) break;
      const c = mulQ(t, H[m - i - 1][m - 1]);
      if (isZ(c)) continue;
      const q = P[m - i - 1];
      for (let k = 0; k < q.length; k += 1) next[k] = subQ(next[k], mulQ(c, q[k]));
    }
    P.push(next);
  }
  return P[n];
}


// ---------------------------------------------------------------------------
// Characteristic polynomial by word-size primes and Chinese remaindering.
// Same Hessenberg recurrence as above, run modulo p on doubles; much faster
// than rational arithmetic once n passes ten or so.

import { prime, invMod } from "./exact.js";

function mulmod(a, b, p, pinv) {
  const v = a * b;
  let r = v - Math.floor(v * pinv) * p;
  if (r < 0) r += p;
  else if (r >= p) r -= p;
  return r;
}

/** Coefficients of det(xI − A) mod p, constant first; A integer rows. */
export function charpolyMod(intRows, p) {
  const n = intRows.length;
  const P = BigInt(p);
  const pinv = 1 / p;
  const H = intRows.map((row) => row.map((v) => {
    const r = v === 0n ? 0 : Number(v % P);
    return r < 0 ? r + p : r;
  }));
  for (let m = 1; m < n - 1; m += 1) {
    let i = m;
    while (i < n && H[i][m - 1] === 0) i += 1;
    if (i === n) continue;
    if (i !== m) {
      [H[i], H[m]] = [H[m], H[i]];
      for (const row of H) [row[i], row[m]] = [row[m], row[i]];
    }
    const tinv = invMod(H[m][m - 1], p);
    for (let r = m + 1; r < n; r += 1) {
      if (H[r][m - 1] === 0) continue;
      const u = mulmod(H[r][m - 1], tinv, p, pinv);
      const nu = p - u;
      for (let j = 0; j < n; j += 1) H[r][j] = (H[r][j] + mulmod(nu, H[m][j], p, pinv)) % p;
      for (let j = 0; j < n; j += 1) H[j][m] = (H[j][m] + mulmod(u, H[j][r], p, pinv)) % p;
    }
  }
  const Pm = [[1]];
  for (let m = 1; m <= n; m += 1) {
    const prev = Pm[m - 1];
    const next = new Array(m + 1).fill(0);
    const h = H[m - 1][m - 1];
    for (let k = 0; k < prev.length; k += 1) {
      next[k + 1] = (next[k + 1] + prev[k]) % p;
      next[k] = (next[k] + p - mulmod(h, prev[k], p, pinv)) % p;
    }
    let t = 1;
    for (let i = 1; i <= m - 1; i += 1) {
      t = mulmod(t, H[m - i][m - i - 1], p, pinv);
      if (t === 0) break;
      const c = mulmod(t, H[m - i - 1][m - 1], p, pinv);
      if (c === 0) continue;
      const q = Pm[m - i - 1];
      for (let k = 0; k < q.length; k += 1) next[k] = (next[k] + p - mulmod(c, q[k], p, pinv)) % p;
    }
    Pm.push(next);
  }
  return Pm[n];
}

/** Exact characteristic polynomial (pairs, constant first) by multi-modular CRT. */
export function* charpolyModular(pairs) {
  const n = pairs.length;
  if (n === 0) return [[1n, 1n]];
  let D = 1n;
  for (const row of pairs) for (const [, d] of row) D = lcm(D, d);
  const ints = pairs.map((row) => row.map(([a, d]) => a * (D / d)));
  // |coefficient of x^(n−k)| ≤ C(n, k)·(√k·M)^k, M the largest entry
  let M = 1n;
  for (const row of ints) for (const v of row) if (abs(v) > M) M = abs(v);
  const lm = M.toString(2).length;
  let need = 0;
  let logBinom = 0;
  for (let k = 1; k <= n; k += 1) {
    logBinom += Math.log2((n - k + 1) / k);
    need = Math.max(need, logBinom + k * (lm + 0.5 * Math.log2(k)));
  }
  need = Math.ceil(need) + 2;
  let coeffs = new Array(n + 1).fill(0n);
  let mod = 1n;
  let bits = 0;
  for (let i = 0; bits < need; i += 1) {
    const p = prime(i);
    const r = charpolyMod(ints, p);
    const Pn = BigInt(p);
    const minv = BigInt(invMod(Number(mod % Pn), p));
    coeffs = coeffs.map((x, k) => {
      const xm = x % Pn;
      const t = (((BigInt(r[k]) - xm) % Pn + Pn) % Pn) * minv % Pn;
      return x + mod * t;
    });
    mod *= Pn;
    bits += Math.log2(p);
    yield { step: "prime", done: Math.floor(bits), of: need };
  }
  const half = mod / 2n;
  // undo the scaling: c_k(A) = c_k(DA) / D^(n−k)
  return coeffs.map((c, k) => reduce(c > half ? c - mod : c, D ** BigInt(n - k)));
}
