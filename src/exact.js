// Exact algorithms. Every routine is a generator: it yields a small progress
// record after each unit of work (a pivot, a prime) so a scheduler can stop,
// let the page breathe, and resume, or drop the work when a slider moves on.
// The answer is the generator's return value; `run(gen)` drives one to the end.

import { gcd, lcm, reduce, abs, bitLength } from "./q.js";

export function run(gen) {
  for (;;) {
    const { value, done } = gen.next();
    if (done) return value;
  }
}

/**
 * Clear denominators row by row. Returns integer rows and each row's scale
 * factor L_i (row_i was multiplied by L_i). Row scaling changes neither the
 * row space nor the solutions of [A | b], and scales the determinant by ΠL_i.
 */
export function integerRows(pairs) {
  const scales = [];
  const rows = pairs.map((row) => {
    let L = 1n;
    for (const [, d] of row) if (d !== 1n) L = lcm(L, d);
    scales.push(L);
    return row.map(([n, d]) => (L === d ? n : n * (L / d)));
  });
  return { rows, scales };
}

/**
 * Fraction-free Gauss–Jordan (Bareiss): integer arithmetic only, every
 * division exact by Sylvester's identity, entries stay minors of the input so
 * they never blow up the way naive fraction elimination does. Pivots on the
 * first `columns` columns; the rest ride along (augmented parts).
 *
 * Returns { rows, pivots, last, sign }: after the run every pivot entry equals
 * `last`, so the reduced row echelon form is rows[k] / last.
 */
export function* bareiss(intRows, columns = intRows[0]?.length ?? 0) {
  const rows = intRows.map((r) => [...r]);
  const m = rows.length;
  const width = m ? rows[0].length : 0;
  const pivots = [];
  let prev = 1n;
  let sign = 1;
  let r = 0;
  for (let col = 0; col < columns && r < m; col += 1) {
    let found = -1;
    // smallest nonzero pivot keeps the intermediate products a little smaller
    for (let i = r; i < m; i += 1) {
      const v = rows[i][col];
      if (v !== 0n && (found < 0 || abs(v) < abs(rows[found][col]))) found = i;
    }
    if (found < 0) continue;
    if (found !== r) {
      [rows[r], rows[found]] = [rows[found], rows[r]];
      sign = -sign;
    }
    const pivotRow = rows[r];
    const p = pivotRow[col];
    for (let i = 0; i < m; i += 1) {
      if (i === r) continue;
      const row = rows[i];
      const f = row[col];
      for (let j = 0; j < width; j += 1) {
        if (j === col) continue;
        const a = row[j];
        const b = pivotRow[j];
        if (a === 0n && (f === 0n || b === 0n)) continue;
        row[j] = (p * a - f * b) / prev;
      }
      row[col] = 0n;
    }
    prev = p;
    pivots.push(col);
    r += 1;
    yield { step: "pivot", done: r, of: Math.min(m, columns) };
  }
  return { rows, pivots, last: prev, sign };
}

/** Reduced row echelon form as rows of [n, d] pairs plus pivot columns. */
export function* rrefPairs(pairs, columns) {
  const { rows: ints } = integerRows(pairs);
  const { rows, pivots, last } = yield* bareiss(ints, columns);
  // rows past the rank are zero in the pivot columns; in augmented columns
  // they keep a (scaled) remainder whose only meaning is zero or not zero
  const out = rows.map((row) => row.map((v) => reduce(v, last)));
  return { rows: out, pivots };
}

export function* determinantBareiss(pairs) {
  const n = pairs.length;
  if (n === 0) return [1n, 1n];
  const { rows: ints, scales } = integerRows(pairs);
  const { pivots, last, sign } = yield* bareiss(ints, n);
  if (pivots.length < n) return [0n, 1n];
  const denom = scales.reduce((a, b) => a * b, 1n);
  return reduce(sign < 0 ? -last : last, denom);
}

/**
 * Solve A x = B (B may have several columns). Returns
 *   { status: "unique" | "infinite" | "none", x, basis }
 * x: one solution as rows of pairs (free unknowns set to 0), or null;
 * basis: null-space vectors of A (each a column of pairs), empty when unique.
 */
export function* solvePairs(aPairs, bPairs) {
  const n = aPairs[0]?.length ?? 0;
  const k = bPairs[0]?.length ?? 0;
  const wide = aPairs.map((row, i) => [...row, ...bPairs[i]]);
  const { rows, pivots } = yield* rrefPairs(wide, n);
  for (let i = pivots.length; i < rows.length; i += 1) {
    if (rows[i].slice(n).some(([v]) => v !== 0n)) return { status: "none", x: null, basis: nullBasis(rows, pivots, n) };
  }
  const x = Array.from({ length: n }, () => Array.from({ length: k }, () => [0n, 1n]));
  pivots.forEach((col, r) => {
    for (let j = 0; j < k; j += 1) x[col][j] = rows[r][n + j];
  });
  const basis = nullBasis(rows, pivots, n);
  return { status: basis.length ? "infinite" : "unique", x, basis };
}

export function nullBasis(rows, pivots, n) {
  const isPivot = new Set(pivots);
  const basis = [];
  for (let f = 0; f < n; f += 1) {
    if (isPivot.has(f)) continue;
    const v = Array.from({ length: n }, () => [0n, 1n]);
    v[f] = [1n, 1n];
    pivots.forEach((col, r) => {
      const [a, b] = rows[r][f];
      v[col] = [-a, b];
    });
    basis.push(v);
  }
  return basis;
}

export function* inversePairs(pairs) {
  const n = pairs.length;
  const wide = pairs.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => [i === j ? 1n : 0n, 1n])]);
  const { rows, pivots } = yield* rrefPairs(wide, n);
  if (pivots.length < n) return null;
  return rows.map((row) => row.slice(n));
}

// ---------------------------------------------------------------------------
// Multi-modular determinant: the determinant of an integer matrix modulo many
// word-sized primes, glued together with the Chinese remainder theorem until
// the product of primes passes twice the Hadamard bound. All the O(n^3) work is
// on plain doubles (primes below 2^26, so a·b < 2^52 stays exact), which is
// far faster than BigInt elimination for anything past a few rows.

const PRIME_TOP = 2 ** 26;
const primes = [];

function isPrime(p) {
  if (p % 2 === 0) return p === 2;
  for (let q = 3; q * q <= p; q += 2) if (p % q === 0) return false;
  return true;
}

/** The i-th prime below 2^26, counting down. */
export function prime(i) {
  let p = primes.length ? primes[primes.length - 1] - 2 : PRIME_TOP - 1;
  while (primes.length <= i) {
    if (isPrime(p)) primes.push(p);
    p -= 2;
  }
  return primes[i];
}

export function invMod(a, p) {
  let [t, newT, r, newR] = [0, 1, p, a];
  while (newR !== 0) {
    const q = Math.floor(r / newR);
    [t, newT] = [newT, t - q * newT];
    [r, newR] = [newR, r - q * newR];
  }
  return t < 0 ? t + p : t;
}

/** Determinant mod p of an integer matrix (rows of BigInt). */
export function detMod(intRows, p) {
  const n = intRows.length;
  const P = BigInt(p);
  const a = new Float64Array(n * n);
  for (let i = 0; i < n; i += 1) {
    const row = intRows[i];
    for (let j = 0; j < n; j += 1) {
      const v = row[j];
      const r = v === 0n ? 0 : Number(v % P);
      a[i * n + j] = r < 0 ? r + p : r;
    }
  }
  const pinv = 1 / p;
  let det = 1;
  for (let col = 0; col < n; col += 1) {
    let found = -1;
    for (let i = col; i < n; i += 1) if (a[i * n + col] !== 0) { found = i; break; }
    if (found < 0) return 0;
    if (found !== col) {
      for (let j = col; j < n; j += 1) {
        const t = a[col * n + j];
        a[col * n + j] = a[found * n + j];
        a[found * n + j] = t;
      }
      det = p - det;
    }
    const base = col * n;
    const pivot = a[base + col];
    det = (det * pivot) % p;
    const inv = invMod(pivot, p);
    for (let i = col + 1; i < n; i += 1) {
      const ri = i * n;
      const f = (a[ri + col] * inv) % p;
      if (f === 0) continue;
      const g = p - f; // add g·pivotRow instead of subtracting f·pivotRow
      for (let j = col + 1; j < n; j += 1) {
        // v < 2^53, so v·(1/p) is off by at most one; fix up instead of fmod
        const v = a[ri + j] + g * a[base + j];
        let r = v - Math.floor(v * pinv) * p;
        if (r < 0) r += p;
        else if (r >= p) r -= p;
        a[ri + j] = r;
      }
    }
  }
  return det;
}

/** log2 of the Hadamard bound Π‖row‖₂, rounded up. */
export function hadamardBits(intRows) {
  let bits = 0;
  for (const row of intRows) {
    let s = 0n;
    for (const v of row) s += v * v;
    if (s === 0n) return 0; // a zero row: determinant is 0
    bits += Math.ceil(bitLength(s) / 2);
  }
  return bits;
}

export function* determinantModular(pairs) {
  const n = pairs.length;
  if (n === 0) return [1n, 1n];
  const { rows, scales } = integerRows(pairs);
  const need = hadamardBits(rows) + 2; // product of primes must exceed 2·bound
  let x = 0n;
  let M = 1n;
  let bits = 0;
  for (let i = 0; bits < need; i += 1) {
    const p = prime(i);
    const r = detMod(rows, p);
    // Garner step: x ≡ r (mod p), keeping x mod M
    const xm = Number(x % BigInt(p));
    const t = (((r - xm) % p + p) % p) * invMod(Number(M % BigInt(p)), p) % p;
    x += M * BigInt(t);
    M *= BigInt(p);
    bits += Math.log2(p);
    yield { step: "prime", done: Math.floor(bits), of: need };
  }
  if (x > M / 2n) x -= M;
  const denom = scales.reduce((a, b) => a * b, 1n);
  return reduce(x, denom);
}

/** True when A is certainly nonsingular (det ≠ 0 mod some prime). Never wrong when true. */
export function provablyNonsingular(pairs, tries = 2) {
  const { rows } = integerRows(pairs);
  for (let i = 0; i < tries; i += 1) if (detMod(rows, prime(i)) !== 0) return true;
  return false;
}

export { gcd };
