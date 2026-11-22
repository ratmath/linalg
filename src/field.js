// Field-generic matrix algorithms, for entries that are not plain rationals:
// ranges (rational intervals), numbers known only through enclosures, or
// anything else that can add, multiply and say whether it is zero. Lifted
// from the calculator app's engine so both use one copy.
//
//   field = { zero, one, add, sub, mul, div, neg,
//             sign(v) -> 0 (certainly zero), 1 (certainly not zero), or null (cannot tell) }
//
// Elimination has to know whether a pivot is zero. When that cannot be told
// (a range that contains 0, or a real too close to 0 to separate), the
// `undecided(row, col)` option is called; it should throw (for example asking
// for more digits) or the algorithm treats the column as having no pivot.

import { Rational, RationalInterval } from "@ratmath/core";

/** Reduced row echelon form. Returns { rows, pivots } with pivots[i] the pivot column of row i. */
export function rref(input, field, { columns = input[0]?.length ?? 0, undecided } = {}) {
  const rows = input.map((row) => [...row]);
  const pivots = [];
  let r = 0;
  for (let col = 0; col < columns && r < rows.length; col += 1) {
    let found = -1;
    let unsure = false;
    for (let i = r; i < rows.length; i += 1) {
      const sign = field.sign(rows[i][col]);
      if (sign === 1) {
        found = i;
        break;
      }
      if (sign === null) unsure = true;
    }
    if (found < 0) {
      if (unsure) undecided?.(r, col);
      continue;
    }
    [rows[r], rows[found]] = [rows[found], rows[r]];
    const pivot = rows[r][col];
    rows[r] = rows[r].map((x, j) => (j === col ? field.one : j < col ? field.zero : field.div(x, pivot)));
    for (let i = 0; i < rows.length; i += 1) {
      if (i === r) continue;
      const f = rows[i][col];
      if (field.sign(f) === 0) continue;
      rows[i] = rows[i].map((x, j) => (j === col ? field.zero : field.sub(x, field.mul(f, rows[r][j]))));
    }
    pivots.push(col);
    r += 1;
  }
  return { rows, pivots };
}

/** Determinant. Small matrices are expanded along rows, which needs no zero tests. */
export function det(m, field, options) {
  const n = m.length;
  if (n === 0) return field.one;
  if (n <= 5) return expand(m, field);
  // elimination, tracking row swaps and pivots
  const rows = m.map((row) => [...row]);
  let result = field.one;
  for (let col = 0; col < n; col += 1) {
    let found = -1;
    let unsure = false;
    for (let i = col; i < n; i += 1) {
      const sign = field.sign(rows[i][col]);
      if (sign === 1) {
        found = i;
        break;
      }
      if (sign === null) unsure = true;
    }
    if (found < 0) {
      if (unsure) options?.undecided?.(col, col);
      return field.zero;
    }
    if (found !== col) {
      [rows[col], rows[found]] = [rows[found], rows[col]];
      result = field.neg(result);
    }
    result = field.mul(result, rows[col][col]);
    for (let i = col + 1; i < n; i += 1) {
      const f = field.div(rows[i][col], rows[col][col]);
      for (let j = col; j < n; j += 1) rows[i][j] = field.sub(rows[i][j], field.mul(f, rows[col][j]));
    }
  }
  return result;
}

function expand(m, field) {
  const n = m.length;
  if (n === 1) return m[0][0];
  if (n === 2) return field.sub(field.mul(m[0][0], m[1][1]), field.mul(m[0][1], m[1][0]));
  let total = field.zero;
  for (let j = 0; j < n; j += 1) {
    if (field.sign(m[0][j]) === 0) continue;
    const minor = m.slice(1).map((row) => row.filter((_, k) => k !== j));
    const term = field.mul(m[0][j], expand(minor, field));
    total = j % 2 === 0 ? field.add(total, term) : field.sub(total, term);
  }
  return total;
}

export const transpose = (m) => (m.length ? m[0].map((_, j) => m.map((row) => row[j])) : []);

export function multiply(a, b, field) {
  return a.map((row) => b[0].map((_, j) => row.reduce((total, x, k) => field.add(total, field.mul(x, b[k][j])), field.zero)));
}

export const identity = (n, field) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? field.one : field.zero)));

/** The inverse, or null when the matrix is singular. */
export function inverse(m, field, options) {
  const n = m.length;
  const wide = m.map((row, i) => [...row, ...identity(n, field)[i]]);
  const { rows, pivots } = rref(wide, field, { ...options, columns: n });
  if (pivots.length < n) return null;
  return rows.map((row) => row.slice(n));
}

/** Core Rationals as a field. */
export const rationalField = {
  zero: Rational.zero,
  one: Rational.one,
  add: (a, b) => a.add(b),
  sub: (a, b) => a.subtract(b),
  mul: (a, b) => a.multiply(b),
  div: (a, b) => a.divide(b),
  neg: (a) => a.negate(),
  sign: (a) => (a.numerator === 0n ? 0 : 1),
};

/**
 * Core RationalIntervals as a (non-exact) field: results are enclosures.
 * A pivot interval that contains 0 but is not [0, 0] is undecided.
 */
export const intervalField = {
  zero: new RationalInterval(0, 0),
  one: new RationalInterval(1, 1),
  add: (a, b) => a.add(b),
  sub: (a, b) => a.subtract(b),
  mul: (a, b) => a.multiply(b),
  div: (a, b) => a.divide(b),
  neg: (a) => a.negate(),
  sign: (a) => (!a.containsZero() ? 1 : a.low.numerator === 0n && a.high.numerator === 0n ? 0 : null),
};

/**
 * Integers modulo a prime p as a field (elements are BigInts in [0, p)).
 * `from(x)` maps a rational (anything with numerator/denominator, a bigint or
 * a number) into the field; a denominator divisible by p has no image.
 */
export function primeField(p) {
  p = BigInt(p);
  if (p < 2n || !isProbablePrime(p)) throw new RangeError(`${p} is not prime`);
  const mod = (x) => ((x % p) + p) % p;
  const inv = (a) => {
    if (a === 0n) throw new RangeError(`division by zero mod ${p}`);
    let [r0, r1, s0, s1] = [p, a, 0n, 1n];
    while (r1) {
      const q = r0 / r1;
      [r0, r1, s0, s1] = [r1, r0 - q * r1, s1, s0 - q * s1];
    }
    return mod(s0);
  };
  return {
    p,
    zero: 0n,
    one: 1n,
    add: (a, b) => mod(a + b),
    sub: (a, b) => mod(a - b),
    mul: (a, b) => mod(a * b),
    div: (a, b) => mod(a * inv(b)),
    neg: (a) => mod(-a),
    inv,
    sign: (a) => (a === 0n ? 0 : 1),
    from(x) {
      if (typeof x === "bigint") return mod(x);
      if (typeof x === "number") return mod(BigInt(x));
      const d = mod(x.denominator);
      if (d === 0n) throw new RangeError(`${x} has a denominator divisible by ${p}, so it has no value mod ${p}`);
      return mod(mod(x.numerator) * inv(d));
    },
  };
}

function isProbablePrime(n) {
  if (n < 4n) return n >= 2n;
  if (n % 2n === 0n) return false;
  const powmod = (b, e, m) => {
    let r = 1n;
    b %= m;
    while (e) {
      if (e & 1n) r = (r * b) % m;
      b = (b * b) % m;
      e >>= 1n;
    }
    return r;
  };
  let d = n - 1n;
  let s = 0;
  while (d % 2n === 0n) {
    d /= 2n;
    s += 1;
  }
  for (const a of [2n, 3n, 5n, 7n, 11n, 13n, 17n, 19n, 23n, 29n, 31n, 37n]) {
    if (a % n === 0n) continue;
    let x = powmod(a, d, n);
    if (x === 1n || x === n - 1n) continue;
    let composite = true;
    for (let i = 1; i < s; i += 1) {
      x = (x * x) % n;
      if (x === n - 1n) {
        composite = false;
        break;
      }
    }
    if (composite) return false;
  }
  return true;
}

/** Real algebraic numbers (AlgebraicReal) as an exact field; zero tests are exact. */
export const algebraicField = {
  get zero() { return ALG.fromRational([0n, 1n]); },
  get one() { return ALG.fromRational([1n, 1n]); },
  add: (a, b) => ALG.add(a, b),
  sub: (a, b) => ALG.sub(a, b),
  mul: (a, b) => ALG.mul(a, b),
  div: (a, b) => ALG.div(a, b),
  neg: (a) => ALG.neg(a),
  sign: (a) => (ALG.sign(a) === 0 ? 0 : 1),
};
import * as ALG from "./algebraic.js";
