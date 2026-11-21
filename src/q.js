// Exact rationals as [numerator, denominator] pairs of BigInt, denominator > 0
// and reduced. Kept as bare pairs (not objects) so the inner loops of the
// exact algorithms do not allocate classes or re-normalize on every step.

import { Rational } from "@ratmath/core";

export const abs = (a) => (a < 0n ? -a : a);

export function gcd(a, b) {
  a = abs(a);
  b = abs(b);
  while (b) [a, b] = [b, a % b];
  return a;
}

export const lcm = (a, b) => (a === 0n || b === 0n ? 0n : abs(a / gcd(a, b) * b));

/** Reduce n/d with d > 0. */
export function reduce(n, d) {
  if (d === 0n) throw new RangeError("Denominator cannot be zero");
  if (d < 0n) {
    n = -n;
    d = -d;
  }
  const g = gcd(n, d);
  return g > 1n ? [n / g, d / g] : [n, d];
}

/**
 * Read one entry. Accepts bigint, integer numbers, decimal numbers (read by
 * their shortest decimal spelling, so 0.1 is 1/10, which is what a person
 * typing 0.1 means), strings in Core's literal grammar ("3/4", "1.25",
 * "0.#3"), and anything with bigint numerator/denominator (Core Rational,
 * Integer via .value, Fraction).
 */
export function toPair(x) {
  if (typeof x === "bigint") return [x, 1n];
  if (typeof x === "number") {
    if (!Number.isFinite(x)) throw new RangeError(`Matrix entries must be finite, got ${x}`);
    if (Number.isInteger(x)) return [BigInt(x), 1n];
    return fromDecimalString(String(x));
  }
  if (typeof x === "string") {
    const r = new Rational(x.trim());
    return [r.numerator, r.denominator];
  }
  if (Array.isArray(x) && x.length === 2 && typeof x[0] === "bigint") return reduce(x[0], x[1]);
  if (x && typeof x === "object") {
    if (typeof x.numerator === "bigint" && typeof x.denominator === "bigint") return reduce(x.numerator, x.denominator);
    if (typeof x.value === "bigint") return [x.value, 1n];
  }
  throw new TypeError(`Cannot read ${String(x)} as an exact rational`);
}

function fromDecimalString(s) {
  // handles "1.25", "-3e-7", "1.5e+21"
  const m = /^(-?)(\d*)(?:\.(\d*))?(?:e([+-]?\d+))?$/i.exec(s);
  if (!m) throw new TypeError(`Cannot read ${s} as a decimal`);
  const [, sign, int, frac = "", exp = "0"] = m;
  let n = BigInt((int || "0") + frac);
  let e = Number(exp) - frac.length;
  let d = 1n;
  if (e >= 0) n *= 10n ** BigInt(e);
  else d = 10n ** BigInt(-e);
  return reduce(sign ? -n : n, d);
}

export const toRational = ([n, d]) => new Rational(n, d);

/** Nearest double to n/d, also for numerators and denominators beyond 2^53. */
export function pairToNumber(n, d) {
  if (n === 0n) return 0;
  const safe = 2n ** 53n;
  if (abs(n) < safe && d < safe) return Number(n) / Number(d);
  // keep about 64 significant bits of the quotient
  const shift = BigInt(bitLength(abs(n)) - bitLength(d) - 64);
  const q = shift >= 0n ? n / (d << shift) : (n << -shift) / d;
  return Number(q) * 2 ** Number(shift);
}

export function bitLength(a) {
  return a === 0n ? 0 : a.toString(2).length;
}

export const addQ = ([a, b], [c, d]) => reduce(a * d + c * b, b * d);
export const subQ = ([a, b], [c, d]) => reduce(a * d - c * b, b * d);
export const mulQ = ([a, b], [c, d]) => reduce(a * c, b * d);
export const divQ = ([a, b], [c, d]) => reduce(a * d, b * c);
export const pairToString = ([n, d]) => (d === 1n ? `${n}` : `${n}/${d}`);
