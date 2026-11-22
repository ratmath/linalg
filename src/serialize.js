// Turn results into plain data (for postMessage, JSON, logs) and back.
// Rationals become "n/d" strings, matrices { $: "matrix", rows }, algebraic
// numbers keep their polynomial and interval so they can still be refined.

import { Rational, RationalInterval } from "@ratmath/core";
import { Matrix } from "./matrix.js";
import { AlgebraicReal } from "./poly.js";
import { ComplexAlgebraic } from "./complex.js";
import { AlgebraicVector } from "./structure.js";
import { toPair } from "./q.js";

const ps = ([n, d]) => (d === 1n ? `${n}` : `${n}/${d}`);
const sp = (s) => toPair(s);

export function toPlain(v) {
  if (v === null || v === undefined || typeof v === "number" || typeof v === "string" || typeof v === "boolean") return v;
  if (typeof v === "bigint") return { $: "bigint", v: `${v}` };
  if (v instanceof Matrix) return { $: "matrix", rows: v.toStrings() };
  if (v instanceof Rational) return { $: "rational", v: v.toString() };
  if (v instanceof RationalInterval) return { $: "interval", lo: v.low.toString(), hi: v.high.toString() };
  if (v instanceof AlgebraicReal) {
    return v.isRational ? { $: "real", rational: ps(v.rational) } : { $: "real", poly: v.poly.map(String), lo: ps(v.lo), hi: ps(v.hi) };
  }
  if (v instanceof ComplexAlgebraic) return { $: "complex", poly: v.poly.map(String), re: v.re.map(ps), im: v.im.map(ps) };
  if (v instanceof AlgebraicVector) return { $: "avector", lambda: toPlain(v.lambda), entries: v.entries.map((p) => p.map(ps)) };
  if (typeof v === "function") return undefined;
  if (Array.isArray(v)) return v.map(toPlain);
  if (typeof v === "object") {
    const o = {};
    for (const [k, x] of Object.entries(v)) {
      const y = toPlain(x);
      if (y !== undefined) o[k] = y;
    }
    return o;
  }
  return v;
}

export function fromPlain(v) {
  if (v === null || typeof v !== "object") return v;
  if (Array.isArray(v)) return v.map(fromPlain);
  switch (v.$) {
    case "bigint": return BigInt(v.v);
    case "matrix": return Matrix.from(v.rows);
    case "rational": return new Rational(v.v);
    case "interval": return new RationalInterval(new Rational(v.lo), new Rational(v.hi));
    case "real":
      return v.rational ? AlgebraicReal.rational(sp(v.rational)) : new AlgebraicReal(v.poly.map(BigInt), sp(v.lo), sp(v.hi));
    case "complex": return new ComplexAlgebraic(v.poly.map(BigInt), v.re.map(sp), v.im.map(sp));
    case "avector": return new AlgebraicVector(fromPlain(v.lambda), v.entries.map((p) => p.map(sp)));
    default: {
      const o = {};
      for (const [k, x] of Object.entries(v)) o[k] = fromPlain(x);
      return o;
    }
  }
}
