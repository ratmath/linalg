// The workbench's little language: exact rational scalars and matrices,
// arithmetic, and the library's exact operations as functions.
//
//   A*B - 2A'      [1, 1/2; 0.#3, 2..1/4]      A^-1      C = rref(A)
//
// Values are Core Rationals, Matrices, Poly (rational coefficients, constant
// first), Rec (labelled fields) and Note (a message). Everything is exact; the
// operations that return irrational numbers are not offered here.
import * as L from "./linalg.js?v=5591b17725";

const { Matrix, Rational, RationalInterval, AlgebraicReal, ComplexAlgebraic, parseNumber } = L;
const ALG = L.algebraic;

/**
 * Which numbers are allowed. domain: "rational" (optionally with real
 * algebraic numbers, and complex ones on top), "ranges" (rational intervals)
 * or "modp" (integers modulo the prime p).
 */
export let CTX = { domain: "rational", algebraic: false, complex: false, p: 7n };
export function setContext(ctx) {
  CTX = { ...CTX, ...ctx, p: BigInt(ctx.p ?? CTX.p) };
  FP_CACHE = null;
}

/** An integer modulo p. */
export class ModP {
  constructor(v, p) {
    this.v = v;
    this.p = p;
  }
  toString() {
    return String(this.v);
  }
}

/** A matrix whose entries live in a field other than Q (algebraic, ranges, mod p). */
export class FMatrix {
  constructor(F, a) {
    this.F = F;
    this.a = a;
    this.m = a.length;
    this.n = a[0]?.length ?? 0;
  }
}

const isA = (v) => v instanceof AlgebraicReal;
const isI = (v) => v instanceof RationalInterval;
const isP = (v) => v instanceof ModP;
const isF = (v) => v instanceof FMatrix;
const isC = (v) => v instanceof ComplexAlgebraic;

const FA = {
  ...L.field.algebraicField,
  tag: "real algebraic",
  of: (x) => (isA(x) ? x : isQ(x) ? ALG.fromRational([x.numerator, x.denominator]) : fail(`a ${kind(x)} is not a real algebraic number`)),
  wrap: (v) => v,
};
const FI = {
  ...L.field.intervalField,
  tag: "ranges",
  of: (x) => (isI(x) ? x : isQ(x) ? new RationalInterval(x, x) : fail(`a ${kind(x)} cannot be used as a range`)),
  wrap: (v) => v,
};
let FP_CACHE = null;
function FP() {
  if (FP_CACHE) return FP_CACHE;
  let F;
  try {
    F = L.field.primeField(CTX.p);
  } catch (e) {
    fail(e.message);
  }
  FP_CACHE = {
    ...F,
    tag: `integers mod ${F.p}`,
    of: (x) => {
      if (isP(x)) return x.v;
      if (isQ(x)) {
        try {
          return F.from(x);
        } catch (e) {
          return fail(e.message);
        }
      }
      return fail(`a ${kind(x)} has no value mod ${F.p}`);
    },
    wrap: (v) => new ModP(v, F.p),
  };
  return FP_CACHE;
}
/** The field for the current setting, used to lift rationals. */
const ctxField = () => (CTX.domain === "modp" ? FP() : CTX.domain === "ranges" ? FI : CTX.algebraic ? FA : null);

export class Poly {
  constructor(coeffs, variable = "x") {
    this.c = coeffs.map(toRational);
    this.variable = variable;
  }
}
export class Rec {
  constructor(fields, note = null) {
    this.fields = fields.filter(([, v]) => v !== undefined);
    this.note = note;
  }
}
export class Note {
  constructor(text, link = null) {
    this.text = text;
    this.link = link;
  }
}
export class EvalError extends Error {}

const fail = (msg) => {
  throw new EvalError(msg);
};

export function toRational(x) {
  if (x instanceof Rational) return x;
  if (typeof x === "bigint") return new Rational(x, 1n);
  if (typeof x === "number") return new Rational(String(x));
  if (x && typeof x.numerator === "bigint") return new Rational(x.numerator, x.denominator);
  if (x && typeof x.value === "bigint") return new Rational(x.value, 1n);
  return fail(`not a rational number: ${String(x)}`);
}

// ---------------------------------------------------------------------------
// Tokens

const NUM = /^(?:\d[\d_]*\.\.\d+\/\d+|\d[\d_]*\.~[\d~]*|(?:\d[\d_]*)?\.?[\d_]*(?:#\d+)?)/;

function tokenize(src) {
  const out = [];
  let i = 0;
  while (i < src.length) {
    const ch = src[i];
    if (ch === " " || ch === "\t" || ch === "\r") { i += 1; continue; }
    if (ch === "\n") { out.push({ t: ";", at: i }); i += 1; continue; }
    const prev = out[out.length - 1];
    const valueBefore = prev && (prev.t === "num" || prev.t === "id" || prev.t === ")" || prev.t === "]" || prev.t === "'");
    if (CTX.domain === "ranges" && (/[0-9.]/.test(ch) || (ch === "-" && !valueBefore && /^-[0-9.][^:]*:/.test(src.slice(i, i + 40))))) {
      const m = src.slice(i).match(/^-?(?:\d[\d_]*\.\.\d+\/\d+|(?:\d[\d_]*)?\.?[\d_]*(?:#\d+)?(?:\/\d+)?):-?(?:\d[\d_]*\.\.\d+\/\d+|(?:\d[\d_]*)?\.?[\d_]*(?:#\d+)?(?:\/\d+)?)/);
      if (m && /\d/.test(m[0].split(":")[0]) && /\d/.test(m[0].split(":")[1])) {
        let v;
        try {
          v = parseNumber(m[0].replace(/_/g, ""));
        } catch (e) {
          fail(`"${m[0]}" is not a range (${e.message})`);
        }
        out.push({ t: "num", v: isI(v) ? v : toRational(v), at: i });
        i += m[0].length;
        continue;
      }
    }
    if (/[0-9.]/.test(ch)) {
      const m = src.slice(i).match(NUM);
      const text = m && m[0];
      if (!text || text === ".") fail(`unexpected "${ch}" at ${i + 1}`);
      let v;
      try {
        v = parseNumber(text.replace(/_/g, ""));
      } catch (e) {
        fail(`"${text}" is not a number (${e.message})`);
      }
      if (v && v.low !== undefined) fail(`"${text}" is a range; turn on Ranges to use them`);
      out.push({ t: "num", v: toRational(v), at: i });
      i += text.length;
      continue;
    }
    if (/[A-Za-z_]/.test(ch)) {
      const m = src.slice(i).match(/^[A-Za-z_][A-Za-z0-9_]*/)[0];
      out.push({ t: "id", v: m, at: i });
      i += m.length;
      continue;
    }
    if ("+-*/^()[],;'=·×−".includes(ch)) {
      const t = ch === "·" || ch === "×" ? "*" : ch === "−" ? "-" : ch;
      out.push({ t, at: i });
      i += 1;
      continue;
    }
    fail(`unexpected "${ch}" at ${i + 1}`);
  }
  out.push({ t: "end", at: src.length });
  return out;
}

// ---------------------------------------------------------------------------
// Parser to a small tree

export function parse(src) {
  const toks = tokenize(src);
  let k = 0;
  const peek = () => toks[k];
  const take = (t) => {
    if (toks[k].t !== t) fail(`expected "${t}" at ${toks[k].at + 1}`);
    return toks[k++];
  };
  const startsPrimary = (tok) => tok.t === "num" || tok.t === "id" || tok.t === "(" || tok.t === "[";

  function statement() {
    if (peek().t === "id" && toks[k + 1].t === "=") {
      const name = take("id").v;
      take("=");
      return { k: "assign", name, e: expr() };
    }
    return expr();
  }
  function expr() {
    let a = term();
    while (peek().t === "+" || peek().t === "-") {
      const op = toks[k++].t;
      a = { k: op, a, b: term() };
    }
    return a;
  }
  function term() {
    let a = unary();
    for (;;) {
      if (peek().t === "*" || peek().t === "/") {
        const op = toks[k++].t;
        a = { k: op, a, b: unary() };
      } else if (startsPrimary(peek())) {
        a = { k: "*", a, b: power() }; // 2A, (A+B)C
      } else return a;
    }
  }
  function unary() {
    if (peek().t === "-") { k += 1; return { k: "neg", a: unary() }; }
    if (peek().t === "+") { k += 1; return unary(); }
    return power();
  }
  function power() {
    const base = postfix();
    if (peek().t === "^") { k += 1; return { k: "^", a: base, b: unary() }; }
    return base;
  }
  function postfix() {
    let a = primary();
    while (peek().t === "'") { k += 1; a = { k: "T", a }; }
    return a;
  }
  function primary() {
    const tok = peek();
    if (tok.t === "num") { k += 1; return { k: "num", v: tok.v }; }
    if (tok.t === "id") {
      k += 1;
      if (peek().t === "(") {
        k += 1;
        const args = [];
        if (peek().t !== ")") {
          args.push(expr());
          while (peek().t === ",") { k += 1; args.push(expr()); }
        }
        take(")");
        return { k: "call", name: tok.v, args };
      }
      return { k: "var", name: tok.v };
    }
    if (tok.t === "(") {
      k += 1;
      const e = expr();
      take(")");
      return e;
    }
    if (tok.t === "[") {
      k += 1;
      const rows = [[]];
      while (peek().t === ";") k += 1;
      if (peek().t !== "]") rows[0].push(expr());
      for (;;) {
        const t = peek().t;
        if (t === ",") { k += 1; rows.at(-1).push(expr()); } else if (t === ";") {
          k += 1;
          if (peek().t === ";" || peek().t === "]") continue;
          rows.push([expr()]);
        } else break;
      }
      take("]");
      return { k: "mat", rows: rows.filter((r) => r.length) };
    }
    if (tok.t === "end") fail("the expression ends too soon");
    return fail(`unexpected "${tok.t}" at ${tok.at + 1}`);
  }

  const tree = statement();
  if (peek().t !== "end") fail(`unexpected "${peek().t === "num" || peek().t === "id" ? String(peek().v) : peek().t}" at ${peek().at + 1}`);
  return tree;
}

// ---------------------------------------------------------------------------
// Values and arithmetic

const isM = (v) => v instanceof Matrix;
const isQ = (v) => v instanceof Rational;
const kind = (v) => (isM(v) ? `${v.rows}×${v.cols} matrix` : isF(v) ? `${v.m}×${v.n} matrix over the ${v.F.tag}` : isQ(v) ? "number" : isA(v) ? "real algebraic number" : isI(v) ? "range" : isP(v) ? `number mod ${v.p}` : isC(v) ? "complex number" : v instanceof Poly ? "polynomial" : "result");
const needM = (v, what) => (isM(v) ? v : fail(`${what} needs a matrix, got a ${kind(v)}`));
const needQ = (v, what) => (isQ(v) ? v : isM(v) && v.rows === 1 && v.cols === 1 ? v.get(0, 0) : fail(`${what} needs a number, got a ${kind(v)}`));
const needInt = (v, what) => {
  const q = needQ(v, what);
  if (q.denominator !== 1n) fail(`${what} needs a whole number, got ${q}`);
  return Number(q.numerator);
};
const needSquare = (A, what) => (A.isSquare ? A : fail(`${what} needs a square matrix, got ${A.rows}×${A.cols}`));
const needIntegerEntries = (A, what) => {
  for (const row of A.toArray()) for (const x of row) if (x.denominator !== 1n) fail(`${what} needs integer entries; ${x} is not one`);
  return A;
};

function add(a, b, sign) {
  if (!pureQ(a) || !pureQ(b)) return gAdd(a, b, sign);
  if (isQ(a) && isQ(b)) return sign > 0 ? a.add(b) : a.subtract(b);
  if (isM(a) && isM(b)) {
    if (a.rows !== b.rows || a.cols !== b.cols) fail(`can't ${sign > 0 ? "add" : "subtract"} a ${a.rows}×${a.cols} and a ${b.rows}×${b.cols} matrix`);
    return sign > 0 ? a.add(b) : a.subtract(b);
  }
  return fail(`can't ${sign > 0 ? "add" : "subtract"} a ${kind(a)} and a ${kind(b)} (for A + k·I write A + k*I(n))`);
}
function mul(a, b) {
  if (!pureQ(a) || !pureQ(b)) return gMul(a, b);
  if (isQ(a) && isQ(b)) return a.multiply(b);
  if (isQ(a) && isM(b)) return b.scale(a);
  if (isM(a) && isQ(b)) return a.scale(b);
  if (isM(a) && isM(b)) {
    if (a.cols !== b.rows) fail(`can't multiply a ${a.rows}×${a.cols} by a ${b.rows}×${b.cols} matrix: the inner sizes differ`);
    return a.multiply(b);
  }
  return fail(`can't multiply a ${kind(a)} and a ${kind(b)}`);
}
function div(a, b) {
  if (!pureQ(a) || !pureQ(b)) return gDiv(a, b);
  if (isQ(b) && b.numerator === 0n) fail("division by zero");
  if (isQ(a) && isQ(b)) return a.divide(b);
  if (isM(a) && isQ(b)) return a.scale(new Rational(1n, 1n).divide(b));
  return fail(`can't divide a ${kind(a)} by a ${kind(b)}; for matrices use inv(B) or solve(A, b)`);
}


// ---------------------------------------------------------------------------
// Generic arithmetic for the other number settings

const pureQ = (v) => isQ(v) || isM(v);
const isScalar = (v) => isQ(v) || isA(v) || isI(v) || isP(v);

/** Field for combining a and b (either may be rational). */
function fieldFor(a, b) {
  for (const v of [a, b]) if (isF(v)) return v.F;
  if (isP(a) || isP(b) || CTX.domain === "modp") return FP();
  if (isI(a) || isI(b)) return FI;
  if (isA(a) || isA(b)) return FA;
  return ctxField() ?? fail("no field");
}
const toF = (F, v) => (isF(v) ? (v.F === F ? v : fail(`can't mix a matrix over the ${v.F.tag} with ${F.tag}`)) : isM(v) ? new FMatrix(F, v.toArray().map((r) => r.map((x) => F.of(x)))) : fail(`a ${kind(v)} is not a matrix`));
const norm = (v) => {
  if (isA(v) && v.isRational) return toRational(v.toRational());
  if (isI(v) && v.low.equals(v.high)) return v.low;
  return v;
};
const wrapS = (F, x) => norm(F.wrap(x));
const noComplex = (a, b) => {
  if (isC(a) || isC(b)) fail("arithmetic with complex numbers isn't built yet; they can be shown and refined");
};

function gAdd(a, b, sign) {
  noComplex(a, b);
  const F = fieldFor(a, b);
  const op = sign > 0 ? F.add : F.sub;
  if (isScalar(a) && isScalar(b)) return wrapS(F, op(F.of(a), F.of(b)));
  if (isScalar(a) || isScalar(b)) fail(`can't ${sign > 0 ? "add" : "subtract"} a ${kind(a)} and a ${kind(b)}`);
  const A = toF(F, a);
  const B = toF(F, b);
  if (A.m !== B.m || A.n !== B.n) fail(`can't ${sign > 0 ? "add" : "subtract"} a ${A.m}×${A.n} and a ${B.m}×${B.n} matrix`);
  return new FMatrix(F, A.a.map((r, i) => r.map((x, j) => op(x, B.a[i][j]))));
}
function gMul(a, b) {
  noComplex(a, b);
  const F = fieldFor(a, b);
  if (isScalar(a) && isScalar(b)) return wrapS(F, F.mul(F.of(a), F.of(b)));
  if (isScalar(a)) {
    const k = F.of(a);
    const B = toF(F, b);
    return new FMatrix(F, B.a.map((r) => r.map((x) => F.mul(k, x))));
  }
  if (isScalar(b)) return gMul(b, a);
  const A = toF(F, a);
  const B = toF(F, b);
  if (A.n !== B.m) fail(`can't multiply a ${A.m}×${A.n} by a ${B.m}×${B.n} matrix: the inner sizes differ`);
  return new FMatrix(F, L.field.multiply(A.a, B.a, F));
}
function gInv(F, x) {
  if (F.sign(x) === 0) fail("division by zero");
  if (F.sign(x) === null) fail("division by a range that contains 0");
  return F.div(F.one, x);
}
function gDiv(a, b) {
  noComplex(a, b);
  const F = fieldFor(a, b);
  if (!isScalar(b)) fail(`can't divide by a ${kind(b)}; for matrices use inv(B) or solve(A, b)`);
  const ib = gInv(F, F.of(b));
  if (isScalar(a)) return wrapS(F, F.mul(F.of(a), ib));
  const A = toF(F, a);
  return new FMatrix(F, A.a.map((r) => r.map((x) => F.mul(x, ib))));
}
function gNeg(a) {
  noComplex(a);
  const F = fieldFor(a, a);
  if (isScalar(a)) return wrapS(F, F.neg(F.of(a)));
  const A = toF(F, a);
  return new FMatrix(F, A.a.map((r) => r.map((x) => F.neg(x))));
}
const gT = (A) => new FMatrix(A.F, L.field.transpose(A.a));
function gPow(a, e) {
  noComplex(a);
  const F = fieldFor(a, a);
  if (isScalar(a)) {
    let x = F.of(a);
    if (e < 0) {
      x = gInv(F, x);
      e = -e;
    }
    let out = F.one;
    while (e > 0) {
      if (e & 1) out = F.mul(out, x);
      e >>= 1;
      if (e) x = F.mul(x, x);
    }
    return wrapS(F, out);
  }
  let A = toF(F, a);
  if (A.m !== A.n) fail(`a matrix power needs a square matrix, got ${A.m}×${A.n}`);
  if (e < 0) {
    const inv = L.field.inverse(A.a, F, { undecided: undecidedPivot });
    if (!inv) return new Note("The matrix is singular, so negative powers do not exist.");
    A = new FMatrix(F, inv);
    e = -e;
  }
  let out = L.field.identity(A.m, F);
  let base = A.a;
  while (e > 0) {
    if (e & 1) out = L.field.multiply(out, base, F);
    e >>= 1;
    if (e) base = L.field.multiply(base, base, F);
  }
  return new FMatrix(F, out);
}
function undecidedPivot() {
  fail("a pivot range contains 0, so it can't be told whether it is zero; narrow the ranges");
}
const fSquare = (A, what) => (A.m === A.n ? A : fail(`${what} needs a square matrix, got ${A.m}×${A.n}`));
const fCols = (F, vs, m) => new FMatrix(F, Array.from({ length: m }, (_, i) => vs.map((v) => v[i])));

/** Implementations of the elimination functions for matrices over any field. */
const GENERIC = {
  det: ([A]) => wrapS(A.F, L.field.det(fSquare(A, "det").a, A.F, { undecided: undecidedPivot })),
  rank: ([A]) => String(L.field.rref(A.a, A.F, { undecided: undecidedPivot }).pivots.length),
  rref: ([A]) => {
    const r = L.field.rref(A.a, A.F, { undecided: undecidedPivot });
    return new Rec([["R", new FMatrix(A.F, r.rows)], ["pivot columns", r.pivots.length ? r.pivots.map((j) => j + 1).join(", ") : "none"]]);
  },
  inv: ([A]) => {
    const inv = L.field.inverse(fSquare(A, "inv").a, A.F, { undecided: undecidedPivot });
    return inv ? new FMatrix(A.F, inv) : new Note("The matrix is singular, so it has no inverse.");
  },
  solve: async ([A, b], ex) => {
    if (A.F === FI) {
      if (b.n !== 1) fail("with ranges, solve takes one right-hand column");
      const r = await ex(L.solveIntervals(A.a, b.a.map((row) => row[0])));
      if (r.status !== "unique") return new Note("Could not enclose the solutions: some matrix in the ranges may be singular. Narrow the ranges.");
      return new Rec([["x", new FMatrix(FI, r.box.map((x) => [x]))]], "Every solution, for every choice of entries inside the ranges, lies in these ranges (and each such matrix is nonsingular).");
    }
    if (b.m !== A.m) fail(`solve: A has ${A.m} rows, b has ${b.m}`);
    const n = A.n;
    const { rows, pivots } = L.field.rref(A.a.map((r, i) => [...r, ...b.a[i]]), A.F, { columns: n, undecided: undecidedPivot });
    for (let i = pivots.length; i < rows.length; i += 1) if (rows[i].slice(n).some((x) => A.F.sign(x) !== 0)) return new Note("No solution: the system is inconsistent.");
    const x = Array.from({ length: n }, () => b.a[0].map(() => A.F.zero));
    pivots.forEach((c, i) => { x[c] = rows[i].slice(n); });
    const free = [...Array(n).keys()].filter((j) => !pivots.includes(j));
    const basis = free.map((f) => {
      const v = Array(n).fill(A.F.zero);
      v[f] = A.F.one;
      pivots.forEach((c, i) => { v[c] = A.F.neg(rows[i][f]); });
      return v;
    });
    if (!free.length) return new Rec([["x", new FMatrix(A.F, x)]], "The solution is unique.");
    return new Rec([["particular x", new FMatrix(A.F, x)], ["plus any combination of", fCols(A.F, basis, n)]], "Infinitely many solutions: x plus any combination of the columns shown.");
  },
  null: ([A]) => {
    const { rows, pivots } = L.field.rref(A.a, A.F, { undecided: undecidedPivot });
    const free = [...Array(A.n).keys()].filter((j) => !pivots.includes(j));
    if (!free.length) return new Note("The null space is {0}: the columns are independent.");
    return fCols(A.F, free.map((f) => {
      const v = Array(A.n).fill(A.F.zero);
      v[f] = A.F.one;
      pivots.forEach((c, i) => { v[c] = A.F.neg(rows[i][f]); });
      return v;
    }), A.n);
  },
  trace: ([A]) => wrapS(A.F, fSquare(A, "trace").a.reduce((t, r, i) => A.F.add(t, r[i]), A.F.zero)),
  T: ([A]) => gT(A),
  entry: ([A, i, j]) => {
    const r = needInt(i, "entry");
    const c = needInt(j, "entry");
    if (r < 1 || r > A.m || c < 1 || c > A.n) fail(`entry(${r}, ${c}) is outside a ${A.m}×${A.n} matrix`);
    return wrapS(A.F, A.a[r - 1][c - 1]);
  },
};

function hcat(parts) {
  if (!parts.every(pureQ)) {
    const F = parts.reduce((f, p) => (isF(p) ? p.F : isP(p) ? FP() : isI(p) ? FI : isA(p) ? FA : f), null) ?? ctxField();
    const ms = parts.map((p) => (isScalar(p) ? [[F.of(p)]] : toF(F, p).a));
    if (ms.some((m) => m.length !== ms[0].length)) fail("entries side by side in [ … ] must have the same number of rows");
    return new FMatrix(F, ms[0].map((_, i) => ms.flatMap((m) => m[i])));
  }
  const ms = parts.map((p) => (isQ(p) ? Matrix.from([[p]]) : needM(p, "a matrix literal")));
  const rows = ms[0].rows;
  if (ms.some((m) => m.rows !== rows)) fail("entries side by side in [ … ] must have the same number of rows");
  return ms.reduce((x, y) => x.augment(y));
}
function vcat(ms) {
  if (!ms.every(pureQ)) {
    const F = ms.find(isF).F;
    const fs = ms.map((m) => toF(F, m));
    if (fs.some((m) => m.n !== fs[0].n)) fail("rows in [ … ] must all have the same length");
    return new FMatrix(F, fs.flatMap((m) => m.a));
  }
  const cols = ms[0].cols;
  if (ms.some((m) => m.cols !== cols)) fail("rows in [ … ] must all have the same length");
  return Matrix.from(ms.flatMap((m) => m.toArray()));
}

/** Columns as a single matrix (handy for bases); empty → a note. */
function columns(vs, rows, empty) {
  if (!vs.length) return new Note(empty);
  return vs.reduce((x, y) => x.augment(y));
}

const permutationMatrix = (perm) => Matrix.from(perm.map((p) => perm.map((_, j) => (j === p ? 1 : 0))));

// ---------------------------------------------------------------------------
// Functions

const EIGEN_PAGE = "#eigen";

/** name → { args, help, run(args, exact) } ; exact(r) settles a Refinement. */
export const FUNCTIONS = {
  det: { args: "A", help: "determinant", run: async ([A], ex) => ex(L.det(needSquare(needM(A, "det"), "det"))) },
  rank: { args: "A", help: "rank", run: async ([A], ex) => toRational(BigInt(await ex(L.rank(needM(A, "rank"))))) },
  trace: { args: "A", help: "sum of the diagonal", run: async ([A], ex) => toRational(await ex(L.trace(needSquare(needM(A, "trace"), "trace")))) },
  rref: {
    args: "A",
    help: "reduced row echelon form",
    run: async ([A], ex) => {
      const r = await ex(L.rref(needM(A, "rref")));
      return new Rec([["R", r.matrix], ["pivot columns", r.pivots.length ? r.pivots.map((j) => j + 1).join(", ") : "none"]]);
    },
  },
  inv: {
    args: "A",
    help: "inverse",
    run: async ([A], ex) => (await ex(L.inverse(needSquare(needM(A, "inv"), "inv")))) ?? new Note("The matrix is singular (determinant 0), so it has no inverse. Try pinv(A)."),
  },
  solve: {
    args: "A, b",
    help: "every solution of A x = b",
    run: async ([A, b], ex) => {
      const s = await ex(L.solve(needM(A, "solve"), needM(b, "solve's right side")));
      if (s.status === "none") return new Note("No solution: the system is inconsistent. lsq(A, b) gives the least-squares answer.");
      if (s.status === "unique") return new Rec([["x", s.x]], "The solution is unique.");
      return new Rec([["particular x", s.x], ["plus any combination of", columns(s.basis, A.cols, "")]], "Infinitely many solutions: x plus any combination of the columns shown.");
    },
  },
  null: {
    args: "A",
    help: "null-space basis (as columns)",
    run: async ([A], ex) => columns(await ex(L.nullspace(needM(A, "null"))), A.cols, "The null space is {0}: the columns are independent."),
  },
  col: { args: "A", help: "column-space basis (pivot columns of A)", run: async ([A]) => columns(L.spaces(needM(A, "col")).columnSpace, A.rows, "The column space is {0}.") },
  row: {
    args: "A",
    help: "row-space basis (nonzero rows of rref)",
    run: async ([A]) => {
      const r = L.spaces(needM(A, "row")).rowSpace;
      return r.length ? vcat(r) : new Note("The row space is {0}.");
    },
  },
  spaces: {
    args: "A",
    help: "rank, nullity and the four subspaces",
    run: async ([A], ex) => {
      A = needM(A, "spaces");
      const s = L.spaces(A);
      const left = L.spaces(A.transpose());
      const nul = await ex(L.nullspace(A));
      const leftNul = await ex(L.nullspace(A.transpose()));
      return new Rec([
        ["rank", String(s.rank)],
        ["nullity", String(s.nullity)],
        ["column space", columns(s.columnSpace, A.rows, "{0}")],
        ["row space", s.rowSpace.length ? vcat(s.rowSpace) : new Note("{0}")],
        ["null space", columns(nul, A.cols, "{0}")],
        ["left null space", columns(leftNul, A.rows, "{0}")],
      ], `rank ${s.rank} + nullity ${s.nullity} = ${A.cols} columns; left nullity ${A.rows - left.rank}.`);
    },
  },
  T: { args: "A", help: "transpose (also A')", run: async ([A]) => (isQ(A) ? A : needM(A, "T").transpose()) },
  I: { args: "n", help: "n×n identity", run: async ([n]) => Matrix.identity(checkSize(needInt(n, "I"))) },
  zeros: { args: "r, c", help: "zero matrix", run: async ([r, c]) => Matrix.zeros(checkSize(needInt(r, "zeros")), checkSize(c === undefined ? needInt(r, "zeros") : needInt(c, "zeros"))) },
  lu: {
    args: "A",
    help: "P A = L U",
    run: async ([A], ex) => {
      const r = await ex(L.lu(needM(A, "lu")));
      return new Rec([["P", permutationMatrix(r.perm)], ["L", r.L], ["U", r.U]], "P·A = L·U, L unit lower triangular.");
    },
  },
  ldl: {
    args: "A",
    help: "A = L D Lᵀ (symmetric A)",
    run: async ([A], ex) => {
      A = needSquare(needM(A, "ldl"), "ldl");
      if (!A.equals(A.transpose())) fail("ldl needs a symmetric matrix");
      const r = await ex(L.ldl(A));
      if (!r) return new Note("No LDLᵀ without pivoting: a zero pivot appears. lu(A) works for any matrix.");
      return new Rec([["L", r.L], ["D", r.D], ["positive definite", r.positiveDefinite ? "yes" : "no"]], "A = L·D·Lᵀ.");
    },
  },
  qr: {
    args: "A",
    help: "A = Q R, square-root-free",
    run: async ([A], ex) => {
      const r = await ex(L.qr(needM(A, "qr")));
      return new Rec([["Q", r.Q], ["R", r.R], ["‖qᵢ‖²", r.normsSquared.map(toRational)]], "A = Q·R with Q's columns orthogonal (not unit length, which would need square roots) and R unit upper triangular.");
    },
  },
  gs: { args: "A", help: "Gram–Schmidt: orthogonal basis of the column space", run: async ([A], ex) => columns(await ex(L.gramSchmidt(needM(A, "gs"))), A.rows, "The column space is {0}.") },
  pinv: { args: "A", help: "Moore–Penrose pseudoinverse", run: async ([A], ex) => ex(L.pseudoinverse(needM(A, "pinv"))) },
  lsq: {
    args: "A, b",
    help: "least-squares solution of A x ≈ b",
    run: async ([A, b], ex) => {
      const r = await ex(L.leastSquares(needM(A, "lsq"), needM(b, "lsq's right side")));
      const res = r.residual;
      const sq = res.toArray().flat().reduce((acc, q) => acc.add(q.multiply(q)), new Rational(0n, 1n));
      return new Rec([["x", r.x], ["residual Ax − b", res], ["‖Ax − b‖²", sq]], "x minimises ‖Ax − b‖, and has the least norm among minimisers.");
    },
  },
  hermite: {
    args: "A",
    help: "Hermite normal form (integer A)",
    run: async ([A]) => {
      const r = L.hermite(needIntegerEntries(needM(A, "hermite"), "hermite"));
      return new Rec([["H", r.H], ["U", r.U]], "U·A = H with U unimodular (integer, determinant ±1).");
    },
  },
  smith: {
    args: "A",
    help: "Smith normal form (integer A)",
    run: async ([A]) => {
      const r = L.smith(needIntegerEntries(needM(A, "smith"), "smith"));
      return new Rec([["S", r.S], ["U", r.U], ["V", r.V], ["invariant factors", r.diagonal.map(toRational)]], "U·A·V = S, each diagonal entry dividing the next.");
    },
  },
  charpoly: { args: "A", help: "det(xI − A)", run: async ([A], ex) => new Poly(await ex(L.charpoly(needSquare(needM(A, "charpoly"), "charpoly")))) },
  minpoly: { args: "A", help: "minimal polynomial", run: async ([A], ex) => new Poly((await ex(L.similarity(needSquare(needM(A, "minpoly"), "minpoly")))).minimalPolynomial) },
  invariants: {
    args: "A",
    help: "invariant factors (similarity class)",
    run: async ([A], ex) => {
      const s = await ex(L.similarity(needSquare(needM(A, "invariants"), "invariants")));
      return new Rec(s.invariantFactors.map((p, i) => [`d${i + 1}`, new Poly(p)]), "Two matrices are similar exactly when these agree. Each divides the next; the last is the minimal polynomial.");
    },
  },
  rcf: { args: "A", help: "rational canonical form", run: async ([A], ex) => (await ex(L.similarity(needSquare(needM(A, "rcf"), "rcf")))).rationalCanonical },
  jordan: {
    args: "A",
    help: "Jordan form J and P with A = P J P⁻¹ (rational eigenvalues)",
    run: async ([A], ex) => {
      const s = await ex(L.similarity(needSquare(needM(A, "jordan"), "jordan")));
      if (!s.J) return new Note("Some eigenvalues are irrational or complex, so J is not a rational matrix. The rational canonical form rcf(A) is always rational; the Jordan page shows the exact irrational form.", "#jordan");
      const blocks = s.jordan.map((j) => `λ = ${j.value.toRational ? j.value.toRational() : j.value}: ${j.blocks.join(", ")}`).join(";  ");
      return new Rec([["J", s.J], ["P", s.P], ["block sizes", blocks]], "A = P·J·P⁻¹, checked exactly.");
    },
  },
  eigen: {
    args: "A",
    help: "rational eigenvalues and their eigenspaces",
    run: async ([A], ex) => {
      A = needSquare(needM(A, "eigen"), "eigen");
      const list = await ex(L.eigenvectors(A));
      const fields = [];
      let other = 0;
      for (const e of list) {
        const rational = e.value && e.value.isRational ? toRational(e.value.toRational()) : null;
        if (!rational) {
          other += e.multiplicity;
          continue;
        }
        fields.push([`λ = ${rational}`, new Rec([["multiplicity", String(e.multiplicity)], ["eigenspace", columns(e.vectors, A.rows, "{0}")]])]);
      }
      const note = other
        ? `${other} of the ${A.rows} eigenvalues (with multiplicity) are irrational or complex. They are exact algebraic numbers, not rationals, so they live on the eigenvalue page.`
        : "Every eigenvalue is rational.";
      return new Rec(fields.length ? fields : [["rational eigenvalues", "none"]], note);
    },
  },
  num: { args: "q", help: "numerator", run: async ([q]) => new Rational(needQ(q, "num").numerator, 1n) },
  den: { args: "q", help: "denominator", run: async ([q]) => new Rational(needQ(q, "den").denominator, 1n) },
  entry: {
    args: "A, i, j",
    help: "entry in row i, column j (from 1)",
    run: async ([A, i, j]) => {
      A = needM(A, "entry");
      const r = needInt(i, "entry");
      const c = needInt(j, "entry");
      if (r < 1 || r > A.rows || c < 1 || c > A.cols) fail(`entry(${r}, ${c}) is outside a ${A.rows}×${A.cols} matrix`);
      return A.get(r - 1, c - 1);
    },
  },
};

// ---------------------------------------------------------------------------
// Functions that need a number setting, and algebraic versions of rational ones

const algOn = () => CTX.domain === "rational" && CTX.algebraic;
const rangesOn = () => CTX.domain === "ranges";
const pairOf = (q) => [q.numerator, q.denominator];

function complexSqrtOfRational(q) {
  // the root of x² − q (q < 0) with positive imaginary part
  const roots = L.eigenvalues([[0, q], [1, 0]]).exact().complex;
  const z = roots.map((r) => r.value).find((v) => v && v.toNumbers().im > 0);
  return z ?? fail("could not certify that square root");
}

function eigenAlgebraic(A, ex) {
  return ex(L.eigenvectors(needSquare(needM(A, "eigen"), "eigen"))).then((list) => {
    const fields = [];
    let hidden = 0;
    for (const e of list) {
      const complex = isC(e.value) || (e.value === null);
      if (complex && !CTX.complex) {
        hidden += e.multiplicity;
        continue;
      }
      if (e.value === null) {
        fields.push([`λ (uncertified)`, new Note("This complex eigenvalue could not be certified.")]);
        continue;
      }
      const rational = !complex && e.value.isRational ? toRational(e.value.toRational()) : null;
      const vec = rational
        ? columns(e.vectors, A.rows, "{0}")
        : e.algebraic ? `(${e.algebraic.toStrings().join(",  ")})  with λ this eigenvalue` : "—";
      fields.push([`λ = ${rational ?? e.value.toString(14)}`, new Rec([
        ["value", rational ?? e.value],
        ["multiplicity", String(e.multiplicity)],
        [rational ? "eigenspace" : "eigenvector over Q(λ)", vec],
      ])]);
    }
    const note = hidden ? `${hidden} complex eigenvalue${hidden > 1 ? "s are" : " is"} hidden; turn on complex numbers to see ${hidden > 1 ? "them" : "it"}.` : "Every eigenvalue is shown exactly: a rational, or a polynomial with an interval holding just that root.";
    return new Rec(fields, note);
  });
}

Object.assign(FUNCTIONS, {
  sqrt: {
    args: "x",
    help: "square root (real algebraic; complex when allowed)",
    needs: algOn,
    needsLabel: "real algebraic numbers",
    run: async ([x]) => {
      if (isQ(x) && x.numerator < 0n) return CTX.complex ? complexSqrtOfRational(x) : fail("√ of a negative number is not real; turn on complex numbers");
      if (!(isQ(x) || isA(x))) fail(`sqrt needs a number, got a ${kind(x)}`);
      if (isA(x) && ALG.sign(x) < 0) return CTX.complex ? ALG.sqrt(x) : fail("√ of a negative number is not real; turn on complex numbers");
      try {
        return ALG.sqrt(FA.of(x));
      } catch (e) {
        return fail(e.message);
      }
    },
  },
  roots: {
    args: "p",
    help: "real (and complex) roots of a polynomial, e.g. roots(charpoly(A))",
    needs: algOn,
    needsLabel: "real algebraic numbers",
    run: async ([p]) => {
      if (!(p instanceof Poly)) fail("roots needs a polynomial, such as charpoly(A) or minpoly(A)");
      const n = p.c.length - 1;
      if (n < 1) return new Note("A constant has no roots.");
      const comp = Matrix.from(Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (j === n - 1 ? p.c[i].negate().divide(p.c[n]) : i === j + 1 ? 1 : 0))));
      const ev = await L.eigenvalues(comp).exact();
      const fields = ev.real.map((r) => [`x = ${r.value.isRational ? r.value.toRational() : r.value.toString(14)}`, r.multiplicity > 1 ? `multiplicity ${r.multiplicity}` : norm(r.value)]);
      if (CTX.complex) for (const c of ev.complex) fields.push([`x = ${c.value ? c.value.toString(12) : "?"}`, c.value ?? "uncertified"]);
      const nc = ev.complex.reduce((t, c) => t + (c.multiplicity ?? 1), 0);
      return new Rec(fields.length ? fields : [["real roots", "none"]], nc && !CTX.complex ? `${nc} complex root${nc > 1 ? "s" : ""} hidden; turn on complex numbers to see them.` : null);
    },
  },
  sv: {
    args: "A",
    help: "singular values, exactly",
    needs: algOn,
    needsLabel: "real algebraic numbers",
    run: async ([A], ex) => {
      const list = await ex(L.singularValues(needM(A, "sv")));
      return new Rec(list.map((s, i) => [`σ${i + 1}`, s.value.isRational ? toRational(s.value.toRational()) : s.value]), "Square roots of the eigenvalues of AᵀA, largest first.");
    },
  },
  digits: {
    args: "x, n",
    help: "n correct decimal digits of a real algebraic number",
    needs: algOn,
    needsLabel: "real algebraic numbers",
    run: async ([x, n]) => {
      const k = needInt(n, "digits");
      if (k < 1 || k > 2000) fail("digits runs from 1 to 2000");
      if (isQ(x)) return x.toDecimal(k);
      if (isA(x)) return x.toString(k);
      if (isC(x)) return x.toString(k);
      return fail(`digits needs a number, got a ${kind(x)}`);
    },
  },
  expm: {
    args: "A",
    help: "exp(A), every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.expm(needSquare(needM(A, "expm"), "expm"), { bits: 64 }))),
  },
  sinm: {
    args: "A",
    help: "sin(A), every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.sinm(needSquare(needM(A, "sinm"), "sinm"), { bits: 64 }))),
  },
  cosm: {
    args: "A",
    help: "cos(A), every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.cosm(needSquare(needM(A, "cosm"), "cosm"), { bits: 64 }))),
  },
  sinhm: {
    args: "A",
    help: "sinh(A), every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.sinhm(needSquare(needM(A, "sinhm"), "sinhm"), { bits: 64 }))),
  },
  coshm: {
    args: "A",
    help: "cosh(A), every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.coshm(needSquare(needM(A, "coshm"), "coshm"), { bits: 64 }))),
  },
  sqrtm: {
    args: "A",
    help: "the principal square root of A, every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.sqrtm(needSquare(needM(A, "sqrtm"), "sqrtm"), { bits: 64 }))),
  },
  logm: {
    args: "A",
    help: "the principal logarithm of A, every entry as a certified range",
    needs: rangesOn,
    needsLabel: "Ranges",
    run: async ([A], ex) => new FMatrix(FI, await ex(L.logm(needSquare(needM(A, "logm"), "logm"), { bits: 64 }))),
  },
  lo: { args: "x", help: "lower end of a range", needs: rangesOn, needsLabel: "Ranges", run: async ([x]) => (isI(x) ? x.low : needQ(x, "lo")) },
  hi: { args: "x", help: "upper end of a range", needs: rangesOn, needsLabel: "Ranges", run: async ([x]) => (isI(x) ? x.high : needQ(x, "hi")) },
  mid: { args: "x", help: "midpoint of a range", needs: rangesOn, needsLabel: "Ranges", run: async ([x]) => (isI(x) ? x.low.add(x.high).divide(new Rational(2n, 1n)) : needQ(x, "mid")) },
  width: { args: "x", help: "width of a range", needs: rangesOn, needsLabel: "Ranges", run: async ([x]) => (isI(x) ? x.high.subtract(x.low) : new Rational(0n, 1n)) },
});
FUNCTIONS.eigen.alg = ([A], ex) => eigenAlgebraic(A, ex);
FUNCTIONS.jordan.alg = async ([A], ex) => {
  const s = await ex(L.similarity(needSquare(needM(A, "jordan"), "jordan")));
  const blocks = new Rec(s.jordan
    .filter((j) => CTX.complex || !isC(j.value))
    .map((j) => [`λ = ${j.value?.isRational ? j.value.toRational() : j.value?.toString(12) ?? "?"}`, `blocks of size ${j.blocks.join(", ")}`]));
  if (s.J) return new Rec([["J", s.J], ["P", s.P], ["blocks", blocks]], "A = P·J·P⁻¹, checked exactly.");
  return new Rec([["Jordan blocks", blocks], ["rational canonical form", s.rationalCanonical]], `J has irrational${CTX.complex ? " or complex" : ""} entries, so it is given by its blocks: each block is λ on the diagonal and 1 just above.`);
};

const ALIASES = { inverse: "inv", transpose: "T", tr: "trace", nullspace: "null", kernel: "null", colspace: "col", rowspace: "row", eye: "I", gramSchmidt: "gs", pseudoinverse: "pinv", leastSquares: "lsq", minimal: "minpoly", eigenvalues: "eigen", eigenvectors: "eigen" };
export const resolveName = (name) => (FUNCTIONS[name] ? name : ALIASES[name]);

function checkSize(n) {
  if (n < 1 || n > 60) fail(`sizes run from 1 to 60, got ${n}`);
  return n;
}

// ---------------------------------------------------------------------------
// Evaluation

/**
 * Evaluate src against `vars` (name → value). Returns { name?, value }.
 * `exact(refinement)` settles a Refinement (pass one that honours a signal).
 */
const liftForCall = (v) => (CTX.domain === "modp" && isM(v) ? toF(FP(), v) : v);
function liftResult(v) {
  if (isQ(v)) return FP().wrap(FP().of(v));
  if (isM(v)) return toF(FP(), v);
  if (v instanceof Rec) return new Rec(v.fields.map(([k, x]) => [k, liftResult(x)]), v.note);
  return v;
}

export async function evaluate(src, vars, exact = async (r) => r.exact()) {
  const tree = parse(src);
  const ex = async (r) => (r instanceof L.Refinement ? exact(r) : r);
  async function ev(n) {
    switch (n.k) {
      case "num": return n.v;
      case "var": {
        if (n.name in vars) return vars[n.name];
        if (resolveName(n.name)) fail(`${n.name} is a function: write ${n.name}(…)`);
        return fail(`nothing is called ${n.name}`);
      }
      case "neg": {
        const a = await ev(n.a);
        return isQ(a) ? a.negate() : isM(a) ? a.scale(-1) : gNeg(a);
      }
      case "+": return add(await ev(n.a), await ev(n.b), 1);
      case "-": return add(await ev(n.a), await ev(n.b), -1);
      case "*": return mul(await ev(n.a), await ev(n.b));
      case "/": return div(await ev(n.a), await ev(n.b));
      case "T": {
        const a = await ev(n.a);
        return isM(a) ? a.transpose() : isF(a) ? gT(a) : a;
      }
      case "^": {
        const a = await ev(n.a);
        const e = needInt(await ev(n.b), "a power");
        if (!pureQ(a) || CTX.domain === "modp") {
          if (Math.abs(e) > 1e6) fail("that power is too large");
          return gPow(a, e);
        }
        if (isQ(a)) {
          if (a.numerator === 0n && e < 0) fail("0 has no negative powers");
          return a.pow(e);
        }
        const A = needSquare(needM(a, "a power"), "a matrix power");
        if (Math.abs(e) > 1e6) fail("that power is too large");
        const r = await ex(L.power(A, e));
        return r ?? new Note("The matrix is singular, so negative powers do not exist.");
      }
      case "mat": {
        const rows = [];
        for (const row of n.rows) {
          const parts = [];
          for (const x of row) parts.push(await ev(x));
          rows.push(hcat(parts));
        }
        return vcat(rows);
      }
      case "call": {
        const name = resolveName(n.name);
        if (!name) fail(`there is no function ${n.name}`);
        const f = FUNCTIONS[name];
        const want = f.args.split(",").length;
        const optional = name === "zeros" ? 1 : 0;
        if (n.args.length > want || n.args.length < want - optional) fail(`${name} takes ${want} argument${want > 1 ? "s" : ""}: ${name}(${f.args})`);
        const args = [];
        for (const a of n.args) args.push(liftForCall(await ev(a)));
        if (f.needs && !f.needs()) fail(`${name} needs the ${f.needsLabel} setting (see "Numbers" above)`);
        if (args.some(isF)) {
          if (!GENERIC[name]) fail(`${name} works on rational matrices; this one is over the ${args.find(isF).F.tag}`);
          const F = args.find(isF).F;
          return GENERIC[name](args.map((a) => (isM(a) ? toF(F, a) : a)), ex);
        }
        if (CTX.algebraic && f.alg) return f.alg(args, ex);
        return f.run(args, ex);
      }
      case "assign": return ev(n.e);
      default: return fail("can't evaluate that");
    }
  }
  let value = await ev(tree);
  if (CTX.domain === "modp") value = liftResult(value);
  value = norm(value);
  if (tree.k === "assign") {
    if (resolveName(tree.name)) fail(`${tree.name} is a function name; pick another name`);
    return { name: tree.name, value };
  }
  return { value };
}

const gNegOrQ = (a) => (isQ(a) ? a.negate() : gNeg(a));

/** Parse one matrix cell: a Core literal or a small scalar expression. */
export function parseCell(text) {
  const s = String(text).trim();
  if (s === "") return CTX.domain === "modp" ? FP().wrap(0n) : new Rational(0n, 1n);
  let simple = null;
  try {
    simple = parseNumber(s);
  } catch {
    simple = null;
  }
  if (simple !== null) {
    if (simple.low !== undefined) {
      if (CTX.domain !== "ranges") fail("ranges like 1:2 need the Ranges setting");
      return norm(simple);
    }
    const q = toRational(simple);
    return CTX.domain === "modp" ? FP().wrap(FP().of(q)) : q;
  }
  if (CTX.domain !== "rational" || CTX.algebraic) {
    const tree = parse(s);
    const walkG = (n) => {
      switch (n.k) {
        case "num": return n.v;
        case "neg": return gNegOrQ(walkG(n.a));
        case "+": return add(walkG(n.a), walkG(n.b), 1);
        case "-": return add(walkG(n.a), walkG(n.b), -1);
        case "*": return mul(walkG(n.a), walkG(n.b));
        case "/": return div(walkG(n.a), walkG(n.b));
        case "^": {
          const e = needInt(walkG(n.b), "a power");
          if (Math.abs(e) > 4096) fail("that power is too large");
          return gPow(walkG(n.a), e);
        }
        case "call": {
          if (n.name === "sqrt" && algOn() && n.args.length === 1) {
            const x = walkG(n.args[0]);
            if ((isQ(x) && x.numerator < 0n) || (isA(x) && ALG.sign(x) < 0)) fail("a matrix entry must be real");
            return norm(ALG.sqrt(FA.of(x)));
          }
          return fail(`a cell can't call ${n.name}(…)`);
        }
        default: return fail("a cell holds a number, like 3/4, 0.#3, 2..1/2, sqrt(2) or 1:2");
      }
    };
    let v = walkG(tree);
    if (isM(v) || isF(v)) fail("a cell holds one number");
    if (CTX.domain === "modp" && isQ(v)) v = FP().wrap(FP().of(v));
    return v;
  }
  const tree = parse(s);
  const walk = (n) => {
    switch (n.k) {
      case "num": return n.v;
      case "neg": return walk(n.a).negate();
      case "+": return walk(n.a).add(walk(n.b));
      case "-": return walk(n.a).subtract(walk(n.b));
      case "*": return walk(n.a).multiply(walk(n.b));
      case "/": {
        const b = walk(n.b);
        if (b.numerator === 0n) fail("division by zero");
        return walk(n.a).divide(b);
      }
      case "^": {
        const b = walk(n.b);
        if (b.denominator !== 1n) fail("powers must be whole numbers");
        const a = walk(n.a);
        if (a.numerator === 0n && b.numerator < 0n) fail("0 has no negative powers");
        if (b.numerator > 4096n || b.numerator < -4096n) fail("that power is too large");
        return a.pow(Number(b.numerator));
      }
      default: return fail("a cell holds a number, like 3/4, 0.#3, 2..1/2 or (1/2)^3");
    }
  };
  return walk(tree);
}

/** A matrix from parsed cells: a rational Matrix when it can be, else an FMatrix. */
export function fromCells(rows) {
  if (rows.flat().every(isQ) && CTX.domain !== "modp") return Matrix.from(rows);
  const F = CTX.domain === "modp" ? FP() : rows.flat().some(isI) ? FI : FA;
  return new FMatrix(F, rows.map((r) => r.map((x) => F.of(x))));
}

/** Text that parseCell reads back to the same value, or null (algebraic numbers). */
export function cellText(x) {
  if (isQ(x)) return x.toString();
  if (isI(x)) return `${x.low}:${x.high}`;
  if (isP(x)) return String(x.v);
  if (typeof x === "bigint") return String(x);
  return null;
}
export { isA, isI, isP, isF, isC };
