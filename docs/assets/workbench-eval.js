// The workbench's little language: exact rational scalars and matrices,
// arithmetic, and the library's exact operations as functions.
//
//   A*B - 2A'      [1, 1/2; 0.#3, 2..1/4]      A^-1      C = rref(A)
//
// Values are Core Rationals, Matrices, Poly (rational coefficients, constant
// first), Rec (labelled fields) and Note (a message). Everything is exact; the
// operations that return irrational numbers are not offered here.
import * as L from "./linalg.js";

const { Matrix, Rational, parseNumber } = L;

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
      if (v && v.low !== undefined) fail(`"${text}" is an interval; the workbench is exact only`);
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
const kind = (v) => (isM(v) ? `${v.rows}×${v.cols} matrix` : isQ(v) ? "number" : v instanceof Poly ? "polynomial" : "result");
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
  if (isQ(a) && isQ(b)) return sign > 0 ? a.add(b) : a.subtract(b);
  if (isM(a) && isM(b)) {
    if (a.rows !== b.rows || a.cols !== b.cols) fail(`can't ${sign > 0 ? "add" : "subtract"} a ${a.rows}×${a.cols} and a ${b.rows}×${b.cols} matrix`);
    return sign > 0 ? a.add(b) : a.subtract(b);
  }
  return fail(`can't ${sign > 0 ? "add" : "subtract"} a ${kind(a)} and a ${kind(b)} (for A + k·I write A + k*I(n))`);
}
function mul(a, b) {
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
  if (isQ(b) && b.numerator === 0n) fail("division by zero");
  if (isQ(a) && isQ(b)) return a.divide(b);
  if (isM(a) && isQ(b)) return a.scale(new Rational(1n, 1n).divide(b));
  return fail(`can't divide a ${kind(a)} by a ${kind(b)}; for matrices use inv(B) or solve(A, b)`);
}

function hcat(parts) {
  const ms = parts.map((p) => (isQ(p) ? Matrix.from([[p]]) : needM(p, "a matrix literal")));
  const rows = ms[0].rows;
  if (ms.some((m) => m.rows !== rows)) fail("entries side by side in [ … ] must have the same number of rows");
  return ms.reduce((x, y) => x.augment(y));
}
function vcat(ms) {
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
        return isQ(a) ? a.negate() : isM(a) ? a.scale(-1) : fail(`can't negate a ${kind(a)}`);
      }
      case "+": return add(await ev(n.a), await ev(n.b), 1);
      case "-": return add(await ev(n.a), await ev(n.b), -1);
      case "*": return mul(await ev(n.a), await ev(n.b));
      case "/": return div(await ev(n.a), await ev(n.b));
      case "T": {
        const a = await ev(n.a);
        return isM(a) ? a.transpose() : a;
      }
      case "^": {
        const a = await ev(n.a);
        const e = needInt(await ev(n.b), "a power");
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
        for (const a of n.args) args.push(await ev(a));
        return f.run(args, ex);
      }
      case "assign": return ev(n.e);
      default: return fail("can't evaluate that");
    }
  }
  const value = await ev(tree);
  if (tree.k === "assign") {
    if (resolveName(tree.name)) fail(`${tree.name} is a function name; pick another name`);
    return { name: tree.name, value };
  }
  return { value };
}

/** Parse one matrix cell: a Core literal or a small scalar expression. */
export function parseCell(text) {
  const s = String(text).trim();
  if (s === "") return new Rational(0n, 1n);
  try {
    const v = parseNumber(s);
    if (v && v.low !== undefined) fail("intervals are not exact entries");
    return toRational(v);
  } catch (e) {
    if (e instanceof EvalError) throw e;
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
