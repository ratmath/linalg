// The public operations. Each takes matrices (Matrix or nested arrays) and
// returns a Refinement: read `.quick` for the float preview, call `.exact()`
// or `await .settle()` for the exact answer.

import { Matrix } from "./matrix.js";
import { Refinement } from "./refine.js";
import { toRational } from "./q.js";
import {
  run,
  rrefPairs,
  determinantModular,
  determinantBareiss,
  solvePairs,
  inversePairs,
  provablyNonsingular,
} from "./exact.js";
import { rrefFloat, detFloat, solveFloat, inverseFloat, guessRational } from "./float.js";

const asMatrix = (m) => Matrix.from(m);

function asColumns(b, rows) {
  const m = Matrix.from(b);
  // a single row given for a column vector: [1, 2, 3] already reads as a column
  if (m.rows !== rows && m.rows === 1 && m.cols === rows) return m.transpose();
  return m;
}

/** Guess an exact matrix from floats; null if any entry has no small fraction near it. */
function guessMatrix(rows, maxDen) {
  const out = [];
  for (const row of rows) {
    const r = [];
    for (const v of row) {
      const g = guessRational(v, maxDen);
      if (!g) return null;
      r.push(g);
    }
    out.push(r);
  }
  return Matrix.fromPairs(out);
}

/**
 * Determinant. Float preview by LU; exact by the multi-modular method
 * (word-size primes + Chinese remaindering), or Bareiss for tiny matrices.
 */
export function det(input) {
  const A = asMatrix(input);
  if (!A.isSquare) throw new RangeError(`det needs a square matrix, got ${A.rows}×${A.cols}`);
  return new Refinement(
    () => detFloat(A.float(), A.rows),
    function* () {
      const small = A.rows <= 3;
      const pair = small ? yield* determinantBareiss(A.pairs()) : yield* determinantModular(A.pairs());
      return { value: toRational(pair), via: small ? "bareiss" : "modular" };
    },
  );
}

/** Rank. */
export function rank(input) {
  const A = asMatrix(input);
  return new Refinement(
    () => {
      const { pivots, unsure } = rrefFloat(A.float(), A.rows, A.cols);
      return { value: pivots.length, trust: unsure ? "unsure" : "good" };
    },
    function* () {
      const { pivots } = yield* rrefPairs(A.pairs(), A.cols);
      return { value: pivots.length, via: "bareiss" };
    },
  );
}

/** Reduced row echelon form: { matrix, pivots }. */
export function rref(input, { columns } = {}) {
  const A = asMatrix(input);
  const cols = columns ?? A.cols;
  return new Refinement(
    () => {
      const { rows, pivots, unsure } = rrefFloat(A.float(), A.rows, A.cols, cols);
      return { value: { matrix: rows, pivots }, trust: unsure ? "unsure" : "good" };
    },
    function* () {
      const { rows, pivots } = yield* rrefPairs(A.pairs(), cols);
      return { value: { matrix: Matrix.fromPairs(rows), pivots }, via: "bareiss" };
    },
  );
}

/**
 * Solve A x = b. b may be a vector or a matrix of right-hand sides.
 * Value: { status: "unique" | "infinite" | "none", x, basis }
 *   x      one solution (free unknowns 0), or null when there is none
 *   basis  null-space vectors of A; x + any combination of them also solves
 * Float stage: x is number[][] and basis number[][]; exact: Matrix and Matrix[].
 *
 * Exact path: when the preview looks unique, guess small fractions from the
 * floats and check A·x = b exactly, with a modular proof that A is
 * nonsingular. That check costs O(n²) exact operations, so nice answers
 * arrive almost at once; otherwise fall back to fraction-free elimination.
 */
export function solve(aInput, bInput, { guessDenominator = 1e6 } = {}) {
  const A = asMatrix(aInput);
  const B = asColumns(bInput, A.rows);
  if (B.rows !== A.rows) throw new RangeError(`solve: A has ${A.rows} rows, b has ${B.rows}`);
  const job = new Refinement(
    () => solveFloat(A.float(), A.rows, A.cols, B.toNumbers()),
    function* () {
      const q = job.quick;
      if (A.isSquare && q.value.status === "unique" && q.trust !== "unsure") {
        const X = guessMatrix(q.value.x, guessDenominator);
        yield { step: "check", done: 0, of: 1 };
        if (X && A.multiply(X).equals(B) && provablyNonsingular(A.pairs())) {
          return { value: { status: "unique", x: X, basis: [] }, via: "check" };
        }
      }
      const { status, x, basis } = yield* solvePairs(A.pairs(), B.pairs());
      return {
        value: {
          status,
          x: x && Matrix.fromPairs(x),
          basis: basis.map((v) => Matrix.fromPairs(v.map((p) => [p]))),
        },
        via: "bareiss",
      };
    },
  );
  return job;
}

/** Inverse, or null when singular. Exact path guesses from the float inverse and checks A·X = I. */
export function inverse(input, { guessDenominator = 1e6 } = {}) {
  const A = asMatrix(input);
  if (!A.isSquare) throw new RangeError(`inverse needs a square matrix, got ${A.rows}×${A.cols}`);
  const job = new Refinement(
    () => {
      const { value, trust } = inverseFloat(A.float(), A.rows);
      return { value, trust };
    },
    function* () {
      const q = job.quick;
      if (q.value && q.trust !== "unsure") {
        const X = guessMatrix(q.value, guessDenominator);
        yield { step: "check", done: 0, of: 1 };
        if (X && A.multiply(X).equals(Matrix.identity(A.rows))) return { value: X, via: "check" };
      }
      const inv = yield* inversePairs(A.pairs());
      return { value: inv && Matrix.fromPairs(inv), via: "bareiss" };
    },
  );
  return job;
}

/** Null-space basis of A (columns), as a Refinement. */
export function nullspace(input) {
  const A = asMatrix(input);
  const zero = Matrix.zeros(A.rows, 1);
  const job = solve(A, zero);
  return new Refinement(
    () => ({ value: job.quick.value.basis, trust: job.quick.trust }),
    function* () {
      const gen = job.stages();
      let last;
      for (const s of gen) {
        if (s.stage === "progress") yield s;
        last = s;
      }
      return { value: last.value.basis, via: last.via };
    },
  );
}

/** Immediate exact product (multiplication needs no refinement). */
export const multiply = (a, b) => asMatrix(a).multiply(b);

export { run };
