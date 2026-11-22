// Teaching forms. These are deliberately not the fast routines: they use the
// textbook algorithm, keep every intermediate state, and say in words what
// each step does and why. Use them to explain; use the main API to compute.
//
// Every walkthrough returns { title, summary, steps, result }, where each step
// is plain data a page can render:
//   { title, text, matrix?: string[][], bar?: column index of the augmentation
//     bar, pivot?: [row, col], rows?: [highlighted rows], table?: { head, rows },
//     note? }

import { Matrix } from "../matrix.js";
import { toPair, reduce, addQ, subQ, mulQ, divQ, abs, pairToNumber, lcm, gcd } from "../q.js";
import { run, determinantBareiss, integerRows, detMod, hadamardBits, prime } from "../exact.js";
import { guessRational } from "../float.js";
import { solveFloat } from "../float.js";
import { charpolyPairs, fromPairs, squareFree, isolate, derivative, signAt, degree } from "../poly.js";
import { formatP } from "../qpoly.js";

export { RowReducer } from "./row-reducer.js";

const Z = [0n, 1n];
const ONE = [1n, 1n];
const isZ = ([n]) => n === 0n;
export const fmt = ([n, d]) => (d === 1n ? `${n}` : `${n}/${d}`);
const fmtSigned = (p) => (p[0] < 0n ? `− ${fmt([-p[0], p[1]])}` : `+ ${fmt(p)}`);
const snap = (M) => M.map((row) => row.map(fmt));
const sub = (k) => String(k).split("").map((c) => "₀₁₂₃₄₅₆₇₈₉"[c]).join("");
const R = (i) => `R${sub(i + 1)}`;
const x = (j) => `x${sub(j + 1)}`;
const pairsOf = (A) => Matrix.from(A).pairs();
const neg = ([n, d]) => [-n, d];
const coefTimes = (c, what) => (c[0] === 1n && c[1] === 1n ? what : `${fmt(c)}·${what}`);

/** The ordinary Gauss–Jordan elimination, one elementary row operation per step. */
export function gaussJordan(input, { bar = null, columns } = {}) {
  const M = pairsOf(input).map((r) => [...r]);
  const m = M.length;
  const n = M[0]?.length ?? 0;
  const cols = columns ?? (bar ?? n);
  const steps = [{
    title: "Start",
    text: bar === null
      ? "We reduce the matrix with three kinds of row operation: swap two rows, multiply a row by a nonzero number, and add a multiple of one row to another. None of them changes the solutions of the system the rows describe."
      : "The matrix is augmented: everything right of the bar is carried along. Row operations on the whole augmented matrix never change the solutions.",
    matrix: snap(M),
    bar,
  }];
  const pivots = [];
  let r = 0;
  for (let col = 0; col < cols && r < m; col += 1) {
    let p = r;
    while (p < m && isZ(M[p][col])) p += 1;
    if (p === m) {
      steps.push({
        title: `Column ${col + 1}: no pivot`,
        text: `Every entry of column ${col + 1} from row ${r + 1} down is 0, so there is nothing to pivot on. ${bar !== null ? `The unknown ${x(col)} will be free.` : "This column will not hold a leading 1."}`,
        matrix: snap(M), bar, rows: [], note: "free",
      });
      continue;
    }
    if (p !== r) {
      [M[p], M[r]] = [M[r], M[p]];
      steps.push({
        title: `Swap ${R(r)} ↔ ${R(p)}`,
        text: `Row ${r + 1} has 0 in column ${col + 1}, so bring up row ${p + 1}, which does not.`,
        matrix: snap(M), bar, pivot: [r, col], rows: [r, p], op: { type: "swap", i: r, j: p },
      });
    }
    const piv = M[r][col];
    if (!(piv[0] === 1n && piv[1] === 1n)) {
      const c = divQ(ONE, piv);
      M[r] = M[r].map((v) => mulQ(v, c));
      steps.push({
        title: `${R(r)} ← ${fmt(c)}·${R(r)}`,
        text: `Divide row ${r + 1} by the pivot ${fmt(piv)} so the pivot becomes 1.`,
        matrix: snap(M), bar, pivot: [r, col], rows: [r], op: { type: "scale", i: r, c: fmt(c) },
      });
    }
    for (let i = 0; i < m; i += 1) {
      if (i === r || isZ(M[i][col])) continue;
      const f = M[i][col];
      M[i] = M[i].map((v, j) => subQ(v, mulQ(f, M[r][j])));
      steps.push({
        title: `${R(i)} ← ${R(i)} ${f[0] < 0n ? "+" : "−"} ${coefTimes([abs(f[0]), f[1]], R(r))}`,
        text: `Clear the ${fmt(f)} in row ${i + 1}, column ${col + 1}, using the pivot row. ${i < r ? "Clearing above the pivot too is what makes the form reduced." : ""}`,
        matrix: snap(M), bar, pivot: [r, col], rows: [i], op: { type: "add", i, j: r, c: fmt(neg(f)) },
      });
    }
    pivots.push(col);
    r += 1;
  }
  steps.push({
    title: "Reduced row echelon form",
    text: `Each nonzero row starts with a 1 (a pivot), each pivot is the only nonzero entry in its column, and pivots move right as we go down. Pivot columns: ${pivots.map((c) => c + 1).join(", ") || "none"}. Rank ${pivots.length}.`,
    matrix: snap(M), bar,
  });
  return { title: "Gauss–Jordan elimination", summary: `rank ${pivots.length}`, steps, result: { matrix: snap(M), pivots, rank: pivots.length, pairs: M } };
}

/** Solve A·x = b, explaining how the reduced form is read. */
export function solveSystem(aInput, bInput) {
  const A = pairsOf(aInput);
  let B = Matrix.from(bInput);
  if (B.rows !== A.length && B.rows === 1) B = B.transpose();
  const b = B.pairs();
  const n = A[0].length;
  const wide = A.map((row, i) => [...row, ...b[i]]);
  const walk = gaussJordan(Matrix.fromPairs(wide), { bar: n });
  const { pairs: M, pivots } = walk.result;
  const steps = walk.steps;
  const bad = M.findIndex((row) => row.slice(0, n).every(isZ) && row.slice(n).some((v) => !isZ(v)));
  let result;
  if (bad >= 0) {
    steps.push({
      title: "No solution",
      text: `Row ${bad + 1} now reads 0 = ${fmt(M[bad][n])}. That is impossible, so the system is inconsistent: no x satisfies every equation.`,
      matrix: snap(M), bar: n, rows: [bad],
    });
    result = { status: "none" };
  } else {
    const free = Array.from({ length: n }, (_, j) => j).filter((j) => !pivots.includes(j));
    const lines = pivots.map((col, r) => {
      const terms = free.filter((f) => !isZ(M[r][f])).map((f) => `${fmtSigned(neg(M[r][f]))}·${x(f)}`);
      return `${x(col)} = ${fmt(M[r][n])}${terms.length ? " " + terms.join(" ") : ""}`;
    });
    free.forEach((f) => lines.push(`${x(f)} is free`));
    steps.push({
      title: free.length ? "Infinitely many solutions" : "The unique solution",
      text: free.length
        ? `Each pivot row expresses one pivot unknown in terms of the free ones. Choose any values for ${free.map(x).join(", ")}: every choice gives a solution.`
        : "Every unknown has a pivot, so each row reads off one unknown directly.",
      matrix: snap(M), bar: n, lines,
    });
    result = { status: free.length ? "infinite" : "unique", lines, free };
  }
  return { title: "Solving a linear system", summary: result.status, steps, result };
}

/** Inverse via [A | I]. */
export function inverse(input) {
  const A = pairsOf(input);
  const n = A.length;
  const wide = A.map((row, i) => [...row, ...Array.from({ length: n }, (_, j) => (i === j ? ONE : Z))]);
  const walk = gaussJordan(Matrix.fromPairs(wide), { bar: n });
  const { pairs: M, pivots } = walk.result;
  const ok = pivots.length === n;
  walk.steps.push(ok
    ? { title: "Read off the inverse", text: "The left half became the identity, so the row operations, applied to I on the right, built A⁻¹ there.", matrix: snap(M.map((r) => r.slice(n))) }
    : { title: "Not invertible", text: `Only ${pivots.length} pivots for ${n} columns: the left half cannot become I. A is singular.`, matrix: snap(M), bar: n });
  return { title: "Inverse by Gauss–Jordan", summary: ok ? "invertible" : "singular", steps: walk.steps, result: ok ? snap(M.map((r) => r.slice(n))) : null };
}

/** Determinant by elimination, tracking how each operation changes it. */
export function determinantByElimination(input) {
  const M = pairsOf(input).map((r) => [...r]);
  const n = M.length;
  let sign = 1;
  const steps = [{ title: "Start", text: "Reduce to upper triangular form. A swap flips the sign of the determinant; adding a multiple of one row to another leaves it unchanged. Then the determinant is the product of the diagonal.", matrix: snap(M) }];
  for (let col = 0; col < n; col += 1) {
    let p = col;
    while (p < n && isZ(M[p][col])) p += 1;
    if (p === n) {
      steps.push({ title: `Column ${col + 1}: no pivot`, text: "A zero column below the diagonal means the rows are dependent: the determinant is 0.", matrix: snap(M) });
      return { title: "Determinant by elimination", summary: "0", steps, result: "0" };
    }
    if (p !== col) {
      [M[p], M[col]] = [M[col], M[p]];
      sign = -sign;
      steps.push({ title: `Swap ${R(col)} ↔ ${R(p)}`, text: `A swap multiplies the determinant by −1 (sign now ${sign > 0 ? "+" : "−"}).`, matrix: snap(M), rows: [col, p], pivot: [col, col] });
    }
    for (let i = col + 1; i < n; i += 1) {
      if (isZ(M[i][col])) continue;
      const f = divQ(M[i][col], M[col][col]);
      M[i] = M[i].map((v, j) => subQ(v, mulQ(f, M[col][j])));
      steps.push({ title: `${R(i)} ← ${R(i)} ${f[0] < 0n ? "+" : "−"} ${coefTimes([abs(f[0]), f[1]], R(col))}`, text: "Adding a multiple of another row does not change the determinant.", matrix: snap(M), rows: [i], pivot: [col, col] });
    }
  }
  let d = sign > 0 ? ONE : [-1n, 1n];
  for (let i = 0; i < n; i += 1) d = mulQ(d, M[i][i]);
  steps.push({ title: "Multiply the diagonal", text: `${sign < 0 ? "−" : ""}(${M.map((r, i) => fmt(r[i])).join(")·(")}) = ${fmt(d)}`, matrix: snap(M) });
  return { title: "Determinant by elimination", summary: fmt(d), steps, result: fmt(d) };
}

/** Cofactor (Laplace) expansion along the first row, for small matrices. */
export function cofactorExpansion(input) {
  const A = pairsOf(input);
  const n = A.length;
  if (n > 5) throw new RangeError("cofactor expansion is shown for at most 5×5");
  const steps = [{ title: "Start", text: "Expand along the first row: det A = Σ (−1)^(1+j) · a₁ⱼ · det(M₁ⱼ), where M₁ⱼ deletes row 1 and column j.", matrix: snap(A) }];
  let total = Z;
  const table = [];
  for (let j = 0; j < n; j += 1) {
    const minor = A.slice(1).map((r) => r.filter((_, k) => k !== j));
    const dm = n === 1 ? ONE : run(determinantBareiss(minor));
    const term = mulQ(mulQ(A[0][j], dm), j % 2 ? [-1n, 1n] : ONE);
    total = addQ(total, term);
    table.push([`${j + 1}`, j % 2 ? "−" : "+", fmt(A[0][j]), fmt(dm), fmt(term)]);
    steps.push({
      title: `Term ${j + 1}`,
      text: `Delete row 1 and column ${j + 1}. The minor has determinant ${fmt(dm)}; with sign ${j % 2 ? "−" : "+"} and entry ${fmt(A[0][j])} the term is ${fmt(term)}.${isZ(A[0][j]) ? " (A zero entry makes the whole term 0, which is why one expands along rows with many zeros.)" : ""}`,
      matrix: snap(A), rows: [0], cols: [j], minor: snap(minor),
    });
  }
  steps.push({ title: "Add the terms", text: `det A = ${fmt(total)}`, table: { head: ["j", "sign", "a₁ⱼ", "det M₁ⱼ", "term"], rows: table } });
  return { title: "Cofactor expansion", summary: fmt(total), steps, result: fmt(total) };
}

/** Fraction-free (Bareiss) elimination: integers throughout, exact divisions. */
export function bareiss(input) {
  const { rows: ints, scales } = integerRows(pairsOf(input));
  const M = ints.map((r) => [...r]);
  const n = M.length;
  const s = (v) => `${v}`;
  const steps = [];
  if (scales.some((c) => c !== 1n)) {
    steps.push({ title: "Clear denominators", text: `Multiply each row by the least common multiple of its denominators (${scales.join(", ")}). The determinant gets multiplied by their product, ${scales.reduce((a, b) => a * b, 1n)}, which we divide out at the end.`, matrix: M.map((r) => r.map(s)) });
  } else {
    steps.push({ title: "Start", text: "All entries are integers. Bareiss elimination keeps them integers: each update is (p·a − b·c) / previous pivot, and that division is always exact.", matrix: M.map((r) => r.map(s)) });
  }
  let prev = 1n;
  let sign = 1n;
  for (let k = 0; k < n; k += 1) {
    let p = k;
    while (p < n && M[p][k] === 0n) p += 1;
    if (p === n) {
      steps.push({ title: "Zero column", text: "No nonzero pivot: the determinant is 0.", matrix: M.map((r) => r.map(s)) });
      return { title: "Bareiss elimination", summary: "0", steps, result: "0" };
    }
    if (p !== k) {
      [M[p], M[k]] = [M[k], M[p]];
      sign = -sign;
      steps.push({ title: `Swap ${R(k)} ↔ ${R(p)}`, text: "Swapping flips the sign.", matrix: M.map((r) => r.map(s)), rows: [k, p] });
    }
    const piv = M[k][k];
    for (let i = k + 1; i < n; i += 1) {
      for (let j = k + 1; j < n; j += 1) M[i][j] = (piv * M[i][j] - M[i][k] * M[k][j]) / prev;
      M[i][k] = 0n;
    }
    steps.push({
      title: `Pivot ${k + 1}: ${piv}`,
      text: `For every entry below and right of the pivot: a ← (${piv}·a − (its entry in the pivot column)·(its entry in the pivot row)) / ${prev}, the previous pivot. Each new entry is a ${k + 2}×${k + 2} minor of the original matrix (Sylvester's identity), so the division leaves no remainder and the numbers stay as small as they can.`,
      matrix: M.map((r) => r.map(s)), pivot: [k, k],
    });
    prev = piv;
  }
  const d = reduce(sign * M[n - 1][n - 1], scales.reduce((a, b) => a * b, 1n));
  steps.push({ title: "Determinant", text: `The last pivot ${M[n - 1][n - 1]} is the determinant of the scaled matrix${sign < 0n ? " (with the sign flipped for the swaps)" : ""}; dividing by the row scales gives ${fmt(d)}.` });
  return { title: "Bareiss elimination", summary: fmt(d), steps, result: fmt(d) };
}

/** Determinant modulo small primes, glued with the Chinese remainder theorem. */
export function modularDeterminant(input, { start = 101 } = {}) {
  const { rows, scales } = integerRows(pairsOf(input));
  const bits = hadamardBits(rows);
  const bound = 2n ** BigInt(bits);
  const steps = [{
    title: "How big can the answer be?",
    text: `Hadamard's inequality: |det| ≤ product of the row lengths ≤ 2^${bits}. So if we know det modulo a number bigger than 2·2^${bits}, we know det exactly. Instead of one big computation we do several small ones, one per prime, where every number stays below the prime.`,
    matrix: rows.map((r) => r.map(String)),
  }];
  const primes = [];
  for (let p = start; primes.reduce((a, q) => a * BigInt(q), 1n) <= 2n * bound; p += 1) {
    let isP = p > 1;
    for (let q = 2; q * q <= p; q += 1) if (p % q === 0) { isP = false; break; }
    if (isP) primes.push(p);
  }
  let x = 0n;
  let M = 1n;
  const table = [];
  for (const p of primes) {
    const r = detMod(rows, p);
    const P = BigInt(p);
    // solve x + M·t ≡ r (mod p)
    let inv = 1n;
    const mm = ((M % P) + P) % P;
    for (let t = 1n; t < P; t += 1n) if ((mm * t) % P === 1n) { inv = t; break; }
    const t = ((((BigInt(r) - x) % P) + P) % P) * inv % P;
    x += M * t;
    M *= P;
    table.push([`${p}`, `${r}`, `${x}`, `${M}`]);
    steps.push({
      title: `Modulo ${p}`,
      text: `Eliminate with arithmetic mod ${p}, where every number stays below ${p}: det ≡ ${r} (mod ${p}). The Chinese remainder theorem combines this with the earlier residues: ${x} is the only number below ${M} that matches all of them.`,
      table: { head: ["prime", "det mod p", "combined", "modulus"], rows: [...table] },
    });
  }
  const signed = x > M / 2n ? x - M : x;
  const d = reduce(signed, scales.reduce((a, b) => a * b, 1n));
  steps.push({ title: "Symmetric representative", text: `The modulus ${M} exceeds twice the bound, so the true determinant is the representative nearest 0: ${signed}${scales.some((c) => c !== 1n) ? `, divided by the row scales: ${fmt(d)}` : ""}.` });
  return { title: "Determinant by the Chinese remainder theorem", summary: fmt(d), steps, result: fmt(d) };
}

/** Continued-fraction guess from floats, then an exact check. */
export function guessAndCheck(aInput, bInput) {
  const A = Matrix.from(aInput);
  let B = Matrix.from(bInput);
  if (B.rows !== A.rows && B.rows === 1) B = B.transpose();
  const fl = solveFloat(A.float(), A.rows, A.cols, B.toNumbers()).value;
  const steps = [{ title: "Float solve", text: "Solve in ordinary floating point first. It is fast but each value carries rounding error.", table: { head: ["unknown", "float value"], rows: (fl.x ?? []).map((r, i) => [x(i), `${r[0]}`]) } }];
  if (fl.status !== "unique") {
    steps.push({ title: "No unique float solution", text: "The guess needs a unique solution; the exact solver would use elimination instead." });
    return { title: "Guess and check", summary: "no guess", steps, result: null };
  }
  const guesses = fl.x.map((r) => guessRational(r[0], 1e6));
  steps.push({
    title: "Nearest simple fractions",
    text: "For each value, run the continued-fraction expansion and stop at the last convergent with denominator up to 10⁶. If the true answer has small denominators, this recovers it exactly.",
    table: { head: ["unknown", "float", "guess"], rows: fl.x.map((r, i) => [x(i), `${r[0]}`, guesses[i] ? fmt(guesses[i]) : "—"]) },
  });
  if (guesses.some((g) => !g)) return { title: "Guess and check", summary: "no guess", steps, result: null };
  const X = Matrix.fromPairs(guesses.map((g) => [g]));
  const AX = A.multiply(X);
  const ok = AX.equals(B);
  steps.push({
    title: ok ? "Check passes" : "Check fails",
    text: ok
      ? "Multiply A by the guessed fractions in exact arithmetic: it gives b exactly. A guess is never trusted on its own; this check is the proof."
      : "Exactly, A·x ≠ b, so the answer does not have small denominators. The library falls back to exact elimination.",
    table: { head: ["row", "A·x (exact)", "b"], rows: AX.toStrings().map((r, i) => [`${i + 1}`, r[0], B.toStrings()[i][0]]) },
  });
  return { title: "Guess and check", summary: ok ? "exact" : "fallback", steps, result: ok ? X.toStrings().map((r) => r[0]) : null };
}

/** Gram–Schmidt orthogonalisation without square roots. */
export function gramSchmidt(input) {
  const A = pairsOf(input);
  const m = A.length;
  const n = A[0].length;
  const col = (j) => A.map((r) => r[j]);
  const dot = (u, v) => u.reduce((s, a, i) => addQ(s, mulQ(a, v[i])), Z);
  const qs = [];
  const steps = [{ title: "Start", text: "Make the columns perpendicular one at a time: from each column subtract its projection onto every earlier result. Projections use (a·q)/(q·q), which needs no square roots, so everything stays exact.", matrix: snap(A) }];
  for (let j = 0; j < n; j += 1) {
    let v = col(j);
    const parts = [];
    for (const [k, q] of qs.entries()) {
      const qq = dot(q, q);
      if (isZ(qq)) continue;
      const c = divQ(dot(col(j), q), qq);
      parts.push(`(${fmt(c)})·q${sub(k + 1)}`);
      v = v.map((a, i) => subQ(a, mulQ(c, q[i])));
    }
    qs.push(v);
    steps.push({
      title: `q${sub(j + 1)}`,
      text: `q${sub(j + 1)} = a${sub(j + 1)}${parts.length ? " − " + parts.join(" − ") : ""}${v.every(isZ) ? ". It is the zero vector: column " + (j + 1) + " depends on the earlier ones." : ""}`,
      matrix: snap(Array.from({ length: m }, (_, i) => qs.map((q) => q[i]))),
      cols: [j],
    });
  }
  const norms = qs.map((q) => fmt(dot(q, q)));
  steps.push({ title: "Squared lengths", text: `‖q‖² = ${norms.join(", ")}. Dividing each q by √‖q‖² would make them orthonormal; that is the only step that leaves the rationals.` });
  return { title: "Gram–Schmidt", summary: `${qs.filter((q) => !q.every(isZ)).length} independent`, steps, result: { Q: snap(Array.from({ length: m }, (_, i) => qs.map((q) => q[i]))), normsSquared: norms } };
}

/** LU decomposition: the multipliers of elimination fill L. */
export function lu(input) {
  const U = pairsOf(input).map((r) => [...r]);
  const n = U.length;
  const L = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? ONE : Z)));
  const perm = Array.from({ length: n }, (_, i) => i);
  const steps = [{ title: "Start", text: "Eliminate below each pivot, and record each multiplier in L at the place it cleared. At the end PA = LU, with P recording any row swaps.", matrix: snap(U), L: snap(L) }];
  for (let k = 0; k < Math.min(n, U[0].length); k += 1) {
    let p = k;
    while (p < n && isZ(U[p][k])) p += 1;
    if (p === n) continue;
    if (p !== k) {
      [U[p], U[k]] = [U[k], U[p]];
      [perm[p], perm[k]] = [perm[k], perm[p]];
      for (let j = 0; j < k; j += 1) [L[p][j], L[k][j]] = [L[k][j], L[p][j]];
      steps.push({ title: `Swap ${R(k)} ↔ ${R(p)}`, text: "Zero pivot, so swap rows (and the multipliers already recorded in L).", matrix: snap(U), L: snap(L), rows: [k, p] });
    }
    for (let i = k + 1; i < n; i += 1) {
      if (isZ(U[i][k])) continue;
      const f = divQ(U[i][k], U[k][k]);
      L[i][k] = f;
      U[i] = U[i].map((v, j) => subQ(v, mulQ(f, U[k][j])));
      steps.push({ title: `ℓ${sub(i + 1)}${sub(k + 1)} = ${fmt(f)}`, text: `${R(i)} ← ${R(i)} − ${fmt(f)}·${R(k)}; the multiplier ${fmt(f)} goes into L.`, matrix: snap(U), L: snap(L), rows: [i], pivot: [k, k] });
    }
  }
  steps.push({ title: "PA = LU", text: `Row order ${perm.map((p) => p + 1).join(", ")}.`, matrix: snap(U), L: snap(L) });
  return { title: "LU decomposition", summary: `P = (${perm.map((p) => p + 1).join(" ")})`, steps, result: { L: snap(L), U: snap(U), perm } };
}

/** Least squares by the normal equations. */
export function leastSquares(aInput, bInput) {
  const A = Matrix.from(aInput);
  let B = Matrix.from(bInput);
  if (B.rows !== A.rows && B.rows === 1) B = B.transpose();
  const At = A.transpose();
  const AtA = At.multiply(A);
  const Atb = At.multiply(B);
  const steps = [
    { title: "The problem", text: "There are more equations than unknowns, so usually no x solves A·x = b. Look instead for the x that makes the error b − A·x as short as possible.", matrix: A.toStrings() },
    { title: "Normal equations", text: "The shortest error is perpendicular to every column of A: Aᵀ(b − A·x) = 0, so AᵀA·x = Aᵀb. AᵀA is square.", matrix: AtA.augment(Atb).toStrings(), bar: AtA.cols },
  ];
  const walk = solveSystem(AtA, Atb);
  steps.push(...walk.steps.slice(1).map((s) => ({ ...s, title: `Solve: ${s.title}` })));
  return { title: "Least squares", summary: walk.result.status, steps, result: walk.result };
}

/** Eigenvalues the way a course does it, then how the library certifies them. */
export function eigen(input) {
  const A = pairsOf(input);
  const n = A.length;
  const cp = run(charpolyPairs(A));
  const steps = [
    { title: "Characteristic polynomial", text: `λ is an eigenvalue when A − λI is singular, i.e. det(λI − A) = 0. Here det(λI − A) = ${formatP(cp, "λ")}.`, matrix: snap(A) },
  ];
  const sf = squareFree(fromPairs(cp));
  steps.push({
    title: "Repeated factors",
    text: `Split off repeated roots with gcd(p, p′): ${sf.map((f) => `(${formatP(f.factor.map((c) => [c, 1n]), "λ")})${f.multiplicity > 1 ? "^" + f.multiplicity : ""}`).join(" · ")}. Each factor now has only simple roots.`,
  });
  const found = [];
  for (const f of sf) {
    const sturm = sturmWalk(f.factor);
    steps.push(...sturm.steps);
    for (const r of sturm.result) found.push({ value: r, multiplicity: f.multiplicity });
  }
  for (const { value, multiplicity } of found) {
    if (!value.isRational) continue;
    const lam = value.rational;
    const shifted = A.map((row, i) => row.map((v, j) => (i === j ? subQ(v, lam) : v)));
    const w = gaussJordan(Matrix.fromPairs(shifted));
    const free = Array.from({ length: n }, (_, j) => j).filter((j) => !w.result.pivots.includes(j));
    steps.push({
      title: `Eigenvectors for λ = ${fmt(lam)}`,
      text: `Row reduce A − ${fmt(lam)}I. ${free.length} free column${free.length === 1 ? "" : "s"}, so the eigenspace has dimension ${free.length} (algebraic multiplicity ${multiplicity}).${free.length < multiplicity ? " Fewer eigenvectors than the multiplicity: A is not diagonalizable." : ""}`,
      matrix: w.result.matrix,
    });
  }
  return { title: "Eigenvalues", summary: found.map((f) => f.value.toString(8)).join(", "), steps, result: found };
}

/** Sturm sequence and bisection for the real roots of a square-free integer polynomial. */
export function sturmWalk(p) {
  const P = p.map((c) => [c, 1n]);
  const chain = [p, derivative(p)];
  while (degree(chain[chain.length - 1]) > 0) {
    const a = chain[chain.length - 2].map((c) => [c, 1n]);
    const b = chain[chain.length - 1].map((c) => [c, 1n]);
    // remainder over Q, negated
    let r = [...a];
    const db = b.length - 1;
    while (r.length - 1 >= db && r.length) {
      const k = r.length - 1 - db;
      const c = divQ(r[r.length - 1], b[db]);
      for (let i = 0; i <= db; i += 1) r[k + i] = subQ(r[k + i], mulQ(c, b[i]));
      while (r.length && isZ(r[r.length - 1])) r.pop();
    }
    if (!r.length) break;
    chain.push(fromPairs(r.map(neg)).map((c) => c)); // primitive keeps sign? recompute sign below
    const last = chain[chain.length - 1];
    const lead = r[r.length - 1];
    if ((last[last.length - 1] > 0n) !== (-lead[0] > 0n)) chain[chain.length - 1] = last.map((c) => -c);
  }
  const V = (pt) => {
    let v = 0;
    let last = 0;
    for (const q of chain) {
      const s = signAt(q, pt);
      if (s === 0) continue;
      if (last && s !== last) v += 1;
      last = s;
    }
    return v;
  };
  const roots = isolate(p);
  const steps = [{
    title: `Sturm sequence for ${formatP(P, "λ")}`,
    text: "Sturm's theorem: start with p and p′, keep taking negated remainders. The number of real roots in (a, b] is V(a) − V(b), where V(x) counts sign changes along the sequence at x. Counting is exact, so we can bisect until each interval holds one root.",
    lines: chain.map((q, i) => `S${sub(i)} = ${formatP(q.map((c) => [c, 1n]), "λ")}`),
  }];
  const table = roots.map((r) => (r.isRational ? [fmt(r.rational), fmt(r.rational), "exact"] : [fmt(r.lo), fmt(r.hi), `${V(r.lo) - V(r.hi)} root`]));
  steps.push({ title: "Isolating intervals", text: `${roots.length} real root${roots.length === 1 ? "" : "s"}. Each interval below holds exactly one; bisecting it further pins the root down to as many digits as wanted.`, table: { head: ["from", "to", "count"], rows: table } });
  const r0 = roots.find((r) => !r.isRational);
  if (r0) {
    const rows = [];
    const lo = [...r0.lo];
    const hi = [...r0.hi];
    let a = lo;
    let b = hi;
    const sa = signAt(p, a);
    for (let k = 0; k < 8; k += 1) {
      const m = reduce(a[0] * b[1] + b[0] * a[1], 2n * a[1] * b[1]);
      const s = signAt(p, m);
      rows.push([fmt(a), fmt(b), fmt(m), s > 0 ? "+" : s < 0 ? "−" : "0"]);
      if (s === sa) a = m;
      else b = m;
    }
    steps.push({ title: "Bisection", text: "Halve the interval and keep the half where p changes sign. Each step gains one bit; the library does this on demand, so the number stays exact and you choose the precision.", table: { head: ["a", "b", "midpoint", "sign of p"], rows } });
  }
  return { title: "Sturm sequences", summary: `${roots.length} real roots`, steps, result: roots };
}

/** Newton's method for a polynomial root: floats versus exact rationals, then a Krawczyk certificate. */
export function newtonKrawczyk(coeffs, { x0 = 1, iterations = 5 } = {}) {
  const p = coeffs.map((c) => toPair(c));
  const dp = p.slice(1).map((c, i) => mulQ(c, [BigInt(i + 1), 1n]));
  const ev = (q, x) => q.reduceRight((s, c) => addQ(mulQ(s, x), c), Z);
  const steps = [{ title: "Newton's method", text: `Root of p(x) = ${formatP(p)}: x ← x − p(x)/p′(x). Floats converge fast but never prove anything; exact rationals prove nothing either and their denominators roughly double every step.` }];
  let xf = Number(x0);
  let xq = toPair(x0);
  const rows = [];
  for (let k = 0; k < iterations; k += 1) {
    const pf = coeffs.reduceRight((s, c) => s * xf + Number(c), 0);
    const dpf = coeffs.slice(1).reduceRight((s, c, i) => s * xf + Number(c) * (i + 1), 0);
    xf -= pf / dpf;
    xq = subQ(xq, divQ(ev(p, xq), ev(dp, xq)));
    rows.push([`${k + 1}`, `${xf}`, fmt(xq).length > 40 ? `${fmt(xq).slice(0, 37)}…` : fmt(xq), `${xq[1].toString().length}`]);
  }
  steps.push({ title: "Iterates", text: "Float and exact Newton side by side.", table: { head: ["step", "float", "exact", "digits in denominator"], rows } });
  // Krawczyk on an interval around the float answer
  const r = 1e-6;
  let X = [toPair(xf - r), toPair(xf + r)];
  const krow = [];
  for (let k = 0; k < 4; k += 1) {
    const m = reduce(X[0][0] * X[1][1] + X[1][0] * X[0][1], 2n * X[0][1] * X[1][1]);
    const y = toPair(1 / pairToNumber(...ev(dp, m)));
    // interval p′(X) by Horner on endpoints (monotone pieces are not assumed)
    let lo = Z;
    let hi = Z;
    for (let i = dp.length - 1; i >= 0; i -= 1) {
      const prods = [mulQ(lo, X[0]), mulQ(lo, X[1]), mulQ(hi, X[0]), mulQ(hi, X[1])];
      const lt = (u, v) => u[0] * v[1] < v[0] * u[1];
      lo = addQ(prods.reduce((u, v) => (lt(u, v) ? u : v)), dp[i]);
      hi = addQ(prods.reduce((u, v) => (lt(u, v) ? v : u)), dp[i]);
    }
    // K = m − y·p(m) + (1 − y·p′(X))·(X − m)
    const c = [subQ(ONE, mulQ(y, y[0] > 0n ? hi : lo)), subQ(ONE, mulQ(y, y[0] > 0n ? lo : hi))];
    const d = [subQ(X[0], m), subQ(X[1], m)];
    const lt = (u, v) => u[0] * v[1] < v[0] * u[1];
    const prods = [mulQ(c[0], d[0]), mulQ(c[0], d[1]), mulQ(c[1], d[0]), mulQ(c[1], d[1])];
    const base = subQ(m, mulQ(y, ev(p, m)));
    const K = [addQ(base, prods.reduce((u, v) => (lt(u, v) ? u : v))), addQ(base, prods.reduce((u, v) => (lt(u, v) ? v : u)))];
    const inside = lt(X[0], K[0]) && lt(K[1], X[1]);
    krow.push([`${k + 1}`, `${pairToNumber(...X[0])} … ${pairToNumber(...X[1])}`, `${pairToNumber(...subQ(K[1], K[0])).toExponential(2)}`, inside ? "K ⊂ X: unique root" : "not yet"]);
    if (inside) X = K;
    else break;
  }
  steps.push({
    title: "Krawczyk certificate",
    text: "Take an interval X around the float answer, its midpoint m and y ≈ 1/p′(m). K(X) = m − y·p(m) + (1 − y·p′(X))(X − m) is computed with exact interval arithmetic. If K(X) lies inside X, the mean value theorem proves X holds exactly one root; repeating shrinks the interval quadratically.",
    table: { head: ["round", "X", "width of K", "verdict"], rows: krow },
  });
  return { title: "Newton and Krawczyk", summary: krow.length && krow[0][3].startsWith("K") ? "certified" : "not certified", steps, result: { float: xf } };
}

/** Matrix powers by repeated squaring. */
export function powerBySquaring(input, k) {
  const A = Matrix.from(input);
  let e = BigInt(k);
  const bitsStr = e.toString(2);
  const steps = [{ title: "Binary exponent", text: `${k} in binary is ${bitsStr}. Square A repeatedly (A, A², A⁴, …) and multiply together the powers for the 1 bits: about ${2 * bitsStr.length} products instead of ${k}.` }];
  let base = A;
  let result = Matrix.identity(A.rows);
  let power = 1n;
  const parts = [];
  while (e > 0n) {
    if (e & 1n) {
      result = result.multiply(base);
      parts.push(`A^${power}`);
      steps.push({ title: `Use A^${power}`, text: `The bit for ${power} is 1: multiply it in. Product so far: ${parts.join("·")}.`, matrix: result.toStrings() });
    }
    e >>= 1n;
    if (e > 0n) {
      base = base.multiply(base);
      power *= 2n;
      steps.push({ title: `Square to A^${power}`, text: "", matrix: base.toStrings() });
    }
  }
  return { title: "Powers by squaring", summary: `A^${k}`, steps, result: result.toStrings() };
}

/** Jordan structure from kernel dimensions of (A − λI)^k. */
export function jordanBlocks(input, lambda) {
  const A = Matrix.from(input);
  const lam = toPair(lambda);
  const n = A.rows;
  const N = A.subtract(Matrix.identity(n).scale(fmt(lam)));
  const dims = [0];
  let Nk = N;
  const steps = [{ title: `N = A − ${fmt(lam)}I`, text: "Jordan blocks for λ are read from how fast the kernels of N, N², N³, … grow. The number of blocks of size at least k is dim ker Nᵏ − dim ker Nᵏ⁻¹.", matrix: N.toStrings() }];
  for (let k = 1; k <= n; k += 1) {
    const r = gaussJordan(Nk).result.rank;
    dims.push(n - r);
    steps.push({ title: `dim ker N${k === 1 ? "" : "^" + k} = ${n - r}`, text: `rank ${r}, so the kernel has dimension ${n} − ${r} = ${n - r}.`, matrix: Nk.toStrings() });
    if (k > 1 && dims[k] === dims[k - 1]) break;
    Nk = Nk.multiply(N);
  }
  const atLeast = dims.slice(1).map((d, i) => d - dims[i]);
  const blocks = [];
  atLeast.forEach((c, i) => {
    const exactly = c - (atLeast[i + 1] ?? 0);
    for (let t = 0; t < exactly; t += 1) blocks.push(i + 1);
  });
  blocks.sort((a, b) => b - a);
  steps.push({ title: "Blocks", text: `Blocks of size ≥ k: ${atLeast.join(", ")}. So the blocks for λ = ${fmt(lam)} have sizes ${blocks.join(", ") || "none (λ is not an eigenvalue)"}.` });
  return { title: "Jordan blocks", summary: blocks.join(", "), steps, result: blocks };
}

/** Partial sums of the exponential series, and why scaling helps. */
export function exponentialSeries(input, { terms = 12 } = {}) {
  const A = Matrix.from(input).toNumbers();
  const n = A.length;
  const norm = Math.max(...A.map((r) => r.reduce((s, v) => s + Math.abs(v), 0)));
  let s = 0;
  while (norm / 2 ** s > 0.5) s += 1;
  const mul = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((acc, v, k) => acc + v * Y[k][j], 0)));
  const B = A.map((r) => r.map((v) => v / 2 ** s));
  let term = Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0)));
  let sum = term.map((r) => [...r]);
  const rows = [];
  for (let k = 1; k <= terms; k += 1) {
    term = mul(term, B).map((r) => r.map((v) => v / k));
    sum = sum.map((r, i) => r.map((v, j) => v + term[i][j]));
    rows.push([`${k}`, sum[0].map((v) => v.toPrecision(10)).join(", "), Math.max(...term.flat().map(Math.abs)).toExponential(1)]);
  }
  const steps = [
    { title: "The series", text: "exp(A) = I + A + A²/2! + A³/3! + …. It always converges, but slowly when A is large, and large terms cancel badly in floating point." },
    { title: `Scale by 2^${s}`, text: `‖A‖ = ${norm}. Use B = A/2^${s}, whose norm is at most 1/2, so the terms shrink at least as fast as 1/(2ᵏ·k!). Then exp(A) = exp(B)^(2^${s}): square ${s} times.` },
    { title: "Partial sums of exp(B), first row", text: "The last column is the size of the term just added; it bounds the error that remains. The library adds a rigorous bound for the whole tail to every entry, as an interval.", table: { head: ["k", "first row of the sum", "term size"], rows } },
  ];
  for (let k = 0; k < s; k += 1) sum = mul(sum, sum);
  steps.push({ title: "Square back", text: `After ${s} squarings: first row ${sum[0].map((v) => v.toPrecision(12)).join(", ")}.` });
  return { title: "The matrix exponential", summary: "series", steps, result: sum };
}


/**
 * The four fundamental subspaces of A, read off the rref of A and of Aᵀ.
 */
export function fourSubspaces(input) {
  const A = Matrix.from(input);
  const m = A.rows;
  const n = A.cols;
  const g = gaussJordan(A);
  const { pivots, pairs: Rm } = g.result;
  const r = pivots.length;
  const free = [...Array(n).keys()].filter((j) => !pivots.includes(j));
  const steps = [{ title: "A", text: `A is ${m}×${n}: it sends vectors with ${n} entries to vectors with ${m} entries. Four subspaces describe it completely; one row reduction finds three of them.`, matrix: A.toStrings() }];
  steps.push({ title: "Row reduce", text: `Row operations do not change which combinations of columns vanish, so the reduced form R has the same null space and row space as A. Pivot columns: ${pivots.map((c) => c + 1).join(", ") || "none"}, so the rank is ${r}.`, matrix: g.result.matrix, cols: pivots });
  steps.push({ title: "Column space", text: r ? `The pivot columns of A itself (not of R) are a basis for the column space, the set of all A·x: columns ${pivots.map((c) => c + 1).join(", ")}. Its dimension is the rank, ${r}.` : "A is zero, so its column space is {0}.", matrix: A.toStrings(), cols: pivots });
  steps.push({ title: "Row space", text: r ? `The ${r} nonzero rows of R are a basis for the row space. They are independent because each has a 1 where the others have 0.` : "The row space is {0}.", matrix: snap(Rm), rows: [...Array(r).keys()] });
  const nullBasis = free.map((f) => {
    const v = Array(n).fill(Z);
    v[f] = ONE;
    pivots.forEach((pc, i) => { v[pc] = neg(Rm[i][f]); });
    return v;
  });
  steps.push({
    title: "Null space",
    text: free.length
      ? `Free columns ${free.map((c) => c + 1).join(", ")} give one solution of A·x = 0 each: set that free variable to 1, the other free ones to 0, and read the pivot variables off R. These ${free.length} vectors are a basis; rank + nullity = ${r} + ${free.length} = ${n}.`
      : `Every column has a pivot, so A·x = 0 only for x = 0: the null space is {0}, and rank ${r} = ${n} columns.`,
    matrix: free.length ? nullBasis[0].map((_, i) => nullBasis.map((v) => fmt(v[i]))) : undefined,
  });
  const gt = gaussJordan(A.transpose());
  const pt = gt.result.pivots;
  const freeT = [...Array(m).keys()].filter((j) => !pt.includes(j));
  const leftBasis = freeT.map((f) => {
    const v = Array(m).fill(Z);
    v[f] = ONE;
    pt.forEach((pc, i) => { v[pc] = neg(gt.result.pairs[i][f]); });
    return v;
  });
  steps.push({
    title: "Left null space",
    text: freeT.length
      ? `The left null space holds the y with yᵀA = 0, the null space of Aᵀ. Row reducing Aᵀ the same way gives ${freeT.length} basis vectors (as columns). Each one says which combination of the equations cancels, which is exactly the condition for A·x = b to be solvable: yᵀb = 0.`
      : `Aᵀ has a pivot in every column, so the left null space is {0}: A·x = b is solvable for every b.`,
    matrix: freeT.length ? leftBasis[0].map((_, i) => leftBasis.map((v) => fmt(v[i]))) : undefined,
  });
  steps.push({ title: "The picture", text: `Row space (dimension ${r}) ⟂ null space (dimension ${n - r}) inside Qⁿ, n = ${n}; column space (dimension ${r}) ⟂ left null space (dimension ${m - r}) inside Qᵐ, m = ${m}. A maps the row space one-to-one onto the column space.` });
  return { title: "The four fundamental subspaces", summary: `rank ${r}`, steps, result: { rank: r, nullity: n - r, leftNullity: m - r } };
}

/**
 * Dixon's p-adic lifting for a square system: invert A once modulo a small
 * prime, peel off the solution's base-p digits, then recover fractions.
 */
export function dixonLifting(aInput, bInput, { p = 97 } = {}) {
  const A = Matrix.from(aInput);
  let B = Matrix.from(bInput);
  if (B.rows !== A.rows && B.rows === 1) B = B.transpose();
  const n = A.rows;
  if (!A.isSquare) throw new RangeError("dixonLifting needs a square matrix");
  // integer system: scale each row by the lcm of its denominators
  const { rows: ints } = integerRows(A.augment(B).pairs());
  const Ai = ints.map((r) => r.slice(0, n));
  const bi = ints.map((r) => r[n]);
  const P = BigInt(p);
  const mod = (x) => ((x % P) + P) % P;
  const steps = [{ title: "The system", text: `Clear denominators row by row so every entry is an integer. The plan: find x modulo ${p}, then modulo ${p}², ${p}³, …, one base-${p} digit at a time, using only arithmetic below ${p} and integer updates.`, matrix: ints.map((r) => r.map(String)), bar: n }];
  // inverse mod p by Gauss–Jordan
  const W = Ai.map((r, i) => [...r.map(mod), ...Ai.map((_, j) => (i === j ? 1n : 0n))]);
  const invP = (a) => { let [r0, r1, s0, s1] = [P, mod(a), 0n, 1n]; while (r1) { const q = r0 / r1; [r0, r1, s0, s1] = [r1, r0 - q * r1, s1, s0 - q * s1]; } return mod(s0); };
  for (let c = 0; c < n; c += 1) {
    const f = W.findIndex((row, i) => i >= c && mod(row[c]) !== 0n);
    if (f < 0) {
      steps.push({ title: `A is singular mod ${p}`, text: `No pivot in column ${c + 1} modulo ${p}. Either A is singular, or ${p} is an unlucky prime dividing det A; the library would just pick another prime.` });
      return { title: "Dixon lifting", summary: "unlucky prime", steps, result: null };
    }
    [W[c], W[f]] = [W[f], W[c]];
    const iv = invP(W[c][c]);
    W[c] = W[c].map((x) => mod(x * iv));
    for (let i = 0; i < n; i += 1) if (i !== c && W[i][c] !== 0n) { const k = W[i][c]; W[i] = W[i].map((x, j) => mod(x - k * W[c][j])); }
  }
  const C = W.map((r) => r.slice(n));
  steps.push({ title: `Invert A modulo ${p}`, text: `Gauss–Jordan with arithmetic mod ${p} gives C with A·C ≡ I (mod ${p}). This is the only elimination in the whole method, and every number in it is below ${p}.`, matrix: C.map((r) => r.map(String)) });
  // Hadamard-style bound to know how many digits we need
  const norm = (v) => v.reduce((a, x) => a + x * x, 0n);
  let hb = 1n;
  for (const r of Ai) hb *= sqrtCeil(norm(r));
  let nb = hb;
  for (let j = 0; j < n; j += 1) {
    let prod = sqrtCeil(norm(bi));
    for (let k = 0; k < n; k += 1) if (k !== j) prod *= sqrtCeil(norm(Ai.map((r) => r[k])));
    if (prod > nb) nb = prod;
  }
  const need = 2n * nb * hb + 1n;
  let r = [...bi];
  let x = Array(n).fill(0n);
  let M = 1n;
  const table = [];
  let k = 0;
  while (M <= need) {
    const digit = C.map((row) => mod(row.reduce((a, c, j) => a + c * r[j], 0n)));
    const Ad = Ai.map((row) => row.reduce((a, c, j) => a + c * digit[j], 0n));
    const nr = r.map((v, i) => (v - Ad[i]) / P);
    x = x.map((v, i) => v + digit[i] * M);
    M *= P;
    k += 1;
    table.push([String(k), `(${digit.join(", ")})`, `(${x.join(", ")})`]);
    if (k <= 3) {
      steps.push({ title: `Digit ${k}`, text: `x_${k} = C·r mod ${p} = (${digit.join(", ")}). Then r ← (r − A·x_${k}) / ${p}, which divides exactly because A·x_${k} ≡ r (mod ${p}). Now r = (${nr.join(", ")}). The solution so far, x mod ${p}${k > 1 ? `^${k}` : ""}: (${x.join(", ")}).`, table: { head: ["k", "digit", `x mod ${p}^k`], rows: [...table] } });
    }
    r = nr;
  }
  steps.push({ title: `${k} digits`, text: `Cramer's rule and Hadamard's bound say numerators and the denominator are at most ${nb} and ${hb}, so ${k} base-${p} digits (modulus above 2·${nb}·${hb}) are enough.${table.some((t) => t[1].includes(String(p - 1))) ? ` A digit of ${p - 1} that keeps repeating is how a negative number looks in base ${p}: −1 = …(${p - 1})(${p - 1})(${p - 1}).` : ""}`, table: { head: ["k", "digit", `x mod ${p}^k`], rows: table } });
  const out = x.map((v) => {
    let [r0, r1] = [M, mod(v) === 0n ? 0n : ((v % M) + M) % M];
    let [t0, t1] = [0n, 1n];
    while (r1 * r1 * 2n > M) { const q = r0 / r1; [r0, r1] = [r1, r0 - q * r1]; [t0, t1] = [t1, t0 - q * t1]; }
    return reduce(t1 < 0n ? -r1 : r1, t1 < 0n ? -t1 : t1);
  });
  steps.push({ title: "Rational reconstruction", text: `Each entry is known modulo ${M}. The extended Euclidean algorithm, stopped halfway, finds the unique fraction with small numerator and denominator matching it: (${out.map(fmt).join(", ")}). A final exact check A·x = b confirms it.`, matrix: out.map((v) => [fmt(v)]) });
  const ok = A.multiply(Matrix.fromPairs(out.map((v) => [v]))).equals(B);
  if (!ok) steps.push({ title: "Check failed", text: "The reconstruction does not satisfy the system; the library falls back to elimination in that case." });
  return { title: "Dixon p-adic lifting", summary: out.map(fmt).join(", "), steps, result: ok ? out.map(fmt) : null };
}

function sqrtCeil(n) {
  if (n < 2n) return n;
  let x = BigInt(Math.ceil(Math.sqrt(Number(n))));
  while (x * x > n) x -= 1n;
  while (x * x < n) x += 1n;
  return x;
}

export const WALKTHROUGHS = {
  gaussJordan, solveSystem, inverse, determinantByElimination, cofactorExpansion, bareiss, modularDeterminant,
  guessAndCheck, gramSchmidt, lu, leastSquares, eigen, sturmWalk, newtonKrawczyk, powerBySquaring, jordanBlocks, exponentialSeries,
  fourSubspaces, dixonLifting,
};

export { gcd, lcm };
