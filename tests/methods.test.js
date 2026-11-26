// The teaching methods, checked against exact arithmetic and the theory.
import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { methods as M, eigenvalues } from "../index.js";

const rng = (seed) => () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
const randInt = (r, lo, hi) => lo + Math.floor(r() * (hi - lo + 1));

describe("doubles as exact rationals", () => {
  test("exactDouble gives the value the double is", () => {
    const [n, d] = M.exactDouble(0.1);
    expect(d).toBe(2n ** 55n);
    expect(n).toBe(3602879701896397n);
    expect(M.exactDouble(-3)).toEqual([-3n, 1n]);
    const r = rng(1);
    for (let i = 0; i < 500; i++) {
      const x = (r() - 0.5) * 10 ** randInt(r, -30, 30);
      const [p, q] = M.exactDouble(x);
      expect(new Rational(p, q).toNumber()).toBe(x);
    }
  });
});

describe("stationary iterations", () => {
  test("Jacobi's iteration matrix is −D⁻¹(L + U), and the error obeys e(k+1) = T e(k) exactly", () => {
    const A = [[5, 1, 2], [1, 6, -2], [2, 1, 7]];
    const run = M.jacobi(A, [1, 2, 3], { iterations: 6 });
    expect(run.iterationMatrix).toEqual([["0", "-1/5", "-2/5"], ["-1/6", "0", "1/3"], ["-2/7", "-1/7", "0"]]);
    expect(run.strictlyDiagonallyDominant).toBe(true);
    expect(run.converges).toBe(true);
  });
  test("for a 2×2 system Gauss–Seidel's spectral radius is the square of Jacobi's", () => {
    const r = rng(2);
    for (let t = 0; t < 40; t++) {
      const A = [[randInt(r, 3, 9), randInt(r, -3, 3)], [randInt(r, -3, 3), randInt(r, 3, 9)]];
      const j = M.jacobi(A, [1, 1], { iterations: 2 }).spectralRadius;
      const g = M.gaussSeidel(A, [1, 1], { iterations: 2 }).spectralRadius;
      expect(g).toBeCloseTo(j * j, 12);
    }
  });
  test("the observed error ratio tends to the spectral radius, and rounding stays at the level of the last bit", () => {
    const run = M.gaussSeidel([[3, -1, 0], [-1, 3, -1], [0, -1, 3]], [1, 2, 3], { iterations: 15, exactIterations: 15 });
    expect(run.ratio).toBeCloseTo(run.spectralRadius, 6);
    for (const row of run.rows) expect(row.rounding).toBeLessThan(1e-15);
    expect(run.rows.at(-1).methodError).toBeLessThan(1e-9);
  });
  test("ω = 1 is Gauss–Seidel; a spectral radius above 1 diverges", () => {
    const A = [[4, -1, 0], [-1, 4, -1], [0, -1, 4]];
    const g = M.gaussSeidel(A, [15, 10, 10], { iterations: 5 });
    const s = M.sor(A, [15, 10, 10], { omega: 1, iterations: 5 });
    expect(s.rows.map((r) => r.x)).toEqual(g.rows.map((r) => r.x));
    expect(M.sor(A, [15, 10, 10], { omega: "11/10" }).spectralRadius).toBeCloseTo(0.1, 12);
    const bad = M.jacobi([[1, 2], [2, 1]], [3, 3], { iterations: 6 });
    expect(bad.converges).toBe(false);
    expect(bad.rows.map((r) => r.error)).toEqual([1, 2, 4, 8, 16, 32, 64]);
  });
  test("a zero on the diagonal and a singular matrix are refused", () => {
    expect(() => M.jacobi([[0, 1], [1, 0]], [1, 1])).toThrow(/diagonal/);
    expect(() => M.gaussSeidel([[1, 2], [2, 4]], [1, 1])).toThrow(/singular/);
  });
});

describe("conjugate gradients", () => {
  test("exact in at most n steps on random symmetric positive definite matrices", () => {
    const r = rng(3);
    for (let t = 0; t < 25; t++) {
      const n = randInt(r, 2, 5);
      const B = Array.from({ length: n }, () => Array.from({ length: n }, () => randInt(r, -3, 3)));
      const A = B.map((_, i) => B.map((_, j) => B.reduce((s, row) => s + row[i] * row[j], 0) + (i === j ? 1 : 0)));
      const b = Array.from({ length: n }, () => randInt(r, -5, 5));
      if (b.every((v) => v === 0)) b[0] = 1;
      const cg = M.conjugateGradient(A, b);
      expect(cg.exactSteps).toBeGreaterThan(0);
      expect(cg.exactSteps).toBeLessThanOrEqual(n);
      expect(cg.rows[cg.exactSteps].error).toBeLessThan(1e-9);
    }
  });
  test("refuses a matrix that is not symmetric, or not positive definite", () => {
    expect(() => M.conjugateGradient([[1, 2], [0, 1]], [1, 1])).toThrow(/symmetric/);
    expect(() => M.conjugateGradient([[1, 2], [2, 1]], [1, 0])).toThrow(/positive definite/);
  });
});

describe("eigenvalue methods", () => {
  test("the power method converges to the dominant eigenvalue at the predicted rate", () => {
    const p = M.powerMethod([[2, 1], [1, 3]], { iterations: 12 });
    expect(p.target).toBeCloseTo((5 + Math.sqrt(5)) / 2, 14);
    expect(p.convergenceFactor).toBeCloseTo((5 - Math.sqrt(5)) / (5 + Math.sqrt(5)), 14);
    // symmetric: the Rayleigh quotient's error falls like the square of the factor
    const e = p.rows.map((r) => r.methodError).filter((v) => v !== null && v > 1e-14);
    expect(e.at(-1) / e.at(-2)).toBeCloseTo(p.convergenceFactor ** 2, 2);
  });
  test("inverse iteration finds the eigenvalue nearest the shift", () => {
    const A = [[4, 1, 0], [1, 2, 1], [0, 1, -3]];
    const all = M.powerMethod(A, { iterations: 2 }).eigenvalues.map((e) => e.re);
    for (const shift of ["-2", "1", "7/2"]) {
      const want = all.reduce((b, v) => (Math.abs(v - new Rational(shift).toNumber()) < Math.abs(b - new Rational(shift).toNumber()) ? v : b));
      expect(M.inverseIteration(A, { shift, iterations: 30 }).rows.at(-1).estimate).toBeCloseTo(want, 10);
    }
    expect(() => M.inverseIteration([[2, 0], [0, 3]], { shift: 2 })).toThrow(/eigenvalue/);
  });
  test("no dominant eigenvalue is reported, not hidden", () => {
    const p = M.powerMethod([[0, 1], [-1, 0]], { iterations: 3 });
    expect(p.target).toBe(null);
    expect(p.note).toMatch(/strictly largest/);
  });
  test("the QR algorithm's diagonal approaches the exact eigenvalues", () => {
    const A = [[2, 1, 0], [1, 3, 1], [0, 1, 4]];
    const q = M.qrAlgorithm(A, { iterations: 60 });
    const exact = eigenvalues(A).exact().real.map((r) => r.value.toNumber()).sort((a, b) => b - a);
    q.rows.at(-1).diagonal.forEach((v, i) => expect(v).toBeCloseTo(exact[i], 9));
    expect(q.rows.at(-1).subdiagonal).toBeLessThan(1e-9);
  });
});

describe("pivoting", () => {
  test("a tiny pivot without row exchanges loses the answer; partial pivoting keeps it", () => {
    const r = M.pivoting([[1e-17, 1], [1, 1]], [1, 2]);
    expect(r.none.error).toBeGreaterThan(0.5);
    expect(r.none.growth).toBeGreaterThan(1e16);
    expect(r.partial.error).toBeLessThan(1e-15);
    expect(r.partial.growth).toBe(1);
  });
});
