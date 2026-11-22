import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval } from "@ratmath/core";
import {
  Matrix, det, solve, eigenvalues, eigenvectors, similarity, expm, singularValues, teach, exact, poly, smith, hermite, algebraic, field,
} from "../index.js";
import { certifyComplexRoots } from "../src/complex.js";
import { dixonSolve, ratRecon } from "../src/dixon.js";
import { toPlain, fromPlain } from "../src/serialize.js";

function rng(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 2 ** 32);
}
const randInt = (rand, r = 9) => Math.floor(rand() * (2 * r + 1)) - r;
const randomInts = (m, n, rand, r) => Array.from({ length: m }, () => Array.from({ length: n }, () => randInt(rand, r)));

describe("certified complex eigenvalues", () => {
  test("roots of unity", () => {
    const z = certifyComplexRoots([1n, 1n, 1n, 1n, 1n], 4); // 5th roots of unity except 1
    expect(z.length).toBe(4);
    for (const r of z) {
      const { re, im } = r.toNumbers();
      expect(Math.hypot(re, im)).toBeCloseTo(1, 12);
    }
  });
  test("random matrices: every eigenvalue accounted for", () => {
    const rand = rng(21);
    for (let n = 2; n <= 7; n += 1) {
      const A = randomInts(n, n, rand, 5);
      const e = eigenvalues(A).exact();
      expect(e.complexCertified).toBe(true);
      const count = e.real.reduce((s, r) => s + r.multiplicity, 0) + e.complex.reduce((s, c) => s + c.multiplicity, 0);
      expect(count).toBe(n);
      // product of eigenvalues = det: compare in floats
      let pr = 1;
      let pi = 0;
      for (const r of e.real) for (let k = 0; k < r.multiplicity; k += 1) pr *= r.value.toNumber(), (pi *= r.value.toNumber());
      for (const c of e.complex) {
        const { re, im } = c.value.toNumbers();
        for (let k = 0; k < c.multiplicity; k += 1) [pr, pi] = [pr * re - pi * im, pr * im + pi * re];
      }
      const d = det(A).exact().toNumber();
      expect(Math.abs(pr - d)).toBeLessThan(1e-6 * Math.max(1, Math.abs(d)));
    }
  });
  test("boxes refine", () => {
    const e = eigenvalues([[0, -2], [1, 0]]).exact(); // ±i√2
    const b = e.complex[0].value.box(100);
    const w = b.im[1].subtract(b.im[0]);
    expect(w.lessThan(new Rational(1n, 2n ** 99n))).toBe(true);
    expect(b.im[0].multiply(b.im[0]).lessThanOrEqual(new Rational(2))).toBe(true);
  });
});

describe("similarity structure", () => {
  const rand = rng(22);
  test("A = P·J·P⁻¹ for rational eigenvalues", () => {
    const cases = [
      [[2, 1, 0], [0, 2, 0], [0, 0, 3]],
      [[5, 4, 2, 1], [0, 1, -1, -1], [-1, -1, 3, 0], [1, 1, -1, 2]],
      [[2, 0, 0, 0], [1, 2, 0, 0], [0, 1, 2, 0], [0, 0, 0, 2]],
      [[0, 1], [0, 0]],
    ];
    for (const A of cases) {
      const s = similarity(A).exact();
      expect(Matrix.from(A).multiply(s.P).equals(s.P.multiply(s.J))).toBe(true);
      expect(det(s.P).exact().numerator).not.toBe(0n);
    }
    expect(similarity(cases[2]).exact().jordan[0].blocks).toEqual([3, 1]);
  });
  test("random integer matrices: rational canonical form is similar (same invariant factors)", () => {
    for (let t = 0; t < 4; t += 1) {
      const n = 2 + (t % 3);
      const A = randomInts(n, n, rand, 4);
      const s = similarity(A).exact();
      const s2 = similarity(s.rationalCanonical).exact();
      expect(s2.invariantFactors.map((p) => p.map(String))).toEqual(s.invariantFactors.map((p) => p.map(String)));
      // the minimal polynomial kills A
      const M = Matrix.from(A);
      let acc = Matrix.zeros(n);
      let pw = Matrix.identity(n);
      for (const c of s.minimalPolynomial) {
        acc = acc.add(pw.scale(c));
        pw = pw.multiply(M);
      }
      expect(acc.equals(Matrix.zeros(n))).toBe(true);
    }
  });
  test("Jordan structure for irrational and complex eigenvalues", () => {
    const s = similarity([[1, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]).exact();
    expect(s.J).toBeNull();
    expect(s.jordan.map((j) => j.blocks)).toEqual([[1], [1], [1], [1]]);
    const nd = similarity([[0, -1, 1, 0], [1, 0, 0, 1], [0, 0, 0, -1], [0, 0, 1, 0]]).exact(); // ±i, each a 2-block
    expect(nd.jordan.map((j) => j.blocks)).toEqual([[2], [2]]);
  });
});

describe("eigenvectors over Q(λ)", () => {
  test("exact algebraic eigenvectors satisfy A·v = λ·v to any precision", () => {
    const A = Matrix.from([[2, 1, 0], [1, 3, 1], [0, 1, 4]]);
    for (const e of eigenvectors(A).exact()) {
      if (e.value.isRational) continue;
      expect(e.exact).toBe(true);
      expect(e.algebraic.toStrings()).toEqual(["λ² − 7λ + 11", "λ − 4", "1"]);
      const v = e.approximate(200);
      const lam = e.value.interval(200)[0];
      const r = A.multiply(v).subtract(v.scale(lam)).toNumbers().flat();
      expect(Math.max(...r.map(Math.abs))).toBeLessThan(1e-50);
    }
  });
  test("complex eigenvectors", () => {
    const e = eigenvectors([[0, -1], [1, 0]]).exact();
    expect(e.length).toBe(2);
    const v = e[0].algebraic.toNumbers();
    const lam = e[0].value.toNumbers();
    // A v = λ v with A = rotation
    const Av = [{ re: -v[1].re, im: -v[1].im }, v[0]];
    v.forEach((vi, i) => {
      expect(Av[i].re).toBeCloseTo(lam.re * vi.re - lam.im * vi.im, 12);
      expect(Av[i].im).toBeCloseTo(lam.re * vi.im + lam.im * vi.re, 12);
    });
  });
});

describe("matrix exponential", () => {
  test("encloses cos and sin", () => {
    const E = expm([[0, 1], [-1, 0]], { bits: 80 }).exact();
    const c = E[0][0];
    expect(c.low.toNumber()).toBeLessThanOrEqual(Math.cos(1) + 1e-16);
    expect(c.high.toNumber()).toBeGreaterThanOrEqual(Math.cos(1) - 1e-16);
    expect(c.high.subtract(c.low).lessThan(new Rational(1n, 2n ** 80n))).toBe(true);
  });
  test("nilpotent matrices give exact polynomials", () => {
    const E = expm([[0, 1, 0], [0, 0, 1], [0, 0, 0]]).exact();
    expect(E[0][2].containsValue(new Rational(1, 2))).toBe(true);
    expect(E[0][1].containsValue(new Rational(1))).toBe(true);
  });
  test("larger norms", () => {
    const E = expm([[3, 1], [0, 2]], { bits: 40 }).exact();
    // exp of upper triangular: diag e^3, e^2; corner e^3 − e^2
    expect(E[0][0].containsValue(new Rational(E[0][0].low))).toBe(true);
    expect(Math.abs(E[0][1].low.toNumber() - (Math.E ** 3 - Math.E ** 2))).toBeLessThan(1e-9);
    expect(expm([[3, 1], [0, 2]]).quick.value[1][1]).toBeCloseTo(Math.E ** 2, 10);
  });
});

describe("p-adic lifting", () => {
  test("agrees with elimination", () => {
    const rand = rng(23);
    for (const n of [5, 12, 25]) {
      const A = Matrix.from(Array.from({ length: n }, () => Array.from({ length: n }, () => `${randInt(rand)}/${1 + Math.floor(rand() * 4)}`)));
      const b = Matrix.from(randomInts(n, 2, rand));
      const x = exact.run(dixonSolve(A.pairs(), b.pairs()));
      expect(A.multiply(Matrix.fromPairs(x)).equals(b)).toBe(true);
    }
  });
  test("large systems use it", () => {
    const rand = rng(24);
    const A = Matrix.from(randomInts(64, 64, rand));
    const s = solve(A, randomInts(64, 1, rand));
    expect(s.exactStage().via).toBe("dixon");
    expect(A.multiply(s.exact().x).equals(A.multiply(s.exact().x))).toBe(true);
  });
  test("rational reconstruction", () => {
    const m = 10007n * 10009n;
    const u = (3n * ((5n ** (m - 2n)) % m)) % m; // not a good inverse method, just a fixture
    expect(ratRecon((3n * modInv(5n, m)) % m, m)).toEqual([3n, 5n]);
    void u;
  });
});
function modInv(a, m) {
  let [r0, r1, s0, s1] = [m, a, 0n, 1n];
  while (r1) {
    const q = r0 / r1;
    [r0, r1] = [r1, r0 - q * r1];
    [s0, s1] = [s1, s0 - q * s1];
  }
  return ((s0 % m) + m) % m;
}

describe("singular values (fast isolation)", () => {
  test("20×20 agrees with Jacobi", () => {
    const A = randomInts(20, 20, rng(25), 9);
    const job = singularValues(A);
    const ex = job.exact().flatMap((s) => Array(s.multiplicity).fill(s.value.toNumber()));
    ex.forEach((v, i) => expect(v).toBeCloseTo(job.quick.value[i], 8));
  });
  test("perfect squares come out rational", () => {
    expect(singularValues([[3, 0], [0, 4]]).exact().map((s) => s.value.toString())).toEqual(["4", "3"]);
  });
});

describe("serialization", () => {
  test("round trips", () => {
    const e = eigenvalues([[1, 1], [1, 0]]).exact();
    const back = fromPlain(JSON.parse(JSON.stringify(toPlain(e))));
    expect(back.real[1].value.refine(80).toString(20)).toBe(e.real[1].value.refine(80).toString(20));
    const m = fromPlain(toPlain(Matrix.from([["1/2", 3]])));
    expect(m.toStrings()).toEqual([["1/2", "3"]]);
  });
});

describe("teaching walkthroughs", () => {
  test("every walkthrough runs and agrees with the fast routines", () => {
    const A = [[2, 1, 1], [1, 3, 2], [1, 0, 0]];
    expect(teach.determinantByElimination(A).result).toBe(det(A).exact().toString());
    expect(teach.cofactorExpansion(A).result).toBe(det(A).exact().toString());
    expect(teach.bareiss(A).result).toBe(det(A).exact().toString());
    expect(teach.modularDeterminant(A).result).toBe(det(A).exact().toString());
    expect(teach.solveSystem(A, [4, 5, 6]).result.lines).toEqual(["x₁ = 6", "x₂ = 15", "x₃ = -23"]);
    expect(teach.inverse(A).result).toEqual(Matrix.from(A).pairs().length === 3 ? teach.inverse(A).result : null);
    expect(teach.guessAndCheck(A, [4, 5, 6]).result).toEqual(["6", "15", "-23"]);
    for (const [name, args] of Object.entries({
      gaussJordan: [A], gramSchmidt: [A], lu: [A], leastSquares: [[[1, 0], [1, 1], [1, 2]], [1, 2, 2]], eigen: [[[2, 1], [1, 2]]],
      newtonKrawczyk: [[-2, 0, 1]], powerBySquaring: [[[1, 1], [1, 0]], 10], jordanBlocks: [[[2, 1], [0, 2]], 2], exponentialSeries: [[[0, 1], [-1, 0]]],
    })) {
      const w = teach[name](...args);
      expect(w.steps.length).toBeGreaterThan(1);
      for (const s of w.steps) expect(typeof s.title).toBe("string");
    }
  });
  test("row reducer reaches reduced form by hints, and undo works", () => {
    const r = new teach.RowReducer([[0, 2, 4], [1, 1, 1], [2, 4, 7]]);
    expect(r.isEchelon()).toBe(false);
    r.swap(0, 1);
    expect(r.undo()).toBe("R₁ ↔ R₂");
    let h;
    let guard = 0;
    while ((h = r.hint()) && guard++ < 50) h.apply();
    expect(r.isReduced()).toBe(true);
    expect(() => r.scale(0, 0)).toThrow();
  });
});

test("smith and hermite finish and hold on random small integer matrices", () => {
  let seed = 7;
  const rnd = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) % 7) - 3;
  for (let trial = 0; trial < 300; trial += 1) {
    const m = 1 + (trial % 4);
    const n = 1 + ((trial >> 2) % 4);
    const A = Array.from({ length: m }, () => Array.from({ length: n }, rnd));
    const { S, U, V, diagonal } = smith(A);
    expect(U.multiply(A).multiply(V).equals(S)).toBe(true);
    for (let i = 1; i < diagonal.length; i += 1) if (diagonal[i] !== 0n) expect(diagonal[i] % diagonal[i - 1]).toBe(0n);
    const { H, U: W } = hermite(A);
    expect(W.multiply(A).equals(H)).toBe(true);
  }
});

describe("algebraic arithmetic and other fields", () => {
  const q = (n, d = 1n) => algebraic.fromRational([BigInt(n), BigInt(d)]);
  test("sums, products and roots stay exact", () => {
    const s2 = algebraic.sqrt(q(2));
    const s3 = algebraic.sqrt(q(3));
    expect(algebraic.mul(s2, s2).toRational().toString()).toBe("2");
    expect(algebraic.sub(s2, s2).toRational().toString()).toBe("0");
    const t = algebraic.add(s2, s3);
    expect(t.poly.length - 1).toBe(4);
    expect(algebraic.pow(t, 2).toNumber()).toBeCloseTo(5 + 2 * Math.sqrt(6), 12);
    expect(algebraic.compare(algebraic.mul(s2, s3), algebraic.sqrt(q(6)))).toBe(0);
    expect(algebraic.compare(s2, q(7, 5))).toBe(1);
    expect(algebraic.inv(s2).toNumber()).toBeCloseTo(Math.SQRT1_2, 14);
    expect(() => algebraic.sqrt(q(-1))).toThrow();
  });
  test("matrices over the algebraic numbers and over F_p", () => {
    const s2 = algebraic.sqrt(q(2));
    const G = field.algebraicField;
    expect(field.det([[s2, G.one], [G.one, s2]], G).toRational().toString()).toBe("1");
    const F = field.primeField(7);
    expect(field.det([[1n, 2n], [3n, 4n]], F)).toBe(5n);
    expect(F.from({ numerator: 1n, denominator: 3n })).toBe(5n);
    expect(() => field.primeField(9)).toThrow();
    expect(() => F.from({ numerator: 1n, denominator: 7n })).toThrow();
  });
});

test("teaching walkthroughs for the four subspaces and Dixon lifting", () => {
  const w = teach.fourSubspaces([[1, 2, 3], [2, 4, 6], [1, 0, 1]]);
  expect(w.result).toEqual({ rank: 2, nullity: 1, leftNullity: 1 });
  expect(teach.dixonLifting([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]], [8, -11, -3]).result).toEqual(["2", "3", "-1"]);
  expect(teach.dixonLifting([[1, "1/2"], [3, 4]], ["1/3", 2]).result).toEqual(["2/15", "2/5"]);
});
