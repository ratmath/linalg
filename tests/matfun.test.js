import { describe, expect, test } from "bun:test";
import { Rational } from "@ratmath/core";
import { cosm, sinm, coshm, sinhm, sqrtm, logm, powm, expm } from "../index.js";

const mid = (x) => (x.low.toNumber() + x.high.toNumber()) / 2;
const width = (E) => Math.max(...E.flat().map((x) => x.high.subtract(x.low).toNumber()));
const contains = (x, q) => x.low.lessThanOrEqual(q) && x.high.greaterThanOrEqual(q);
const R = [[0, -1], [1, 0]];

/** Interval product of two matrices of RationalIntervals. */
const times = (X, Y) => X.map((row) => Y[0].map((_, j) => row.reduce((acc, x, k) => (acc ? acc.add(x.multiply(Y[k][j])) : x.multiply(Y[k][j])), null)));

describe("sin and cos of a matrix", () => {
  test("a rotation generator gives cosh 1 and sinh 1", () => {
    const C = cosm(R).exact();
    const S = sinm(R).exact();
    expect(mid(C[0][0])).toBeCloseTo(Math.cosh(1), 14);
    expect(mid(S[1][0])).toBeCloseTo(Math.sinh(1), 14);
    expect(contains(C[0][1], new Rational(0))).toBe(true);
    expect(width(C)).toBeLessThan(2 ** -53);
  });
  test("sin² + cos² encloses the identity", () => {
    const A = [[1, 2, 0], ["1/2", -1, 3], [0, 1, "3/4"]];
    const C = cosm(A).exact();
    const S = sinm(A).exact();
    const CC = times(C, C);
    const SS = times(S, S);
    for (let i = 0; i < 3; i += 1) for (let j = 0; j < 3; j += 1) expect(contains(CC[i][j].add(SS[i][j]), new Rational(i === j ? 1 : 0))).toBe(true);
  });
  test("a nilpotent matrix: cos is the identity up to the proven width", () => {
    const C = cosm([[0, 1], [0, 0]]).exact();
    expect(contains(C[0][0], new Rational(1))).toBe(true);
    expect(contains(C[0][1], new Rational(0))).toBe(true);
  });
  test("cosh and sinh add up to exp", () => {
    const A = [[1, 2], [3, -1]];
    const E = expm(A).exact();
    const C = coshm(A).exact();
    const S = sinhm(A).exact();
    for (let i = 0; i < 2; i += 1) for (let j = 0; j < 2; j += 1) expect(mid(C[i][j]) + mid(S[i][j])).toBeCloseTo(mid(E[i][j]), 10);
  });
  test("bits narrows the enclosure and the preview agrees", () => {
    const r = sinm(R, { bits: 150 });
    expect(r.quick.value[1][0]).toBeCloseTo(Math.sinh(1), 12);
    expect(width(r.exact())).toBeLessThan(2 ** -150);
  });
  test("a non-square matrix is refused", () => {
    expect(() => cosm([[1, 2, 3]])).toThrow(RangeError);
  });
});

describe("square root", () => {
  test("a rational root is returned exactly", () => {
    const S = sqrtm([[1, 1], [0, 1]]).exact();
    expect(S[0][1].low.toString()).toBe("1/2");
    expect(S[0][1].high.toString()).toBe("1/2");
    expect(sqrtm([[4, 0], [0, 9]]).exact()[1][1].low.toString()).toBe("3");
  });
  test("the square of the enclosure contains A", () => {
    const A = [[2, 1], [1, 3]];
    const S = sqrtm(A, { bits: 80 }).exact();
    const SS = times(S, S);
    for (let i = 0; i < 2; i += 1) for (let j = 0; j < 2; j += 1) expect(contains(SS[i][j], new Rational(A[i][j]))).toBe(true);
    expect(width(S)).toBeLessThan(2 ** -80);
    expect(mid(S[0][0])).toBeCloseTo(1.3763819204711736, 14);
  });
  test("a quarter turn has the eighth turn as its principal root", () => {
    const S = sqrtm(R).exact();
    expect(mid(S[0][0])).toBeCloseTo(Math.SQRT1_2, 14);
    expect(mid(S[0][1])).toBeCloseTo(-Math.SQRT1_2, 14);
  });
  test("an eigenvalue on the negative axis or at zero is refused, exactly", () => {
    expect(() => sqrtm([[-1, 0], [0, 2]]).exact()).toThrow(/eigenvalue/);
    expect(() => sqrtm([[0, 1], [0, 0]]).exact()).toThrow(RangeError);
    expect(() => logm([[1, 2], [3, 4]]).exact()).toThrow(RangeError); // det < 0, so one eigenvalue is negative
  });
  test("a non-normal 4×4 matrix", () => {
    const A = [[9, 1, 0, 2], [0, 7, 3, 0], ["1/2", 0, 8, 1], [0, 2, 0, 10]];
    const S = sqrtm(A).exact();
    const SS = times(S, S);
    for (let i = 0; i < 4; i += 1) for (let j = 0; j < 4; j += 1) expect(contains(SS[i][j], new Rational(A[i][j]))).toBe(true);
  });
});

describe("logarithm and powers", () => {
  test("log of a shear is nilpotent", () => {
    const L = logm([[1, 1], [0, 1]]).exact();
    expect(contains(L[0][1], new Rational(1))).toBe(true);
    expect(contains(L[0][0], new Rational(0))).toBe(true);
  });
  test("log of a quarter turn is π/2 times the generator", () => {
    const L = logm(R, { bits: 100 }).exact();
    expect(mid(L[1][0])).toBeCloseTo(Math.PI / 2, 14);
    expect(width(L)).toBeLessThan(2 ** -100);
  });
  test("exp(log A) is A", () => {
    const A = [[5, 1, 0], [2, 6, 1], [0, 1, 4]];
    const L = logm(A).exact().map((row) => row.map(mid));
    const E = expm(L).quick.value;
    for (let i = 0; i < 3; i += 1) for (let j = 0; j < 3; j += 1) expect(E[i][j]).toBeCloseTo(A[i][j], 9);
  });
  test("real powers", () => {
    const P = powm([[8, 0], [0, 27]], "1/3").exact();
    expect(contains(P[0][0], new Rational(2))).toBe(true);
    expect(contains(P[1][1], new Rational(3))).toBe(true);
    const H = powm([[2, 1], [1, 3]], "-1/2").exact();
    expect(mid(H[0][0])).toBeCloseTo(0.7608452130361228, 13);
  });
});
