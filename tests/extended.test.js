// Exact results beyond the rationals, and the review fixes.
import { describe, test, expect } from "bun:test";
import * as L from "../index.js";

const A = L.algebraic;
const f = (x) => A.format(x);
const F = L.complexField;
const mul = (X, Y) => L.field.multiply(X, Y, F);
const same = (X, Y) => X.every((row, i) => row.every((x, j) => A.equals(x, Y[i][j])));
const transpose = (X) => X[0].map((_, j) => X.map((row) => row[j]));
const identity = (n) => L.algebraicMatrix(Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => (i === j ? 1 : 0))));
const hilbert = (n) => Array.from({ length: n }, (_, i) => Array.from({ length: n }, (_, j) => `1/${i + j + 1}`));

describe("review fixes", () => {
  test("the determinant preview of an ill-conditioned matrix is small, not 0", () => {
    const d = L.det(hilbert(12));
    expect(d.quick.value).toBeGreaterThan(1e-80);
    expect(d.quick.value).toBeLessThan(1e-75);
    expect(d.quick.trust).toBe("unsure");
    expect(d.exact().toNumber()).toBeCloseTo(2.6377806512535473e-78, 90);
    expect(L.det([[1, 2], [2, 4]]).quick.value).toBe(0);
  });
});

describe("eigenspaces over Q(λ)", () => {
  test("a repeated irrational eigenvalue gets a full exact basis", () => {
    const B = [[0, 2, 0, 0], [1, 0, 0, 0], [0, 0, 0, 2], [0, 0, 1, 0]];
    const ev = L.eigenvectors(B).exact();
    expect(ev.map((e) => f(e.value))).toEqual(["−√2", "√2"]);
    for (const e of ev) {
      expect(e.exact).toBe(true);
      expect(e.basis.length).toBe(2);
      // A·v = λ·v in floats for each basis vector
      const lam = e.value.toNumber();
      for (const b of e.basis) {
        const v = b.toNumbers();
        B.forEach((row, i) => expect(row.reduce((s, x, j) => s + x * v[j], 0)).toBeCloseTo(lam * v[i], 10));
      }
    }
  });
  test("a repeated complex eigenvalue too", () => {
    const R = [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]];
    const ev = L.eigenvectors(R).exact();
    expect(ev.length).toBe(2);
    expect(ev.every((e) => e.exact && e.basis.length === 2)).toBe(true);
  });
});

describe("Jordan form over algebraic numbers", () => {
  test.each([
    ["defective with ±√2", [[0, 0, 0, -4], [1, 0, 0, 0], [0, 1, 0, 4], [0, 0, 1, 0]]],
    ["complex eigenvalues", [[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]]],
    ["mixed rational and irrational", [[1, 1, 0], [1, 0, 0], [0, 0, 3]]],
    ["a rational Jordan block", [[2, 1], [0, 2]]],
  ])("%s: A·P = P·J", (_, M) => {
    const { J, P } = L.jordan(M).exact().algebraicForm();
    expect(same(mul(L.algebraicMatrix(M), P), mul(P, J))).toBe(true);
    // P is invertible
    expect(A.isZero(L.field.det(P, F))).toBe(false);
  });
  test("the block structure is read off J", () => {
    const { J } = L.jordan([[0, 0, 0, -4], [1, 0, 0, 0], [0, 1, 0, 4], [0, 0, 1, 0]]).exact().algebraicForm();
    expect(J.map((row) => row.map(f).join(" "))).toEqual(["−√2 1 0 0", "0 −√2 0 0", "0 0 √2 1", "0 0 0 √2"]);
  });
});

describe("norms and condition numbers", () => {
  const M = [[1, 2], [3, 4]];
  test("the four norms", () => {
    expect(f(L.norm(M, 1).exact())).toBe("6");
    expect(f(L.norm(M, "inf").exact())).toBe("7");
    expect(f(L.norm(M, "fro").exact())).toBe("√30");
    expect(L.norm(M, 2).exact().toNumber()).toBeCloseTo(5.464985704219043, 12);
    expect(L.norm(M, 2).quick.value).toBeCloseTo(5.464985704219043, 10);
    expect(() => L.norm(M, 3).exact()).toThrow();
  });
  test("condition numbers", () => {
    expect(f(L.cond(M, 1).exact())).toBe("21");
    expect(f(L.cond(M, 2).exact())).toBe("(15 + √221)/2");
    expect(L.cond([[1, 2], [2, 4]]).exact()).toBe(Infinity);
    expect(L.cond(hilbert(4), "inf").exact().toNumber()).toBeCloseTo(28375, 6);
  });
});

describe("symmetric matrices", () => {
  const T = [[2, -1, 0], [-1, 2, -1], [0, -1, 2]];
  test("definiteness", () => {
    expect(L.definiteness(T).exact()).toEqual({ kind: "positive-definite", positive: 3, negative: 0, zero: 0 });
    expect(L.definiteness([[1, 2], [2, 1]]).exact().kind).toBe("indefinite");
    expect(L.definiteness([[1, 1], [1, 1]]).exact()).toEqual({ kind: "positive-semidefinite", positive: 1, negative: 0, zero: 1 });
    expect(L.definiteness([[-2, 0], [0, -3]]).exact().kind).toBe("negative-definite");
    expect(() => L.definiteness([[1, 2], [3, 4]])).toThrow();
  });
  test("Cholesky: L·Lᵀ = A with exact square roots", () => {
    const { L: C } = L.cholesky([[4, 2], [2, 3]]).exact();
    expect(C.map((row) => row.map(f))).toEqual([["2", "0"], ["1", "√2"]]);
    const big = L.cholesky(T).exact().L;
    expect(same(mul(big, transpose(big)), L.algebraicMatrix(T))).toBe(true);
    expect(L.cholesky([[1, 2], [2, 1]]).exact()).toBe(null);
    expect(L.cholesky([[4, 2], [2, 3]]).quick.value.L[1][1]).toBeCloseTo(Math.SQRT2, 12);
  });
  test("orthogonal diagonalization: QᵀQ = I and A·Q = Q·D", () => {
    for (const M of [T, [[2, 1], [1, 2]], [[1, 0, 0], [0, 1, 0], [0, 0, 5]], [[2, 1, 1], [1, 2, 1], [1, 1, 2]]]) {
      const { values, Q } = L.symmetricEigen(M).exact();
      const n = M.length;
      expect(same(mul(transpose(Q), Q), identity(n))).toBe(true);
      const D = identity(n).map((row, i) => row.map((x, j) => (i === j ? values[i] : x)));
      expect(same(mul(L.algebraicMatrix(M), Q), mul(Q, D))).toBe(true);
    }
    expect(L.symmetricEigen(T).exact().values.map(f)).toEqual(["2 − √2", "2", "2 + √2"]);
  });
});

describe("integer matrices", () => {
  test("the integer kernel is a basis of all integer solutions", () => {
    expect(L.integerKernel([[1, 2, 3], [4, 5, 6]]).toStrings()).toEqual([["-1"], ["2"], ["-1"]]);
    // 2x + 4y = 0: the generator is (2, -1) up to sign, not (4, -2)
    const k = L.integerKernel([[2, 4]]).toStrings().flat();
    expect([["2", "-1"], ["-2", "1"]]).toContainEqual(k);
    expect(L.integerKernel([[1, 0], [0, 1]]).cols).toBe(0);
  });
  test("the image lattice and LLL", () => {
    expect(L.integerImage([[2, 4], [0, 6]]).toStrings()).toEqual([["2", "0"], ["0", "6"]]);
    const r = L.lllReduce([[1, 1, 1], [-1, 0, 2], [3, 5, 6]]);
    expect(L.det(r).exact().abs().toString()).toBe(L.det([[1, 1, 1], [-1, 0, 2], [3, 5, 6]]).exact().abs().toString());
    expect(() => L.lllReduce([["1/2", 1]])).toThrow();
  });
});

describe("matrices of algebraic numbers", () => {
  test("determinant, inverse and product with √2 and i", () => {
    const r2 = A.sqrt(A.fromRational(2));
    const i = A.complex(0, 1);
    const Mx = L.algebraicMatrix([[r2, 1], [1, r2]]);
    expect(f(L.field.det(Mx, F))).toBe("1");
    const inv = L.field.inverse(Mx, F);
    expect(same(mul(Mx, inv), identity(2))).toBe(true);
    const C = L.algebraicMatrix([[i, 1], [1, i]]);
    expect(f(L.field.det(C, F))).toBe("−2");
  });
});
