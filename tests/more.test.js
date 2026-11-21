import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval } from "@ratmath/core";
import {
  Matrix, det, charpoly, eigenvalues, eigenvectors, singularValues, svd, lu, ldl, qr, gramSchmidt, power, trace,
  pseudoinverse, leastSquares, spaces, hermite, smith, certifyRoot, solveIntervals, poly,
} from "../index.js";

function rng(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 2 ** 32);
}
const randInt = (rand, r = 9) => Math.floor(rand() * (2 * r + 1)) - r;
const randomInts = (m, n, rand, r) => Array.from({ length: m }, () => Array.from({ length: n }, () => randInt(rand, r)));
const randomRats = (m, n, rand) => Array.from({ length: m }, () => Array.from({ length: n }, () => `${randInt(rand)}/${1 + Math.floor(rand() * 4)}`));
const symmetric = (n, rand) => {
  const a = randomInts(n, n, rand, 6);
  return a.map((row, i) => row.map((_, j) => (i <= j ? a[i][j] : a[j][i])));
};
const permMatrix = (perm) => Matrix.from(perm.map((p) => perm.map((_, j) => (j === p ? 1 : 0))));

describe("characteristic polynomial", () => {
  const rand = rng(11);
  test("matches det(xI − A) at sample points", () => {
    for (let n = 1; n <= 6; n += 1) {
      const A = Matrix.from(randomRats(n, n, rand));
      const c = charpoly(A).exact();
      expect(c.length).toBe(n + 1);
      expect(c[n].equals(Rational.one)).toBe(true);
      for (const x of [0, 1, -2, "3/7"]) {
        const X = new Rational(String(x));
        const value = c.reduceRight((acc, k) => acc.multiply(X).add(k), Rational.zero);
        const shifted = Matrix.identity(n).scale(X).subtract(A);
        expect(value.equals(det(shifted).exact())).toBe(true);
      }
    }
  });
  test("trace and constant term", () => {
    const A = [[1, 2, 3], [4, 5, 6], [7, 8, 10]];
    const c = charpoly(A).exact();
    expect(c[2].negate().equals(trace(A))).toBe(true);
    expect(c[0].negate().equals(det(A).exact())).toBe(true);
  });
});

describe("eigenvalues", () => {
  const rand = rng(12);
  test("symmetric matrices: all real, agree with the QR preview", () => {
    for (let n = 2; n <= 7; n += 1) {
      const A = symmetric(n, rand);
      const job = eigenvalues(A);
      const { real, complexCount } = job.exact();
      expect(complexCount).toBe(0);
      const flat = real.flatMap((r) => Array(r.multiplicity).fill(r.value.toNumber()));
      const preview = job.quick.value.map((z) => z.re);
      expect(flat.length).toBe(n);
      flat.forEach((v, i) => expect(v).toBeCloseTo(preview[i], 7));
    }
  });
  test("triangular matrices have their exact diagonal as eigenvalues", () => {
    const { real } = eigenvalues([[2, 5, 7], [0, "1/3", 1], [0, 0, 2]]).exact();
    expect(real.map((r) => [r.value.toRational().toString(), r.multiplicity])).toEqual([["1/3", 1], ["2", 2]]);
  });
  test("irrational eigenvalues refine to any precision", () => {
    const { real } = eigenvalues([[1, 1], [1, 0]]).exact(); // golden ratio
    const phi = real[1].value;
    expect(phi.isRational).toBe(false);
    const [lo, hi] = phi.interval(200);
    expect(hi.subtract(lo).lessThan(new Rational(1n, 2n ** 199n))).toBe(true);
    // φ² = φ + 1 changes sign across the interval
    const g = (x) => x.multiply(x).subtract(x).subtract(Rational.one);
    expect(g(lo).multiply(g(hi)).lessThan(Rational.zero)).toBe(true);
  });
  test("rotation: complex pair counted exactly", () => {
    const e = eigenvalues([[0, -1], [1, 0]]).exact();
    expect(e.real.length).toBe(0);
    expect(e.complexCount).toBe(2);
    expect(e.complex.map((z) => Math.abs(z.im))).toEqual([1, 1]);
  });
});

describe("eigenvectors", () => {
  test("rational eigenvalues: exact eigenspaces", () => {
    const A = Matrix.from([[2, 0, 0], [0, 2, 0], [1, 0, 3]]);
    for (const e of eigenvectors(A).exact()) {
      expect(e.exact).toBe(true);
      const lam = e.value.toRational();
      for (const v of e.vectors) expect(A.multiply(v).equals(v.scale(lam))).toBe(true);
    }
    expect(eigenvectors(A).exact().map((e) => e.vectors.length)).toEqual([2, 1]);
  });
  test("irrational eigenvalues: residual shrinks as bits grow", () => {
    const A = Matrix.from([[2, 1, 0], [1, 3, 1], [0, 1, 4]]);
    const e = eigenvectors(A).exact().find((x) => !x.exact);
    const resid = (bits) => {
      const v = e.approximate(bits);
      const lam = e.value.toNumber();
      const Av = A.multiply(v).toNumbers().flat();
      return Math.max(...Av.map((x, i) => Math.abs(x - lam * v.toNumbers()[i][0])));
    };
    expect(resid(60)).toBeLessThan(1e-12);
  });
});

describe("singular values and SVD", () => {
  const rand = rng(13);
  test("agree with the Jacobi preview, include zeros", () => {
    for (const [m, n] of [[3, 2], [2, 3], [4, 4]]) {
      const A = randomInts(m, n, rand, 5);
      const job = singularValues(A);
      const flat = job.exact().flatMap((s) => Array(s.multiplicity).fill(s.value.toNumber()));
      expect(flat.length).toBe(Math.min(m, n));
      flat.forEach((v, i) => expect(v).toBeCloseTo(job.quick.value[i], 8));
    }
    expect(singularValues([[1, 2], [2, 4]]).exact().map((s) => s.value.toString())).toEqual(["5", "0"]);
  });
  test("preview reconstructs A", () => {
    const A = randomInts(4, 3, rand, 5);
    const { U, s, V } = svd(A).quick.value;
    A.forEach((row, i) => row.forEach((a, j) => {
      expect(s.reduce((acc, sk, k) => acc + U[i][k] * sk * V[j][k], 0)).toBeCloseTo(a, 9);
    }));
  });
});

describe("factorizations", () => {
  const rand = rng(14);
  test("PA = LU for square, singular and rectangular", () => {
    for (const [m, n] of [[4, 4], [3, 5], [5, 3]]) {
      for (const A of [Matrix.from(randomRats(m, n, rand)), Matrix.from(randomInts(m, 2, rand)).multiply(randomInts(2, n, rand))]) {
        const { perm, L, U } = lu(A).exact();
        expect(permMatrix(perm).multiply(A).equals(L.multiply(U))).toBe(true);
        U.toArray().forEach((row, i) => row.forEach((x, j) => j < i && expect(x.numerator).toBe(0n)));
      }
    }
  });
  test("LDLᵀ", () => {
    for (let n = 1; n <= 5; n += 1) {
      const A = Matrix.from(symmetric(n, rand));
      const r = ldl(A).exact();
      if (!r) continue;
      const D = Matrix.from(r.D.map((d, i) => r.D.map((_, j) => (i === j ? d : 0))));
      expect(r.L.multiply(D).multiply(r.L.transpose()).equals(A)).toBe(true);
    }
    expect(ldl([[0, 1], [1, 0]]).exact()).toBeNull();
    expect(ldl([[4, 2], [2, 3]]).exact().positiveDefinite).toBe(true);
  });
  test("square-root-free QR", () => {
    const A = Matrix.from(randomRats(4, 3, rand));
    const { Q, R, normsSquared } = qr(A).exact();
    expect(Q.multiply(R).equals(A)).toBe(true);
    const G = Q.transpose().multiply(Q);
    G.toArray().forEach((row, i) => row.forEach((x, j) => expect(x.equals(i === j ? normsSquared[i] : Rational.zero)).toBe(true)));
    const { Q: q, R: r } = qr(A).quick.value;
    A.toNumbers().forEach((row, i) => row.forEach((a, j) => expect(q[i].reduce((s, x, k) => s + x * r[k][j], 0)).toBeCloseTo(a, 10)));
  });
  test("Gram–Schmidt drops dependent columns", () => {
    expect(gramSchmidt([[1, 2, 1], [1, 2, 0], [0, 0, 1]]).exact().length).toBe(2);
  });
});

describe("powers", () => {
  test("Fibonacci and negative powers", () => {
    expect(power([[1, 1], [1, 0]], 90).exact().get(0, 1).toString()).toBe("2880067194370816120");
    const A = Matrix.from([[2, 1], [1, 1]]);
    expect(power(A, -3).exact().multiply(power(A, 3).exact()).equals(Matrix.identity(2))).toBe(true);
    expect(power(A, 0).exact().equals(Matrix.identity(2))).toBe(true);
    expect(power([[1, 2], [2, 4]], -1).exact()).toBeNull();
    expect(power(A, 5).quick.value[0][0]).toBeCloseTo(power(A, 5).exact().get(0, 0).toNumber(), 9);
  });
});

describe("least squares and the pseudoinverse", () => {
  const rand = rng(15);
  test("Penrose conditions on low-rank matrices", () => {
    for (let t = 0; t < 4; t += 1) {
      const A = Matrix.from(randomInts(4, 2, rand)).multiply(randomRats(2, 5, rand));
      const P = pseudoinverse(A).exact();
      expect(A.multiply(P).multiply(A).equals(A)).toBe(true);
      expect(P.multiply(A).multiply(P).equals(P)).toBe(true);
      expect(A.multiply(P).transpose().equals(A.multiply(P))).toBe(true);
      expect(P.multiply(A).transpose().equals(P.multiply(A))).toBe(true);
      const q = pseudoinverse(A).quick.value;
      P.toNumbers().forEach((row, i) => row.forEach((x, j) => expect(q[i][j]).toBeCloseTo(x, 8)));
    }
  });
  test("residual is orthogonal to the columns", () => {
    const A = Matrix.from(randomRats(6, 3, rand));
    const b = Matrix.from(randomInts(6, 1, rand));
    const { x, residual } = leastSquares(A, b).exact();
    expect(A.transpose().multiply(residual).equals(Matrix.zeros(3, 1))).toBe(true);
    expect(leastSquares(A, b).quick.value.x[0][0]).toBeCloseTo(x.get(0, 0).toNumber(), 8);
    const line = leastSquares([[1, 0], [1, 1], [1, 2]], [1, 2, 2]).exact().x;
    expect(line.toStrings()).toEqual([["7/6"], ["1/2"]]);
  });
  test("spaces", () => {
    const s = spaces([[1, 2, 3], [2, 4, 6], [1, 0, 1]]);
    expect([s.rank, s.nullity, s.columnSpace.length, s.rowSpace.length]).toEqual([2, 1, 2, 2]);
  });
});

describe("integer normal forms", () => {
  const rand = rng(16);
  test("Hermite: U·A = H, unimodular U, echelon H", () => {
    for (let t = 0; t < 5; t += 1) {
      const A = Matrix.from(randomInts(3, 4, rand));
      const { H, U, pivots } = hermite(A);
      expect(U.multiply(A).equals(H)).toBe(true);
      expect(Math.abs(det(U).exact().toNumber())).toBe(1);
      pivots.forEach((c, r) => {
        const p = H.get(r, c);
        expect(p.greaterThan(Rational.zero)).toBe(true);
        for (let i = 0; i < r; i += 1) {
          const x = H.get(i, c);
          expect(x.greaterThanOrEqual(Rational.zero) && x.lessThan(p)).toBe(true);
        }
      });
    }
  });
  test("Smith: U·A·V = S with dividing diagonal", () => {
    for (let t = 0; t < 5; t += 1) {
      const A = Matrix.from(randomInts(4, 3, rand));
      const { S, U, V, diagonal } = smith(A);
      expect(U.multiply(A).multiply(V).equals(S)).toBe(true);
      expect(Math.abs(det(U).exact().toNumber())).toBe(1);
      expect(Math.abs(det(V).exact().toNumber())).toBe(1);
      for (let i = 0; i + 1 < diagonal.length; i += 1) {
        if (diagonal[i] !== 0n) expect(diagonal[i + 1] % diagonal[i]).toBe(0n);
        else expect(diagonal[i + 1]).toBe(0n);
      }
    }
    expect(smith([[2, 4, 4], [-6, 6, 12], [10, 4, 16]]).diagonal).toEqual([2n, 2n, 156n]);
  });
  test("rejects fractions", () => {
    expect(() => smith([["1/2"]])).toThrow();
  });
});

describe("nonlinear support", () => {
  const two = new RationalInterval(2, 2);
  test("certify the intersection of a circle and a line", () => {
    const f = ([x, y]) => [x.multiply(x).add(y.multiply(y)).subtract(new Rational(4)), x.subtract(y)];
    const jacobian = ([X, Y]) => [[X.multiply(two), Y.multiply(two)], [new RationalInterval(1, 1), new RationalInterval(-1, -1)]];
    const job = certifyRoot({ f, jacobian, guess: [1, 1], bits: 80 });
    expect(job.quick.value[0]).toBeCloseTo(Math.SQRT2, 12);
    const { status, box } = job.exact();
    expect(status).toBe("unique");
    const x = box[0];
    expect(x.low.multiply(x.low).lessThan(two.low)).toBe(true);
    expect(x.high.multiply(x.high).greaterThan(two.low)).toBe(true);
    expect(x.high.subtract(x.low).lessThan(new Rational(1n, 2n ** 70n))).toBe(true);
  });
  test("enclose an interval linear system", () => {
    const { status, box } = solveIntervals([[["199/100", "201/100"], 1], [1, 3]], [4, 5]).exact();
    expect(status).toBe("unique");
    // x = 7/(3a − 1) for a in [1.99, 2.01]
    expect(box[0].containsValue(new Rational(7, 5))).toBe(true);
    expect(box[0].low.toNumber()).toBeLessThanOrEqual(7 / 5.03);
    expect(box[0].high.toNumber()).toBeGreaterThanOrEqual(7 / 4.97);
    expect(box[0].high.toNumber() - box[0].low.toNumber()).toBeLessThan(0.02);
  });
});

describe("polynomial tools", () => {
  test("square-free split and real roots", () => {
    // (x − 1)²(x² − 2)(2x + 1)
    const mul = (a, b) => {
      const o = new Array(a.length + b.length - 1).fill(0n);
      a.forEach((x, i) => b.forEach((y, j) => (o[i + j] += x * y)));
      return o;
    };
    const p = mul(mul(mul([-1n, 1n], [-1n, 1n]), [-2n, 0n, 1n]), [1n, 2n]);
    const parts = poly.squareFree(p);
    expect(parts.map((f) => f.multiplicity).sort()).toEqual([1, 2]);
    const roots = parts.flatMap((f) => poly.isolate(f.factor)).map((r) => r.toNumber()).sort((a, b) => a - b);
    [-Math.SQRT2, -0.5, 1, Math.SQRT2].forEach((v, i) => expect(roots[i]).toBeCloseTo(v, 14));
  });
  test("simplest fraction between", () => {
    expect(poly.simplestBetween([3n, 10n], [1n, 3n])).toEqual([1n, 3n]);
    expect(poly.simplestBetween([-7n, 5n], [-6n, 5n])).toEqual([-4n, 3n]);
  });
});
