import { describe, expect, test } from "bun:test";
import { Rational, RationalInterval } from "@ratmath/core";
import { Matrix, det, rank, rref, solve, inverse, nullspace, live, Refinement, exact, field } from "../index.js";

// deterministic pseudo-random rationals
function rng(seed) {
  let s = seed >>> 0;
  return () => ((s = (s * 1664525 + 1013904223) >>> 0) / 2 ** 32);
}
function randomMatrix(rows, cols, rand, { denom = 5, range = 9 } = {}) {
  return Array.from({ length: rows }, () =>
    Array.from({ length: cols }, () => `${Math.floor(rand() * (2 * range + 1)) - range}/${1 + Math.floor(rand() * denom)}`));
}
/** Product of a random full-rank m×r and r×n matrix: rank ≤ r. */
function lowRank(m, n, r, rand) {
  return Matrix.from(randomMatrix(m, r, rand)).multiply(randomMatrix(r, n, rand));
}

// naive cofactor determinant with Core Rationals as a reference
function refDet(rows) {
  const n = rows.length;
  if (n === 1) return rows[0][0];
  let total = Rational.zero;
  for (let j = 0; j < n; j += 1) {
    const minor = rows.slice(1).map((r) => r.filter((_, k) => k !== j));
    const t = rows[0][j].multiply(refDet(minor));
    total = j % 2 ? total.subtract(t) : total.add(t);
  }
  return total;
}

describe("Matrix", () => {
  test("reads numbers, decimals, strings, bigints and Core values", () => {
    const m = Matrix.from([[1, 0.1, "3/4"], [2n, new Rational(5, 6), "-1.25"]]);
    expect(m.toStrings()).toEqual([["1", "1/10", "3/4"], ["2", "5/6", "-5/4"]]);
    expect(m.get(1, 1).equals(new Rational(5, 6))).toBe(true);
    expect(m.toNumbers()[0][1]).toBe(0.1);
  });
  test("a flat array is a column", () => {
    expect(Matrix.from([1, 2, 3]).rows).toBe(3);
  });
  test("arithmetic", () => {
    const a = Matrix.from([[1, "1/2"], ["1/3", 1]]);
    expect(a.multiply(Matrix.identity(2)).equals(a)).toBe(true);
    expect(a.add(a).equals(a.scale(2))).toBe(true);
    expect(a.subtract(a).equals(Matrix.zeros(2))).toBe(true);
    expect(a.transpose().get(0, 1).equals(new Rational(1, 3))).toBe(true);
  });
  test("huge entries convert to doubles accurately", () => {
    const m = Matrix.from([[`${10n ** 400n + 1n}/${10n ** 399n}`]]);
    expect(m.toNumbers()[0][0]).toBeCloseTo(10, 12);
  });
});

describe("exact determinant", () => {
  const rand = rng(1);
  test("Bareiss, modular and cofactor agree on random matrices", () => {
    for (let n = 1; n <= 6; n += 1) {
      for (let t = 0; t < 5; t += 1) {
        const rows = randomMatrix(n, n, rand);
        const ref = refDet(Matrix.from(rows).toArray());
        const pairs = Matrix.from(rows).pairs();
        const b = exact.run(exact.determinantBareiss(pairs));
        const m = exact.run(exact.determinantModular(pairs));
        expect(new Rational(...b).equals(ref)).toBe(true);
        expect(new Rational(...m).equals(ref)).toBe(true);
        expect(det(rows).exact().equals(ref)).toBe(true);
      }
    }
  });
  test("singular matrices give 0", () => {
    expect(det(lowRank(6, 6, 4, rand)).exact().numerator).toBe(0n);
    expect(det([[1, 2], [2, 4]]).exact().numerator).toBe(0n);
  });
  test("Hilbert matrix of order 8", () => {
    const h = Array.from({ length: 8 }, (_, i) => Array.from({ length: 8 }, (_, j) => `1/${i + j + 1}`));
    expect(det(h).exact().toString()).toBe("1/365356847125734485878112256000000");
  });
  test("preview is close and flagged", () => {
    const d = det([[2, 1], [1, 3]]);
    expect(d.quick.value).toBeCloseTo(5, 12);
    expect(d.quick.trust).toBe("good");
  });
  test("empty matrix has determinant 1", () => {
    expect(det([]).exact().equals(Rational.one)).toBe(true);
  });
});

describe("rank and rref", () => {
  const rand = rng(2);
  test("rank of low-rank products", () => {
    for (const [m, n, r] of [[5, 7, 3], [7, 5, 2], [6, 6, 6], [4, 4, 1]]) {
      const A = lowRank(m, n, r, rand);
      expect(rank(A).exact()).toBe(r);
      expect(rank(A).quick.value).toBe(r);
    }
  });
  test("rref matches a hand example", () => {
    const { matrix, pivots } = rref([[1, 2, 1, 4], [2, 4, 0, 6], [1, 2, 2, 5]]).exact();
    expect(pivots).toEqual([0, 2]);
    expect(matrix.toStrings()).toEqual([["1", "2", "0", "3"], ["0", "0", "1", "1"], ["0", "0", "0", "0"]]);
  });
  test("exact rref equals naive field rref", () => {
    for (let t = 0; t < 10; t += 1) {
      const A = lowRank(4, 6, 1 + (t % 4), rand);
      const ours = rref(A).exact();
      const ref = field.rref(A.toArray(), field.rationalField);
      expect(ours.pivots).toEqual(ref.pivots);
      expect(ours.matrix.equals(ref.rows)).toBe(true);
    }
  });
});

describe("solve", () => {
  const rand = rng(3);
  test("unique solutions, guessed and checked", () => {
    const s = solve([[2, 1, 1], [1, 3, 2], [1, 0, 0]], [4, 5, 6]);
    expect(s.quick.value.status).toBe("unique");
    const st = s.exactStage();
    expect(st.via).toBe("check");
    expect(st.value.x.toStrings()).toEqual([["6"], ["15"], ["-23"]]);
  });
  test("random systems: A·x = b exactly", () => {
    for (let n = 1; n <= 8; n += 1) {
      const A = Matrix.from(randomMatrix(n, n, rand, { denom: 7 }));
      const b = Matrix.from(randomMatrix(n, 2, rand));
      const { status, x } = solve(A, b).exact();
      if (det(A).exact().numerator === 0n) continue;
      expect(status).toBe("unique");
      expect(A.multiply(x).equals(b)).toBe(true);
    }
  });
  test("ugly answers fall back to elimination", () => {
    const h = Array.from({ length: 9 }, (_, i) => Array.from({ length: 9 }, (_, j) => `1/${i + j + 1}`));
    const s = solve(h, Array(9).fill("1/7"));
    expect(s.quick.trust).not.toBe("good");
    expect(s.exactStage().via).toBe("bareiss");
    expect(Matrix.from(h).multiply(s.exact().x).equals(Array(9).fill("1/7"))).toBe(true);
  });
  test("a wrong guess is rejected (near-singular, guess fits floats)", () => {
    // x = 1/1000003 has a denominator above the guess limit
    const s = solve([[1000003]], [1]);
    expect(s.exact().x.toStrings()).toEqual([["1/1000003"]]);
    expect(s.exactStage().via).toBe("bareiss");
  });
  test("infinitely many and no solutions", () => {
    const A = lowRank(4, 5, 3, rand);
    const x0 = Matrix.from(randomMatrix(5, 1, rand));
    const b = A.multiply(x0);
    const r = solve(A, b).exact();
    expect(r.status).toBe("infinite");
    expect(r.basis.length).toBe(2);
    expect(A.multiply(r.x).equals(b)).toBe(true);
    for (const v of r.basis) expect(A.multiply(v).equals(Matrix.zeros(4, 1))).toBe(true);

    const none = solve([[1, 1], [1, 1]], [1, 2]);
    expect(none.exact().status).toBe("none");
    expect(none.quick.value.status).toBe("none");
  });
  test("nullspace", () => {
    const basis = nullspace([[1, 2, 3], [2, 4, 6]]).exact();
    expect(basis.map((v) => v.toStrings().flat())).toEqual([["-2", "1", "0"], ["-3", "0", "1"]]);
  });
});

describe("inverse", () => {
  const rand = rng(4);
  test("A · A⁻¹ = I for random matrices", () => {
    for (let n = 1; n <= 7; n += 1) {
      const A = Matrix.from(randomMatrix(n, n, rand));
      const inv = inverse(A).exact();
      if (!inv) continue;
      expect(A.multiply(inv).equals(Matrix.identity(n))).toBe(true);
    }
  });
  test("singular gives null", () => {
    expect(inverse([[1, 2], [2, 4]]).exact()).toBeNull();
    expect(inverse(lowRank(5, 5, 4, rand)).exact()).toBeNull();
  });
});

describe("refinement", () => {
  test("stages run float, progress, exact", () => {
    const kinds = [...det(Matrix.from(randomMatrix(6, 6, rng(5)))).stages()].map((s) => s.stage);
    expect(kinds[0]).toBe("float");
    expect(kinds.at(-1)).toBe("exact");
    expect(kinds).toContain("progress");
  });
  test("settle resolves to the exact answer and reports stages", async () => {
    const seen = [];
    const v = await det([[1, 2], [3, 4]]).settle({ onStage: (s) => seen.push(s.stage) });
    expect(v.toString()).toBe("-2");
    expect(seen).toEqual(["float", "exact"]);
  });
  test("settle aborts", async () => {
    const c = new AbortController();
    const job = det(Matrix.from(randomMatrix(30, 30, rng(6))));
    const p = job.settle({ signal: c.signal, sliceMs: 0 });
    c.abort();
    await expect(p).rejects.toThrow("aborted");
    expect(job.settled).toBe(false);
  });
  test("live: quick at once, only the latest question finishes", async () => {
    const seen = [];
    const update = live((t) => solve([[1, t], [0, 1]], [1, 1]), (s) => seen.push([s.stage, s.exact ? s.value.x.toStrings().flat().join(",") : null]), { idleMs: 5 });
    update(1);
    update(2);
    update(3);
    expect(seen.map((s) => s[0])).toEqual(["float", "float", "float"]);
    await new Promise((r) => setTimeout(r, 60));
    expect(seen.filter((s) => s[0] === "exact")).toEqual([["exact", "-2,1"]]);
  });
  test("custom refinements", () => {
    const r = new Refinement(() => ({ value: 0.333, trust: "rough" }), function* () {
      yield { step: "think", done: 1, of: 1 };
      return { value: new Rational(1, 3), via: "hand" };
    });
    expect(r.quick.value).toBe(0.333);
    expect(r.exact().toString()).toBe("1/3");
  });
});

describe("generic fields", () => {
  test("interval elimination encloses the exact solution", () => {
    const I = (a, b) => new RationalInterval(a, b);
    const A = [[I("199/100", "201/100"), I(1, 1)], [I(1, 1), I(3, 3)]];
    const { rows } = field.rref(A.map((r, i) => [...r, I(i === 0 ? 4 : 5, i === 0 ? 4 : 5)]), field.intervalField, { columns: 2 });
    // exact solution at a11 = 2 is x = 7/5, y = 6/5
    expect(rows[0][2].containsValue(new Rational(7, 5))).toBe(true);
    expect(rows[1][2].containsValue(new Rational(6, 5))).toBe(true);
  });
  test("undecided pivots call back", () => {
    const I = (a, b) => new RationalInterval(a, b);
    expect(() => field.rref([[I(-1, 1), I(1, 1)]], field.intervalField, { undecided: () => { throw new Error("more digits"); } })).toThrow("more digits");
  });
});

describe("speed", () => {
  test("40×40 exact determinant and solve within a couple of seconds", () => {
    const A = Matrix.from(randomMatrix(40, 40, rng(7)));
    const t0 = performance.now();
    const d = det(A).exact();
    const t1 = performance.now();
    const { x } = solve(A, randomMatrix(40, 1, rng(8))).exact();
    const t2 = performance.now();
    expect(d.numerator).not.toBe(0n);
    expect(x.rows).toBe(40);
    // loose bounds; the point is the order of magnitude
    expect(t1 - t0).toBeLessThan(2000);
    expect(t2 - t1).toBeLessThan(4000);
    console.log(`40×40: det ${(t1 - t0).toFixed(0)} ms, solve ${(t2 - t1).toFixed(0)} ms`);
  });
});
