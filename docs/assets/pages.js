// Mounts the demos on each docs page and handles page switching.
import { L, h, demo, rowReducerUI, show, matrixTable, parsePoly } from "./explorer.js?v=df9e75031d";
import { workbench, workbenchFunctions } from "./workbench.js?v=df9e75031d";
import { parseCell, setContext, CTX } from "./workbench-eval.js?v=df9e75031d";
import { mountRunboxes } from "./runbox.js?v=df9e75031d";

const M = (rows) => rows.map((r) => r.map(String));
const T = L.teach;

const DEMOS = {

  subspaces: (el) => demo(el, {
    title: "Four subspaces",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[1, 2, 0, 1], [2, 4, 1, 4], [3, 6, 1, 5]]) }],
    fast: ({ A }) => {
      const s = L.spaces(A);
      const left = L.nullspace(L.Matrix.from(A).transpose()).exact();
      const nul = L.nullspace(A).exact();
      const cols = (vs) => (vs.length ? vs.reduce((x, y) => x.augment(y)) : "{0}");
      return {
        rank: String(s.rank),
        "column space": cols(s.columnSpace),
        "row space": s.rowSpace.length ? L.Matrix.from(s.rowSpace.flatMap((r) => r.toArray())) : "{0}",
        "null space": cols(nul),
        "left null space": cols(left),
      };
    },
    teach: ({ A }) => T.fourSubspaces(A),
  }),

  sturm: (el) => demo(el, {
    title: "Counting and isolating real roots",
    inputs: [{ id: "p", kind: "poly", label: "p (coefficients, highest first)", value: "1, 0, -4, 1" }],
    fast: ({ p }) => {
      const roots = L.algebraic.realRoots(p.map((c) => { const q = new L.Rational(String(c).trim()); return [q.numerator, q.denominator]; }));
      return { "real roots": roots.length ? roots.map((r) => `${r.value.isRational ? r.value.toRational() : r.value.toString(16)}${r.multiplicity > 1 ? ` (×${r.multiplicity})` : ""}`) : "none" };
    },
    teach: ({ p }) => T.sturmWalk(integerPoly(p)),
  }),

  dixonSteps: (el) => demo(el, {
    title: "Dixon lifting, digit by digit",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[3, 1, 2], [1, 4, 1], [2, 1, 5]]) }, { id: "b", kind: "vector", label: "b", value: "1, 2, 3" }],
    fast: ({ A, b }) => L.solve(A, b),
    teach: ({ A, b }) => T.dixonLifting(A, b),
  }),

  dixonRace(el) {
    const id = "race-n";
    const input = h("input", { type: "range", id, min: "20", max: "160", step: "10", value: "80" });
    const label = h("output", { for: id }, "");
    const out = h("div", { class: "race" });
    const btn = h("button", { type: "button" }, "Race");
    const setLabel = () => (label.textContent = `${input.value}×${input.value}`);
    input.addEventListener("input", setLabel);
    setLabel();
    btn.addEventListener("click", async () => {
      const n = Number(input.value);
      btn.disabled = true;
      let seed = n * 7919;
      const rnd = () => ((seed = (seed * 16807) % 2147483647) % 199) - 99;
      const A = Array.from({ length: n }, () => Array.from({ length: n }, rnd));
      const b = Array.from({ length: n }, rnd);
      out.replaceChildren(h("p", { class: "muted pending" }, "solving…"));
      await new Promise((r) => setTimeout(r, 30));
      const t0 = performance.now();
      const job = L.solve(A, b);
      const stage = job.exactStage();
      const t1 = performance.now();
      const M2 = L.Matrix.from(A);
      L.exact.run(L.exact.solvePairs(M2.pairs(), L.Matrix.from(b).pairs()));
      const t2 = performance.now();
      const x0 = stage.value.x.get(0, 0);
      out.replaceChildren(h("table", { class: "data" },
        h("tr", {}, h("th", {}, "method"), h("th", {}, "time")),
        h("tr", {}, h("td", {}, `solve() → ${stage.via === "dixon" ? "Dixon lifting" : stage.via}`), h("td", {}, `${Math.round(t1 - t0)} ms`)),
        h("tr", {}, h("td", {}, "fraction-free elimination"), h("td", {}, `${Math.round(t2 - t1)} ms`))),
        h("p", { class: "small muted" }, `x₁ has a ${x0.denominator.toString().length}-digit denominator. Both answers are exact and agree.`));
      btn.disabled = false;
    });
    el.replaceChildren(h("div", { class: "demo-card" }, h("div", { class: "demo-title" }, "Race: lifting against elimination"),
      h("label", { for: id, class: "slider" }, "Random integer system, entries −99…99: ", label), input, h("div", { class: "actions" }, btn), out));
  },

  worker(el) {
    const id = "wk-n";
    const input = h("input", { type: "range", id, min: "40", max: "200", step: "20", value: "140" });
    const label = h("output", { for: id }, "");
    const setLabel = () => (label.textContent = `${input.value}×${input.value}`);
    input.addEventListener("input", setLabel);
    setLabel();
    const dot = h("div", { class: "wk-dot" });
    const track = h("div", { class: "wk-track" }, dot);
    const out = h("div", {});
    let t = 0;
    let longest = 0;
    let last = performance.now();
    const tick = (now) => {
      longest = Math.max(longest, now - last);
      last = now;
      t += 1;
      dot.style.transform = `translateX(${(Math.sin(t / 20) + 1) * 45}%)`;
      if (el.isConnected) requestAnimationFrame(tick);
    };
    requestAnimationFrame(tick);
    let client = null;
    const getClient = () => {
      if (!client) client = L.workerClient(new Worker(new URL("./linalg-worker.js", import.meta.url), { type: "module" }));
      return client;
    };
    const system = () => {
      const n = Number(input.value);
      let seed = n * 104729;
      const rnd = () => ((seed = (seed * 16807) % 2147483647) % 199) - 99;
      return [Array.from({ length: n }, () => Array.from({ length: n }, rnd)), Array.from({ length: n }, rnd)];
    };
    const runWith = (how) => async () => {
      const [A, b] = system();
      out.replaceChildren(h("p", { class: "muted pending" }, "computing…"));
      await new Promise((r) => setTimeout(r, 50));
      longest = 0;
      last = performance.now();
      const t0 = performance.now();
      try {
        if (how === "sync") L.solve(A, b).exact();
        else if (how === "slices") await L.solve(A, b).settle({ sliceMs: 8 });
        else await getClient().solve(A, b).settle();
      } catch (e) {
        out.replaceChildren(h("p", { class: "err" }, `Failed: ${e.message}`));
        return;
      }
      const ms = performance.now() - t0;
      await new Promise((r) => requestAnimationFrame(() => r()));
      const pause = Math.max(longest, performance.now() - last);
      out.replaceChildren(h("table", { class: "data" },
        h("tr", {}, h("th", {}, "where"), h("th", {}, "time"), h("th", {}, "longest freeze")),
        h("tr", {}, h("td", {}, how === "sync" ? "exact() on the page" : how === "slices" ? "settle() in 8 ms slices" : "in a Worker"), h("td", {}, `${Math.round(ms)} ms`), h("td", {}, `${Math.round(pause)} ms`))),
        h("p", { class: "small muted" }, "Watch the dot: it stops while the page's own thread is busy."));
    };
    el.replaceChildren(h("div", { class: "demo-card" }, h("div", { class: "demo-title" }, "Where the exact work runs"),
      h("label", { for: id, class: "slider" }, "Exact solve of a random integer system: ", label), input, track,
      h("div", { class: "actions" },
        h("button", { type: "button", onclick: runWith("sync") }, "On the page"),
        h("button", { type: "button", onclick: runWith("slices") }, "In slices"),
        h("button", { type: "button", onclick: runWith("worker") }, "In a worker")),
      out));
  },

  fieldModp: (el) => demo(el, {
    title: "Integers mod p",
    inputs: [{ id: "A", kind: "grid", label: "A (rationals map into F_p)", value: M([[1, 2, 3], [4, 5, 6], [7, 8, 10]]) }, { id: "p", kind: "text", label: "prime p", value: "7" }],
    fast: ({ A, p }) => {
      const F = L.field.primeField(BigInt(p.trim()));
      const rows = A.map((r) => r.map((c) => F.from(new L.Rational(c))));
      const r = L.field.rref(rows, F);
      const inv = rows.length === rows[0].length ? L.field.inverse(rows, F) : null;
      const str = (m) => L.Matrix.from(m.map((row) => row.map((x) => String(x))));
      return {
        "A mod p": str(rows),
        det: rows.length === rows[0].length ? String(L.field.det(rows, F)) : "not square",
        rank: String(r.pivots.length),
        "A⁻¹ mod p": inv ? str(inv) : "singular mod p",
      };
    },
  }),

  fieldRanges: (el) => demo(el, {
    title: "Ranges",
    inputs: [{ id: "A", kind: "grid", label: "A (entries like 1:1.1 or 3)", value: [["2:2.1", "1"], ["1", "3:3.05"]] }],
    fast: ({ A }) => withCells({ domain: "ranges" }, () => {
      const rows = A.map((r) => r.map((c) => {
        const v = parseCell(c);
        return v instanceof L.RationalInterval ? v : new L.RationalInterval(v, v);
      }));
      const F = L.field.intervalField;
      const fmt = (x) => (x.low.equals(x.high) ? x.low.toString() : `${x.low.toDecimal(6)} … ${x.high.toDecimal(6)}`);
      let inv = null;
      let why = null;
      try {
        inv = L.field.inverse(rows, F, { undecided: () => { throw new Error("a pivot range contains 0"); } });
      } catch (e) {
        why = e.message;
      }
      return {
        det: fmt(L.field.det(rows, F)),
        "A⁻¹ encloses": inv ? inv.map((r) => r.map(fmt).join("    ")) : `undecided: ${why}`,
      };
    }),
  }),

  fieldAlg: (el) => demo(el, {
    title: "Real algebraic numbers",
    inputs: [{ id: "A", kind: "grid", label: "A (entries like 2, 1/2, sqrt(2), 1 + sqrt(3))", value: [["sqrt(2)", "1"], ["1", "sqrt(3)"]] }],
    fast: ({ A }) => withCells({ domain: "rational", algebraic: true }, () => {
      const F = L.field.algebraicField;
      const rows = A.map((r) => r.map((c) => {
        const v = parseCell(c);
        return v instanceof L.AlgebraicReal ? v : L.algebraic.fromRational([v.numerator, v.denominator]);
      }));
      const d = L.field.det(rows, F);
      const inv = L.field.inverse(rows, F);
      const fmt = (x) => (x.isRational ? String(x.toRational()) : `${x.toString(12)} (root of degree ${x.poly.length - 1})`);
      return {
        det: fmt(d),
        "A⁻¹": inv ? inv.map((r) => r.map(fmt).join("    ")) : "singular",
      };
    }),
  }),
  workbench,
  workbenchFunctions,
  hero(el) {
    // a slider moves one entry; the preview tracks it, the exact answer settles
    const id = "hero-t";
    const out = h("div", { class: "hero-out" });
    const input = h("input", { type: "range", id, min: "-40", max: "40", value: "7", step: "1" });
    const label = h("output", { for: id }, "");
    const preview = h("div", { class: "stage float" });
    const exactBox = h("div", { class: "stage exact" });
    out.append(preview, exactBox);
    const matFor = (t) => [[2, `${t}/10`, 1], [1, 3, 2], [1, 0, 1]];
    const update = L.live((t) => L.solve(matFor(t), [1, 2, 3]), (s) => {
      const box = s.exact ? exactBox : preview;
      const head = h("div", { class: "stage-head" }, h("span", { class: `chip ${s.exact ? "exact" : "float"}` }, s.exact ? "exact" : "preview"), s.exact ? h("span", { class: "meta" }, s.via) : h("span", { class: `trust ${s.trust}` }, s.trust));
      const x = s.value.x;
      const vals = s.exact ? x.toStrings().map((r) => r[0]) : x.map((r) => Number(r[0].toPrecision(12)));
      box.replaceChildren(head, h("div", { class: "xs" }, vals.map((v, i) => h("div", {}, h("span", { class: "muted" }, `x${"₁₂₃"[i]} = `), h("span", { class: "num" }, v)))));
      if (!s.exact) exactBox.classList.add("stale");
      else exactBox.classList.remove("stale");
    }, { idleMs: 120 });
    const mat = h("div", {});
    const go = () => {
      const t = Number(input.value);
      label.textContent = `a₁₂ = ${t}/10`;
      mat.replaceChildren(matrixTable(M(matFor(t)).map((r, i) => [...r, String([1, 2, 3][i])]), { bar: 3, cols: [1] }));
      update(t);
    };
    input.addEventListener("input", go);
    el.replaceChildren(h("div", { class: "hero-demo" }, h("div", { class: "hero-left" }, h("label", { for: id, class: "slider" }, "Drag to change one entry ", label), input, mat), out));
    go();
  },

  rref: (el) => demo(el, {
    title: "Reduced row echelon form",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[1, 2, 1, 4], [2, 4, 0, 6], [1, 2, 2, 5]]) }],
    fast: ({ A }) => L.rref(A),
    teach: ({ A }) => T.gaussJordan(A),
  }),
  solve: (el) => demo(el, {
    title: "Solve A·x = b",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]]) }, { id: "b", kind: "vector", label: "b", value: "8, -11, -3" }],
    fast: ({ A, b }) => L.solve(A, b),
    teach: ({ A, b }) => T.solveSystem(A, b),
  }),
  rowreducer: (el) => rowReducerUI(el, M([[0, 2, 4], [1, 1, 1], [2, 4, 7]])),

  detCofactor: (el) => demo(el, {
    title: "Cofactor expansion",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[3, 1, 4], [1, 5, 9], [2, 6, 5]]) }],
    fast: ({ A }) => L.det(A),
    teach: ({ A }) => T.cofactorExpansion(A),
  }),
  detElim: (el) => demo(el, {
    title: "Elimination",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[0, 2, 1], [1, 1, 1], [2, 1, 3]]) }],
    fast: ({ A }) => L.det(A),
    teach: ({ A }) => T.determinantByElimination(A),
  }),
  detBareiss: (el) => demo(el, {
    title: "Fraction-free (Bareiss)",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, "1/2", 1, 3], [1, 1, 1, 0], [2, 1, 3, 1], [1, -1, 0, 2]]) }],
    fast: ({ A }) => L.det(A),
    teach: ({ A }) => T.bareiss(A),
  }),
  detModular: (el) => demo(el, {
    title: "Small primes and the Chinese remainder theorem",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[3, 1, 4, 1], [5, 9, 2, 6], [5, 3, 5, 8], [9, 7, 9, 3]]) }],
    fast: ({ A }) => L.det(A),
    teach: ({ A }) => T.modularDeterminant(A),
  }),

  inverse: (el) => demo(el, {
    title: "Inverse by [A | I]",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 1, 1], [1, 3, 2], [1, 0, 0]]) }],
    fast: ({ A }) => L.inverse(A),
    teach: ({ A }) => T.inverse(A),
  }),
  guess: (el) => demo(el, {
    title: "Guess from floats, check exactly",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[4, -2, 1], [3, 6, -4], [2, 1, 8]]) }, { id: "b", kind: "vector", label: "b", value: "12, -25, 32" }],
    fast: ({ A, b }) => L.solve(A, b),
    teach: ({ A, b }) => T.guessAndCheck(A, b),
  }),

  eigen: (el) => demo(el, {
    title: "Eigenvalues",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 1, 0], [1, 3, 1], [0, 1, 4]]) }],
    fast: ({ A }) => L.eigenvalues(A),
    teach: ({ A }) => T.eigen(A),
    render: (v, kind) => (kind === "float" ? show((v ?? []).map((z) => (z.im ? `${Number(z.re.toPrecision(10))} ${z.im < 0 ? "−" : "+"} ${Number(Math.abs(z.im).toPrecision(10))}i` : Number(z.re.toPrecision(12))))) : show({
      real: v.real.map((r) => `${r.value.toString(18)}${r.multiplicity > 1 ? ` (×${r.multiplicity})` : ""}`),
      complex: v.complex.map((c) => (c.value ? c.value.toString(14) : `${c.approx.re} + ${c.approx.im}i (float)`)),
      "characteristic polynomial": polyFromRationals(v.charpoly),
    })),
  }),
  refine(el) {
    const id = "bits";
    const input = h("input", { type: "range", id, min: "4", max: "400", value: "40" });
    const out = h("div", {});
    const ev = L.eigenvalues([[1, 1], [1, 0]]).exact().real[1].value; // golden ratio
    const go = () => {
      const bits = Number(input.value);
      const [lo, hi] = ev.interval(bits);
      out.replaceChildren(
        h("p", {}, `Interval of width below 2^−${bits}, for the root of ${"x² − x − 1"} between 1 and 2:`),
        h("div", { class: "lines" }, h("div", {}, `low  ${lo.toDecimal ? lo.toDecimal(Math.ceil(bits * 0.302) + 2) : lo.toString()}`), h("div", {}, `high ${hi.toDecimal ? hi.toDecimal(Math.ceil(bits * 0.302) + 2) : hi.toString()}`)),
        h("p", { class: "muted small" }, `As exact fractions: ${lo.toString().length + hi.toString().length} characters for both endpoints.`),
      );
    };
    input.addEventListener("input", go);
    el.replaceChildren(h("div", { class: "demo-card" }, h("div", { class: "demo-title" }, "Sharpen the golden ratio"), h("label", { for: id, class: "slider" }, "Bits of precision"), input, out));
    go();
  },
  eigvec: (el) => demo(el, {
    title: "Eigenvectors, exact over Q(λ)",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 1, 0], [1, 3, 1], [0, 1, 4]]) }],
    fast: ({ A }) => L.eigenvectors(A),
    render: (v, kind) => (kind === "float"
      ? show((v ?? []).map((e) => `λ ≈ ${Number(e.value.toPrecision(10))}:  (${e.vector.map((x) => Number(x.toPrecision(8))).join(", ")})`))
      : show(v.map((e) => `λ = ${e.value.toString(12)}:  ${e.algebraic ? `(${e.algebraic.toStrings().join(",  ")})` : e.vectors.map((m) => `(${m.toStrings().map((r) => r[0]).join(", ")})`).join("  and  ")}`))),
  }),

  jordan: (el) => demo(el, {
    title: "Jordan and rational canonical forms",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[5, 4, 2, 1], [0, 1, -1, -1], [-1, -1, 3, 0], [1, 1, -1, 2]]) }],
    fast: ({ A }) => L.similarity(A),
    render: (v, kind) => (kind === "float" ? show((v.eigenvalues ?? []).map((z) => Number(z.re.toPrecision(10)))) : show({
      "invariant factors": v.invariantFactors.map(polyFromRationals),
      "minimal polynomial": polyFromRationals(v.minimalPolynomial),
      "Jordan blocks": v.jordan.map((j) => `λ = ${j.value.toString(10)}: sizes ${j.blocks.join(", ")}`),
      J: v.J,
      P: v.P,
      "rational canonical form": v.rationalCanonical,
    })),
  }),
  jordanSteps: (el) => demo(el, {
    title: "Blocks from kernel dimensions",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 0, 0, 0], [1, 2, 0, 0], [0, 1, 2, 0], [0, 0, 0, 2]]) }, { id: "lambda", kind: "text", label: "λ", value: "2" }],
    teach: ({ A, lambda }) => T.jordanBlocks(A, lambda.trim()),
    teachLabel: "Step by step",
  }),

  gram: (el) => demo(el, {
    title: "Gram–Schmidt and QR",
    inputs: [{ id: "A", kind: "matrix", label: "A (columns are the vectors)", value: M([[1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 1, 1]]) }],
    fast: ({ A }) => L.qr(A),
    teach: ({ A }) => T.gramSchmidt(A),
  }),
  lsq: (el) => demo(el, {
    title: "Fit a line: least squares",
    inputs: [{ id: "A", kind: "matrix", label: "A (1, xᵢ)", value: M([[1, 0], [1, 1], [1, 2], [1, 3]]) }, { id: "b", kind: "vector", label: "b (yᵢ)", value: "1, 2, 2, 4" }],
    fast: ({ A, b }) => L.leastSquares(A, b),
    teach: ({ A, b }) => T.leastSquares(A, b),
  }),
  pinv: (el) => demo(el, {
    title: "Pseudoinverse",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[1, 2], [2, 4], [3, 6]]) }],
    fast: ({ A }) => L.pseudoinverse(A),
  }),
  svd: (el) => demo(el, {
    title: "Singular values",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[3, 0], [4, 5]]) }],
    fast: ({ A }) => L.singularValues(A),
    render: (v, kind) => (kind === "float" ? show(v) : show(v.map((s) => `${s.value.toString(16)}${s.multiplicity > 1 ? ` (×${s.multiplicity})` : ""}`))),
  }),

  lu: (el) => demo(el, {
    title: "PA = LU",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[2, 1, 1], [4, -6, 0], [-2, 7, 2]]) }],
    fast: ({ A }) => L.lu(A),
    teach: ({ A }) => T.lu(A),
  }),
  ldl: (el) => demo(el, {
    title: "A = L·D·Lᵀ",
    inputs: [{ id: "A", kind: "matrix", label: "A (symmetric)", value: M([[4, 2, 2], [2, 5, 3], [2, 3, 6]]) }],
    fast: ({ A }) => L.ldl(A),
  }),
  power: (el) => demo(el, {
    title: "Powers by squaring",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[1, 1], [1, 0]]) }, { id: "k", kind: "text", label: "k", value: "30" }],
    fast: ({ A, k }) => L.power(A, BigInt(k.trim() || "0")),
    teach: ({ A, k }) => T.powerBySquaring(A, BigInt(k.trim() || "0")),
  }),

  integer: (el) => demo(el, {
    title: "Hermite and Smith normal forms",
    inputs: [{ id: "A", kind: "matrix", label: "A (integers)", value: M([[2, 4, 4], [-6, 6, 12], [10, 4, 16]]) }],
    fast: ({ A }) => {
      const H = L.hermite(A);
      const S = L.smith(A);
      return { "Hermite H": H.H, "U (U·A = H)": H.U, "Smith S": S.S, "invariant factors": S.diagonal.map(String) };
    },
  }),

  expm: (el) => demo(el, {
    title: "exp(A), certified",
    inputs: [{ id: "A", kind: "matrix", label: "A", value: M([[0, 1], [-1, 0]]) }],
    fast: ({ A }) => L.expm(A, { bits: 60 }),
    teach: ({ A }) => T.exponentialSeries(A),
  }),

  newton: (el) => demo(el, {
    title: "Newton, then a certificate",
    inputs: [{ id: "p", kind: "poly", label: "p (coefficients, highest first)", value: "1, 0, -2" }, { id: "x0", kind: "text", label: "start", value: "1" }],
    teach: ({ p, x0 }) => T.newtonKrawczyk(p, { x0: x0.trim() || "1" }),
  }),
  certify: (el) => demo(el, {
    title: "Circle meets line",
    inputs: [{ id: "r", kind: "text", label: "x² + y² =", value: "4" }, { id: "s", kind: "text", label: "y = slope·x, slope", value: "1/2" }],
    fast: ({ r, s }) => {
      const { Rational, RationalInterval } = coreFrom();
      const R2 = new Rational(r.trim());
      const S = new Rational(s.trim());
      const two = new RationalInterval(2, 2);
      return L.certifyRoot({
        f: ([x, y]) => [x.multiply(x).add(y.multiply(y)).subtract(R2), y.subtract(S.multiply(x))],
        jacobian: ([X, Y]) => [[X.multiply(two), Y.multiply(two)], [new RationalInterval(S.negate(), S.negate()), new RationalInterval(1, 1)]],
        guess: [1, 1],
        bits: 60,
      });
    },
  }),
  intervals: (el) => demo(el, {
    title: "A system with uncertain coefficients",
    inputs: [{ id: "A", kind: "grid", label: "A (write a range as lo:hi)", value: M([["1.99:2.01", 1], [1, "2.95:3.05"]]) }, { id: "b", kind: "vector", label: "b", value: "4, 5" }],
    fast: ({ A, b }) => L.solveIntervals(A.map((r) => r.map(rangeOf)), b.map(rangeOf)),
  }),
};

function rangeOf(s) {
  const t = String(s).trim();
  if (t.includes(":")) return t.split(":").map((x) => x.trim());
  return t;
}
function integerPoly(p) {
  const qs = p.map((c) => new L.Rational(String(c).trim()));
  const g = (a, b) => (b ? g(b, a % b) : a);
  const lcm = qs.reduce((l, q) => (l * q.denominator) / g(l, q.denominator), 1n);
  const ints = qs.map((q) => q.numerator * (lcm / q.denominator));
  while (ints.length > 1 && ints[ints.length - 1] === 0n) ints.pop();
  return ints;
}
/** Parse cells under a number setting without disturbing the workbench's. */
function withCells(ctx, fn) {
  const saved = { ...CTX };
  setContext(ctx);
  try {
    return fn();
  } finally {
    setContext(saved);
  }
}
function coreFrom() {
  const one = L.det([[1]]).exact();
  const Rational = one.constructor;
  const RationalInterval = L.expm([[0]]).exact()[0][0].constructor;
  return { Rational, RationalInterval };
}
function polyFromRationals(cs) {
  const sup = (k) => String(k).split("").map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c]).join("");
  const parts = [];
  for (let i = cs.length - 1; i >= 0; i -= 1) {
    const c = cs[i];
    if (c.numerator === 0n) continue;
    const neg = c.numerator < 0n;
    const a = neg ? c.negate() : c;
    const one = a.numerator === 1n && a.denominator === 1n;
    const mono = i === 0 ? a.toString() : `${one ? "" : a.toString()}x${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (neg ? " − " : " + ") : neg ? "−" : "") + mono);
  }
  return parts.join("");
}

// ---------------------------------------------------------------------------
// Each page is its own HTML file; mount the demos on this one.

// Old single-page links (index.html#eigen) go to the page's own file.
const old = location.hash.slice(1);
const here = document.querySelector("section.page")?.dataset.page;
if (here === "intro" && old && old !== "intro" && document.querySelector(`#toc-select option[value="${old}"]`)) location.replace(`${old}.html`);

for (const el of document.querySelectorAll("[data-demo]")) {
  try {
    DEMOS[el.dataset.demo](el);
  } catch (e) {
    el.replaceChildren(h("p", { class: "err" }, `Demo failed: ${e.message}`));
  }
}
document.getElementById("toc-select")?.addEventListener("change", (e) => {
  location.href = e.target.value === "intro" ? "index.html" : `${e.target.value}.html`;
});
mountRunboxes(document, { scope: L });
