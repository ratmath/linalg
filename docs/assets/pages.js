// Mounts the demos on each docs page and handles page switching.
import { L, h, demo, rowReducerUI, show, matrixTable, parsePoly } from "./explorer.js";
import { workbench, workbenchFunctions } from "./workbench.js";

const M = (rows) => rows.map((r) => r.map(String));
const T = L.teach;

const DEMOS = {
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
// Page switching on the hash

const mounted = new Set();
function showPage() {
  const id = (location.hash || "#intro").slice(1);
  const page = document.querySelector(`section.page[data-page="${id}"]`) ?? document.querySelector('section.page[data-page="intro"]');
  document.querySelectorAll("section.page").forEach((s) => (s.hidden = s !== page));
  document.querySelectorAll("nav.toc a").forEach((a) => a.classList.toggle("on", a.getAttribute("href") === `#${page.dataset.page}`));
  const sel = document.getElementById("toc-select");
  if (sel) sel.value = page.dataset.page;
  for (const el of page.querySelectorAll("[data-demo]")) {
    if (mounted.has(el)) continue;
    mounted.add(el);
    try {
      DEMOS[el.dataset.demo](el);
    } catch (e) {
      el.replaceChildren(h("p", { class: "err" }, `Demo failed: ${e.message}`));
    }
  }
  if (id !== (location.hash || "#intro").slice(1)) return;
  window.scrollTo({ top: 0 });
}
window.addEventListener("hashchange", showPage);
document.getElementById("toc-select")?.addEventListener("change", (e) => (location.hash = e.target.value));
showPage();
