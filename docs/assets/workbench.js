// The exact workbench: named matrices you type rationals into, a command line
// for expressions, buttons for the library's exact operations, and a log.
import * as L from "./linalg.js";
import { h, matrixTable } from "./explorer.js";
import { evaluate, parseCell, FUNCTIONS, Poly, Rec, Note, EvalError, setContext, fromCells, cellText, FMatrix, ModP } from "./workbench-eval.js";

const { Matrix, Rational, AlgebraicReal, ComplexAlgebraic, RationalInterval } = L;
const STORE = "ratmath-linalg-workbench-v1";

const GROUPS = [
  ["Elimination", ["det", "rank", "rref", "inv", "solve", "null", "spaces"]],
  ["Factorizations", ["lu", "ldl", "qr", "gs", "pinv", "lsq"]],
  ["Structure", ["charpoly", "minpoly", "invariants", "rcf", "jordan", "eigen"]],
  ["Integer", ["hermite", "smith"]],
  ["Basics", ["trace", "T", "col", "row", "I", "zeros", "entry", "num", "den"]],
  ["Algebraic", ["sqrt", "roots", "sv", "digits"], (c) => c.domain === "rational" && c.algebraic],
  ["Ranges", ["expm", "lo", "hi", "mid", "width"], (c) => c.domain === "ranges"],
];
const GENERIC_NAMES = ["det", "rank", "rref", "inv", "solve", "null", "trace", "T", "entry", "I", "zeros"];

const EXAMPLES = {
  rational: ["A^-1", "solve(A, b)", "rref([A, b])", "A*B - 2B", "B'B", "spaces(B)", "jordan([2,0,0; 1,2,0; 0,1,2])", "eigen(A)", "1/3 + 0.#3 * 2..1/2"],
  algebraic: ["eigen(A)", "sqrt(2) + sqrt(3)", "(sqrt(2) + sqrt(3))^2", "roots(charpoly(A))", "sv(B)", "digits(sqrt(2), 60)", "[sqrt(2), 1; 1, sqrt(2)]^-1", "jordan(A)"],
  complex: ["eigen([0,-1; 1,0])", "roots(charpoly([0,-1; 1,0]))", "sqrt(-2)", "eigen(A)"],
  ranges: ["det([1:1.1, 2; 3, 4])", "inv([1:1.1, 2; 3, 4])", "solve([1:1.1, 2; 3, 4], [1; 2])", "expm([0,-1; 1,0])", "1:2 * (-1:3)", "width(det([1:1.01, 2; 3, 4]))"],
  modp: ["det(A)", "A^-1", "A*A^-1", "solve(A, b)", "rank(A)", "null([1,2; 3,6])", "1/3"],
};

const START = {
  A: [["2", "1", "-1"], ["-3", "-1", "2"], ["-2", "1", "2"]],
  B: [["1", "2"], ["1/2", "0.#3"], ["3", "-1"]],
  b: [["8"], ["-11"], ["-3"]],
};

let uid = 0;
const nid = (p) => `${p}-wb-${++uid}`;

function load() {
  try {
    const raw = localStorage.getItem(STORE);
    return raw ? JSON.parse(raw) : null;
  } catch {
    return null;
  }
}
function save(state) {
  try {
    localStorage.setItem(STORE, JSON.stringify(state));
  } catch {
    // storage is a convenience only
  }
}

// ---------------------------------------------------------------------------
// Formatting

function fmtQ(q, mode) {
  if (mode === "mixed") return q.toMixedString();
  if (mode === "decimal") {
    try {
      return q.toRepeatingDecimal().replace(/#0$/, "");
    } catch {
      return q.toString();
    }
  }
  return q.toString();
}

function polyText(p, mode) {
  const sup = (k) => String(k).split("").map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c]).join("");
  const parts = [];
  for (let i = p.c.length - 1; i >= 0; i -= 1) {
    const c = p.c[i];
    if (c.numerator === 0n) continue;
    const neg = c.numerator < 0n;
    const a = neg ? c.negate() : c;
    const one = a.numerator === 1n && a.denominator === 1n;
    const coef = fmtQ(a, mode);
    const wrapped = a.denominator === 1n || mode === "decimal" ? coef : `(${coef})`;
    const mono = i === 0 ? coef : `${one ? "" : wrapped}${p.variable}${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (neg ? " − " : " + ") : neg ? "−" : "") + mono);
  }
  return parts.join("") || "0";
}

function fmtI(x, mode) {
  const lo = fmtQ(x.low, mode);
  const hi = fmtQ(x.high, mode);
  if (lo.length + hi.length < 36) return `${lo}:${hi}`;
  const l = x.low.toNumber();
  const u = x.high.toNumber();
  const r = x.high.subtract(x.low).toNumber() / 2;
  return `${Number(((l + u) / 2).toPrecision(15))} ± ${r === 0 ? "0" : r.toExponential(1)}`;
}
function fmtAny(x, mode) {
  if (x instanceof Rational) return fmtQ(x, mode);
  if (x instanceof AlgebraicReal) return x.isRational ? fmtQ(x.toRational(), mode) : x.toString(10);
  if (x instanceof RationalInterval) return x.low.equals(x.high) ? fmtQ(x.low, mode) : fmtI(x, mode);
  if (x instanceof ModP) return String(x.v);
  if (typeof x === "bigint") return String(x);
  if (x instanceof ComplexAlgebraic) return x.toString(10);
  return String(x);
}
const polyOf = (p) => {
  const sup = (k) => String(k).split("").map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c]).join("");
  const parts = [];
  for (let i = p.length - 1; i >= 0; i -= 1) {
    const c = BigInt(p[i]);
    if (c === 0n) continue;
    const a = c < 0n ? -c : c;
    parts.push((parts.length ? (c < 0n ? " − " : " + ") : c < 0n ? "−" : "") + (i === 0 ? `${a}` : `${a === 1n ? "" : a}x${i > 1 ? sup(i) : ""}`));
  }
  return parts.join("");
};

function render(v, mode) {
  if (v instanceof FMatrix) {
    return h("div", { class: "wb-mat" }, matrixTable(v.a.map((r) => r.map((x) => fmtAny(x, mode)))), h("span", { class: "meta" }, `${v.m}×${v.n} · ${v.F.tag}`));
  }
  if (v instanceof AlgebraicReal) {
    if (v.isRational) return h("span", { class: "num wb-big" }, fmtQ(v.toRational(), mode));
    return h("span", { class: "wb-alg" }, h("span", { class: "num wb-big" }, v.toString(20)), h("small", { title: `isolated in [${v.lo[0]}/${v.lo[1]}, ${v.hi[0]}/${v.hi[1]}]` }, ` the root of ${polyOf(v.poly)} near ${v.toNumber().toPrecision(6)}`));
  }
  if (v instanceof ComplexAlgebraic) return h("span", { class: "wb-alg" }, h("span", { class: "num wb-big" }, v.toString(16)), h("small", {}, ` the only root of ${polyOf(v.poly)} in a certified box`));
  if (v instanceof RationalInterval) return h("span", { class: "num wb-big", title: `${v.low} : ${v.high}` }, fmtAny(v, mode));
  if (v instanceof ModP) return h("span", { class: "num wb-big" }, String(v.v), h("small", { class: "muted" }, ` (mod ${v.p})`));
  if (v instanceof Matrix) {
    return h("div", { class: "wb-mat" }, matrixTable(v.toArray().map((r) => r.map((q) => fmtQ(q, mode)))), h("span", { class: "meta" }, `${v.rows}×${v.cols}`));
  }
  if (v instanceof Rational) return h("span", { class: "num wb-big" }, fmtQ(v, mode));
  if (v instanceof Poly) return h("span", { class: "num" }, polyText(v, mode));
  if (v instanceof Note) return h("p", { class: "wb-note" }, v.text, v.link ? [" ", h("a", { href: v.link }, "Open that page.")] : null);
  if (v instanceof Rec) {
    const dl = h("dl", { class: "kv" });
    for (const [k, x] of v.fields) dl.append(h("dt", {}, k), h("dd", {}, render(x, mode)));
    return h("div", { class: "wb-rec" }, dl, v.note ? h("p", { class: "wb-note" }, v.note) : null);
  }
  if (Array.isArray(v)) return h("span", { class: "num" }, v.map((x) => fmtAny(x, mode)).join(",  "));
  return h("span", {}, String(v));
}

// ---------------------------------------------------------------------------
// The workbench

export function workbench(mount) {
  const saved = load();
  const state = {
    mats: saved?.mats ?? structuredClone(START),
    scalars: saved?.scalars ?? {},
    mode: saved?.mode ?? "fraction",
    numbers: saved?.numbers ?? { domain: "rational", algebraic: false, complex: false, p: "7" },
    memo: {},
    history: [],
    recall: saved?.recall ?? [],
  };
  const persist = () => save({ mats: state.mats, scalars: state.scalars, mode: state.mode, numbers: state.numbers, recall: past.slice(-50) });
  const applyNumbers = () => {
    try {
      setContext({ ...state.numbers, p: /^\d+$/.test(state.numbers.p) ? state.numbers.p : "7" });
    } catch {
      setContext({ ...state.numbers, p: "7" });
    }
  };
  applyNumbers();

  // ---- matrices
  const matsBox = h("div", { class: "wb-mats" });
  const scalarsBox = h("div", { class: "wb-scalars" });

  function cellValue(text) {
    try {
      return { q: parseCell(text) };
    } catch (e) {
      return { err: e.message };
    }
  }
  function matrixOf(name) {
    const rows = state.mats[name];
    const qs = rows.map((r) => r.map((c) => cellValue(c)));
    const bad = qs.flat().find((x) => x.err);
    if (bad) throw new EvalError(`${name} has an entry that is not a number: ${bad.err}`);
    return fromCells(qs.map((r) => r.map((x) => x.q)));
  }
  function vars() {
    const v = {};
    for (const name of Object.keys(state.mats)) v[name] = () => matrixOf(name);
    for (const [name, text] of Object.entries(state.scalars)) v[name] = () => parseCell(text);
    for (const [name, val] of Object.entries(state.memo)) v[name] = () => val;
    // resolve lazily so one bad matrix only breaks expressions that use it
    return new Proxy(v, {
      has: (t, k) => k in t,
      get: (t, k) => (typeof t[k] === "function" ? t[k]() : undefined),
    });
  }
  const nextName = () => {
    for (const n of "ABCDEFGHKMNPQRSUVWXY") if (!(n in state.mats) && !(n in state.scalars)) return n;
    return `M${Object.keys(state.mats).length + 1}`;
  };

  function matCard(name) {
    let rows = state.mats[name];
    const card = h("div", { class: "wb-card" });
    const grid = h("div", { class: "grid" });
    const area = h("textarea", { class: "wb-text", rows: "4", spellcheck: "false", "aria-label": `${name} as text`, hidden: true });
    const size = h("span", { class: "meta" });
    const msg = h("div", { class: "err small" });
    const commit = () => {
      state.mats[name] = rows;
      persist();
    };
    const validate = () => {
      const bad = [];
      grid.querySelectorAll("input").forEach((inp) => {
        const r = cellValue(inp.value);
        inp.classList.toggle("bad", !!r.err);
        inp.title = r.err ?? fmtAny(r.q, "fraction");
        if (r.err) bad.push(r.err);
      });
      msg.textContent = bad.length ? bad[0] : "";
    };
    const draw = () => {
      size.textContent = `${rows.length}×${rows[0].length}`;
      grid.innerHTML = "";
      grid.style.gridTemplateColumns = `repeat(${rows[0].length}, minmax(3.6rem, 5.2rem))`;
      rows.forEach((row, i) => row.forEach((cell, j) => {
        const inp = h("input", { value: cell, "aria-label": `${name} row ${i + 1}, column ${j + 1}`, spellcheck: "false", autocomplete: "off" });
        inp.addEventListener("input", () => {
          rows[i][j] = inp.value;
          commit();
          validate();
        });
        grid.append(inp);
      }));
      validate();
    };
    const resize = (fn) => () => {
      fn();
      commit();
      draw();
    };
    const textToggle = h("button", { class: "ghost", type: "button", onclick: () => {
      if (area.hidden) {
        area.value = rows.map((r) => r.join("  ")).join("\n");
        area.hidden = false;
        grid.hidden = true;
        textToggle.textContent = "Grid";
        area.focus();
      } else {
        area.hidden = true;
        grid.hidden = false;
        textToggle.textContent = "Text";
        draw();
      }
    } }, "Text");
    area.addEventListener("input", () => {
      const lines = area.value.split(/\n|;/).map((l) => l.trim()).filter(Boolean).map((l) => l.split(/[\s,]+/).filter(Boolean));
      if (!lines.length) return;
      const w = Math.max(...lines.map((l) => l.length));
      rows = lines.map((l) => [...l, ...Array(w - l.length).fill("0")]);
      commit();
      size.textContent = `${rows.length}×${rows[0].length}`;
      const bad = rows.flat().map(cellValue).find((x) => x.err);
      msg.textContent = bad ? bad.err : "";
    });
    const remove = h("button", { class: "ghost", type: "button", "aria-label": `Remove ${name}`, onclick: () => {
      delete state.mats[name];
      persist();
      drawMats();
    } }, "Remove");
    const head = h("div", { class: "wb-card-head" },
      h("span", { class: "wb-name" }, name), size,
      h("span", { class: "wb-spacer" }),
      textToggle, remove);
    const sizer = h("div", { class: "sizer" },
      h("button", { class: "ghost", type: "button", onclick: resize(() => rows.push(rows[0].map(() => "0"))) }, "+ row"),
      h("button", { class: "ghost", type: "button", onclick: resize(() => rows.length > 1 && rows.pop()) }, "− row"),
      h("button", { class: "ghost", type: "button", onclick: resize(() => rows.forEach((r) => r.push("0"))) }, "+ col"),
      h("button", { class: "ghost", type: "button", onclick: resize(() => rows[0].length > 1 && rows.forEach((r) => r.pop())) }, "− col"),
      h("button", { class: "ghost", type: "button", title: "Use this matrix in the command line", onclick: () => insert(name) }, `Use ${name}`),
    );
    card.append(head, grid, area, sizer, msg);
    draw();
    return card;
  }

  function drawMats() {
    matsBox.replaceChildren(...Object.keys(state.mats).map(matCard),
      h("button", { class: "ghost wb-add", type: "button", onclick: () => {
        const n = nextName();
        state.mats[n] = [["1", "0"], ["0", "1"]];
        persist();
        drawMats();
      } }, "+ New matrix"));
    const shown = (text) => {
      try {
        return fmtAny(parseCell(text), state.mode);
      } catch {
        return `${text} (not valid in this setting)`;
      }
    };
    const chip = (k, label, forget) => h("span", { class: "wb-chip" }, `${k} = ${label}`, h("button", { class: "wb-x", type: "button", "aria-label": `Forget ${k}`, onclick: () => {
      forget();
      persist();
      drawMats();
    } }, "×"));
    const chips = [
      ...Object.entries(state.scalars).map(([k, v]) => chip(k, shown(v), () => delete state.scalars[k])),
      ...Object.entries(state.memo).map(([k, v]) => chip(k, v instanceof FMatrix ? `${v.m}×${v.n} matrix (${v.F.tag})` : fmtAny(v, state.mode), () => delete state.memo[k])),
    ];
    scalarsBox.replaceChildren(...(chips.length ? [h("span", { class: "meta" }, "Saved:"), ...chips] : []));
  }

  // ---- command line
  const cmdId = nid("cmd");
  const cmd = h("input", { id: cmdId, class: "text wb-cmd", spellcheck: "false", autocomplete: "off", placeholder: "Try  A^-1   or   C = rref([A, b])" });
  const runBtn = h("button", { type: "button" }, "Run");
  const cancelBtn = h("button", { class: "ghost", type: "button", hidden: true }, "Cancel");
  const status = h("span", { class: "meta", "aria-live": "polite" });
  const log = h("div", { class: "wb-log", "aria-live": "polite" });
  let past = state.recall;
  let pastAt = past.length;
  let running = null;

  function insert(text, run = false) {
    const s = cmd.selectionStart ?? cmd.value.length;
    const e = cmd.selectionEnd ?? cmd.value.length;
    cmd.value = run ? text : cmd.value.slice(0, s) + text + cmd.value.slice(e);
    cmd.focus();
    const at = run ? cmd.value.length : s + text.length;
    cmd.setSelectionRange(at, at);
    if (run) go();
  }

  const argsFor = (name) => {
    const want = FUNCTIONS[name].args.split(",").map((s) => s.trim());
    const mats = Object.keys(state.mats);
    const square = mats.find((m) => state.mats[m].length === state.mats[m][0].length) ?? mats[0] ?? "A";
    const col = mats.find((m) => state.mats[m][0].length === 1) ?? "b";
    return want.map((a) => (a === "A" ? square : a === "b" ? col : a === "n" ? "3" : a === "r" ? "2" : a === "c" ? "3" : a === "i" || a === "j" ? "1" : a === "q" ? "6/4" : a)).join(", ");
  };

  const palette = h("div", { class: "wb-palette" });
  const drawPalette = () => palette.replaceChildren(...GROUPS.filter(([, , on]) => !on || on(state.numbers)).map(([label, names]) => h("div", { class: "wb-group" },
    h("span", { class: "wb-group-label" }, label),
    h("div", { class: "wb-buttons" }, names.map((n) => {
      const off = state.numbers.domain === "modp" && !GENERIC_NAMES.includes(n);
      return h("button", { class: off ? "ghost dim" : "ghost", type: "button", title: `${n}(${FUNCTIONS[n].args}): ${FUNCTIONS[n].help}${off ? " (rational matrices only)" : ""}`, onclick: () => insert(`${n}(${argsFor(n)})`, true) }, n);
    })),
  )));
  drawPalette();

  const examples = h("div", { class: "wb-examples" });
  const drawExamples = () => {
    const c = state.numbers;
    const list = c.domain === "rational" ? (c.complex ? EXAMPLES.complex : c.algebraic ? EXAMPLES.algebraic : EXAMPLES.rational) : EXAMPLES[c.domain];
    examples.replaceChildren(h("span", { class: "meta" }, "Examples:"), ...list.map((x) => h("button", { class: "wb-example", type: "button", onclick: () => insert(x, true) }, x)));
  };
  drawExamples();

  // ---- which numbers are allowed
  const numbersChanged = () => {
    applyNumbers();
    persist();
    drawSettings();
    drawPalette();
    drawExamples();
    drawMats();
  };
  const settings = h("div", { class: "wb-settings", role: "group", "aria-label": "Numbers allowed" });
  function drawSettings() {
    const c = state.numbers;
    const seg = (value, label, hint) => h("button", { type: "button", class: c.domain === value ? "on" : "ghost", "aria-pressed": String(c.domain === value), title: hint, onclick: () => {
      c.domain = value;
      numbersChanged();
    } }, label);
    const check = (key, label, disabled) => {
      const id = nid(key);
      const box = h("input", { type: "checkbox", id, checked: c[key] && !disabled, disabled });
      box.addEventListener("change", () => {
        c[key] = box.checked;
        if (key === "algebraic" && !box.checked) c.complex = false;
        if (key === "complex" && box.checked) c.algebraic = true;
        numbersChanged();
      });
      return h("label", { for: id, class: disabled ? "wb-check off" : "wb-check" }, box, label);
    };
    const pId = nid("p");
    const pInput = h("input", { id: pId, class: "small", value: c.p, inputmode: "numeric", "aria-label": "prime p" });
    pInput.addEventListener("change", () => {
      c.p = pInput.value.trim();
      numbersChanged();
    });
    let pNote = null;
    if (c.domain === "modp") {
      try {
        L.field.primeField(BigInt(/^\d+$/.test(c.p) ? c.p : "0"));
      } catch (e) {
        pNote = h("span", { class: "err small" }, `${c.p} is not prime; using 7`);
      }
    }
    const rational = c.domain === "rational";
    settings.replaceChildren(...[
      h("span", { class: "meta" }, "Numbers"),
      h("div", { class: "wb-seg" },
        seg("rational", "Rational", "Exact rationals, optionally with algebraic and complex numbers"),
        seg("ranges", "Ranges", "Rational intervals like 1:2, and certified enclosures"),
        seg("modp", "Mod p", "Integers modulo a prime")),
      rational ? check("algebraic", "real algebraic (√2, eigenvalues)", false) : null,
      rational ? check("complex", "complex", !c.algebraic && false) : null,
      c.domain === "modp" ? h("label", { for: pId, class: "wb-check" }, "p =", pInput) : null,
      pNote,
    ].filter(Boolean));
  }
  drawSettings();

  function entry(item) {
    const out = h("div", { class: "wb-out" });
    const drawOut = () => {
      if (item.error) out.replaceChildren(h("p", { class: "err" }, item.error));
      else out.replaceChildren(render(item.value, state.mode));
    };
    drawOut();
    item.redraw = drawOut;
    const keepCells = item.value instanceof Matrix ? item.value.toStrings() : item.value instanceof FMatrix ? item.value.a.map((row) => row.map(cellText)) : null;
    const keep = keepCells && keepCells.flat().every((c) => c !== null) && !item.name
      ? h("button", { class: "ghost", type: "button", onclick: () => {
        const n = nextName();
        state.mats[n] = keepCells;
        persist();
        drawMats();
        keep.replaceWith(h("span", { class: "meta" }, `kept as ${n}`));
      } }, "Keep as matrix")
      : null;
    return h("div", { class: "wb-entry" },
      h("div", { class: "wb-entry-head" },
        h("button", { class: "wb-src", type: "button", title: "Put back in the command line", onclick: () => insert(item.src, false) }, item.src),
        h("span", { class: "wb-spacer" }),
        item.ms !== undefined ? h("span", { class: "meta" }, item.name ? `${item.name} saved · ${fmtMs(item.ms)}` : fmtMs(item.ms)) : null,
        keep,
        h("button", { class: "ghost", type: "button", onclick: () => insert(item.src, true) }, "Again")),
      out);
  }
  const fmtMs = (ms) => (ms < 1 ? "<1 ms" : `${Math.round(ms)} ms`);
  const drawLog = () => log.replaceChildren(...state.history.slice().reverse().map(entry));

  async function go() {
    const src = cmd.value.trim();
    if (!src) return;
    if (running) running.abort();
    const ctrl = new AbortController();
    running = ctrl;
    past = [...past.filter((p) => p !== src), src];
    pastAt = past.length;
    runBtn.disabled = true;
    cancelBtn.hidden = false;
    status.textContent = "computing…";
    status.classList.add("pending");
    const t0 = performance.now();
    const item = { src };
    try {
      const r = await evaluate(src, vars(), (job) => job.settle({ signal: ctrl.signal, sliceMs: 12 }));
      item.value = r.value;
      item.ms = performance.now() - t0;
      if (r.name) {
        item.name = r.name;
        const v = r.value;
        delete state.scalars[r.name];
        delete state.memo[r.name];
        const cells = v instanceof Matrix ? v.toStrings() : v instanceof FMatrix ? v.a.map((row) => row.map(cellText)) : null;
        if (cells && cells.flat().every((c) => c !== null)) {
          state.mats[r.name] = cells;
        } else if (cells) {
          delete state.mats[r.name];
          state.memo[r.name] = v;
        } else if (cellText(v) !== null) {
          delete state.mats[r.name];
          state.scalars[r.name] = cellText(v);
        } else if (v instanceof AlgebraicReal || v instanceof ComplexAlgebraic) {
          delete state.mats[r.name];
          state.memo[r.name] = v;
        } else {
          throw new EvalError(`only a matrix or a number can be saved as ${r.name}`);
        }
        drawMats();
      }
    } catch (e) {
      if (e.name === "AbortError") {
        status.textContent = "cancelled";
        status.classList.remove("pending");
        if (running === ctrl) { running = null; runBtn.disabled = false; cancelBtn.hidden = true; }
        return;
      }
      item.error = e instanceof EvalError || e instanceof RangeError || e instanceof TypeError ? e.message : `Error: ${e.message}`;
      delete item.value;
    }
    if (running !== ctrl) return;
    running = null;
    runBtn.disabled = false;
    cancelBtn.hidden = true;
    status.classList.remove("pending");
    status.textContent = "";
    state.history.push(item);
    if (state.history.length > 50) state.history.shift();
    persist();
    drawLog();
    cmd.select();
  }

  cmd.addEventListener("keydown", (e) => {
    if (e.key === "Enter") {
      e.preventDefault();
      go();
    } else if (e.key === "ArrowUp" && past.length) {
      e.preventDefault();
      pastAt = Math.max(0, pastAt - 1);
      cmd.value = past[pastAt];
    } else if (e.key === "ArrowDown" && past.length) {
      e.preventDefault();
      pastAt = Math.min(past.length, pastAt + 1);
      cmd.value = past[pastAt] ?? "";
    }
  });
  runBtn.addEventListener("click", go);
  cancelBtn.addEventListener("click", () => running?.abort());

  // ---- number format
  const modes = [["fraction", "7/3"], ["mixed", "2..1/3"], ["decimal", "2.#3"]];
  const modeBox = h("div", { class: "wb-modes", role: "radiogroup", "aria-label": "Show numbers as" },
    h("span", { class: "meta" }, "Show numbers as"),
    modes.map(([m, label]) => {
      const b = h("button", { class: m === state.mode ? "on" : "ghost", type: "button", role: "radio", "aria-checked": String(m === state.mode), onclick: () => {
        state.mode = m;
        persist();
        modeBox.querySelectorAll("button").forEach((x) => {
          const on = x === b;
          x.className = on ? "on" : "ghost";
          x.setAttribute("aria-checked", String(on));
        });
        state.history.forEach((x) => x.redraw?.());
        drawMats();
      } }, label);
      return b;
    }));

  const reset = h("button", { class: "ghost", type: "button", onclick: () => {
    state.mats = structuredClone(START);
    state.scalars = {};
    state.history = [];
    persist();
    drawMats();
    drawLog();
  } }, "Start over");

  drawMats();
  drawLog();
  mount.replaceChildren(h("div", { class: "wb" },
    h("div", { class: "wb-top" }, settings, h("span", { class: "wb-spacer" }), reset),
    h("div", { class: "wb-top" }, modeBox),
    h("div", { class: "wb-layout" },
      h("div", { class: "wb-left" }, h("h3", {}, "Matrices"), matsBox, scalarsBox),
      h("div", { class: "wb-right" },
        h("label", { for: cmdId, class: "wb-label" }, "Expression"),
        h("div", { class: "wb-cmdrow" }, cmd, runBtn, cancelBtn),
        h("div", { class: "wb-statusrow" }, status),
        examples,
        palette,
        log)),
  ));
}

/** The syntax sheet, built from the function table so it never drifts. */
export function workbenchFunctions(mount) {
  const t = h("table", { class: "ref" }, h("tr", {}, h("th", {}, "Function"), h("th", {}, "What it gives")));
  for (const [label, names] of GROUPS) {
    t.append(h("tr", { class: "wb-sub" }, h("td", { colspan: "2" }, label)));
    for (const n of names) t.append(h("tr", {}, h("td", {}, `${n}(${FUNCTIONS[n].args})`), h("td", {}, FUNCTIONS[n].help)));
  }
  mount.replaceChildren(h("div", { class: "tablewrap" }, t));
}
