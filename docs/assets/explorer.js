// UI pieces for the docs: matrix editors, the fast panel (preview then exact),
// the step-by-step player, and the hands-on row reducer.
import * as L from "./linalg.js?v=df9e75031d";

const h = (tag, attrs = {}, ...kids) => {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (v === null || v === undefined || v === false) continue;
    if (k === "class") el.className = v;
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else el.setAttribute(k, v === true ? "" : v);
  }
  for (const kid of kids.flat()) if (kid !== null && kid !== undefined && kid !== false) el.append(kid instanceof Node ? kid : document.createTextNode(String(kid)));
  return el;
};
export { h };

let uid = 0;
const nextId = (p) => `${p}-${++uid}`;

// ---------------------------------------------------------------------------
// Rendering values

export function matrixTable(rows, { bar = null, pivot = null, hl = [], cols = [], cls = "" } = {}) {
  const t = h("table", { class: `mx ${cls}` });
  rows.forEach((row, i) => {
    const tr = h("tr", { class: hl.includes(i) ? "hl" : null });
    row.forEach((cell, j) => {
      const c = [];
      if (bar !== null && j === bar) c.push("bar");
      if (pivot && pivot[0] === i && pivot[1] === j) c.push("pivot");
      if (cols.includes(j)) c.push("hlc");
      tr.append(h("td", { class: c.join(" ") || null }, cell));
    });
    t.append(tr);
  });
  return h("div", { class: "mxwrap" }, t);
}

const isMatrix = (v) => v instanceof L.Matrix;
const num = (x) => (Math.abs(x) < 1e-12 && x !== 0 ? x.toExponential(3) : Number(x.toPrecision(10)).toString());

/** Best-effort display of any library value. */
export function show(v, depth = 0) {
  if (v === null || v === undefined) return h("span", { class: "muted" }, "none");
  if (typeof v === "number") return h("span", { class: "num" }, num(v));
  if (typeof v === "string" || typeof v === "boolean") return h("span", {}, String(v));
  if (typeof v === "bigint") return h("span", { class: "num" }, String(v));
  if (isMatrix(v)) return matrixTable(v.toStrings());
  if (v instanceof L.AlgebraicReal) {
    return h("span", { class: "alg", title: v.isRational ? "exact rational" : `root of ${polyText(v.poly)} in an isolating interval` }, v.toString(16), v.isRational ? "" : h("small", {}, ` root of ${polyText(v.poly)}`));
  }
  if (v instanceof L.ComplexAlgebraic) return h("span", { class: "alg" }, v.toString(12), h("small", {}, ` root of ${polyText(v.poly)}`));
  if (v instanceof L.AlgebraicVector) return h("span", { class: "alg" }, `(${v.toStrings().join(",  ")})`);
  if (v && typeof v.numerator === "bigint") return h("span", { class: "num" }, v.toString());
  if (v && v.low && v.high) {
    const lo = v.low.toNumber();
    const hi = v.high.toNumber();
    return h("span", { class: "num", title: `[${v.low.toString()}, ${v.high.toString()}]` }, `${num((lo + hi) / 2)} ± ${((hi - lo) / 2).toExponential(1)}`);
  }
  if (Array.isArray(v)) {
    if (v.length && Array.isArray(v[0]) && v[0].every((x) => typeof x === "number")) return matrixTable(v.map((r) => r.map(num)));
    if (v.length && Array.isArray(v[0]) && v[0].every((x) => x && x.low)) return matrixTable(v.map((r) => r.map((x) => show(x).textContent)));
    if (v.every((x) => typeof x === "number")) return h("span", { class: "num" }, v.map(num).join(",  "));
    if (v.every((x) => x && typeof x.numerator === "bigint")) return h("span", { class: "num" }, v.map(String).join(",  "));
    return h("div", { class: "list" }, v.map((x) => h("div", {}, show(x, depth + 1))));
  }
  if (typeof v === "object") {
    const dl = h("dl", { class: "kv" });
    for (const [k, x] of Object.entries(v)) {
      if (typeof x === "function") continue;
      dl.append(h("dt", {}, k), h("dd", {}, show(x, depth + 1)));
    }
    return dl;
  }
  return h("span", {}, String(v));
}

export function polyText(p) {
  const sup = (k) => String(k).split("").map((c) => "⁰¹²³⁴⁵⁶⁷⁸⁹"[c]).join("");
  const parts = [];
  for (let i = p.length - 1; i >= 0; i -= 1) {
    const c = BigInt(p[i]);
    if (c === 0n) continue;
    const a = c < 0n ? -c : c;
    const mono = i === 0 ? `${a}` : `${a === 1n ? "" : a}x${i > 1 ? sup(i) : ""}`;
    parts.push((parts.length ? (c < 0n ? " − " : " + ") : c < 0n ? "−" : "") + mono);
  }
  return parts.join("") || "0";
}

// ---------------------------------------------------------------------------
// Inputs

/** Editable grid. value: rows of strings. */
export function matrixEditor(value, { onChange, resizable = true, bar = null } = {}) {
  let rows = value.map((r) => r.map(String));
  const wrap = h("div", { class: "editor" });
  const grid = h("div", { class: "grid" });
  const fire = () => onChange?.(rows);
  const draw = () => {
    grid.innerHTML = "";
    grid.style.gridTemplateColumns = `repeat(${rows[0].length}, minmax(3.2rem, 4.6rem))`;
    rows.forEach((row, i) => row.forEach((cell, j) => {
      const id = nextId("cell");
      const inp = h("input", { id, value: cell, "aria-label": `row ${i + 1}, column ${j + 1}`, class: bar !== null && j === bar ? "barl" : null, inputmode: "text", spellcheck: "false" });
      inp.addEventListener("input", () => {
        rows[i][j] = inp.value;
        fire();
      });
      grid.append(inp);
    }));
  };
  draw();
  wrap.append(grid);
  if (resizable) {
    const btn = (label, fn) => h("button", { class: "ghost", type: "button", onclick: () => { fn(); draw(); fire(); } }, label);
    wrap.append(h("div", { class: "sizer" },
      btn("+ row", () => rows.push(rows[0].map(() => "0"))),
      btn("− row", () => rows.length > 1 && rows.pop()),
      btn("+ col", () => rows.forEach((r) => r.push("0"))),
      btn("− col", () => rows[0].length > 1 && rows.forEach((r) => r.pop())),
    ));
  }
  return { el: wrap, get: () => rows.map((r) => [...r]), set(v) { rows = v.map((r) => r.map(String)); draw(); } };
}

function textEditor(value, { onChange, label }) {
  const id = nextId("txt");
  const inp = h("input", { id, value, class: "text", spellcheck: "false" });
  inp.addEventListener("input", () => onChange?.(inp.value));
  return { el: h("label", { class: "field", for: id }, h("span", {}, label), inp), get: () => inp.value };
}

// ---------------------------------------------------------------------------
// The fast panel: preview now, exact when it settles

function fastPanel() {
  const quick = h("div", { class: "stage-body" });
  const exact = h("div", { class: "stage-body" });
  const qMeta = h("span", { class: "meta" });
  const eMeta = h("span", { class: "meta" });
  const el = h("div", { class: "fast" },
    h("div", { class: "stage float" }, h("div", { class: "stage-head" }, h("span", { class: "chip float" }, "preview"), qMeta), quick),
    h("div", { class: "stage exact" }, h("div", { class: "stage-head" }, h("span", { class: "chip exact" }, "exact"), eMeta), exact),
  );
  let controller = null;
  return {
    el,
    run(job, render = show) {
      controller?.abort();
      if (!(job instanceof L.Refinement)) {
        quick.replaceChildren(h("span", { class: "muted" }, "computed directly"));
        qMeta.textContent = "";
        exact.replaceChildren(render(job));
        eMeta.textContent = "";
        return;
      }
      const t0 = performance.now();
      const q = job.quick;
      const tq = performance.now() - t0;
      quick.replaceChildren(render(q.value, "float"));
      qMeta.replaceChildren(h("span", { class: `trust ${q.trust}` }, q.trust), ` ${tq.toFixed(1)} ms`);
      exact.replaceChildren(h("span", { class: "muted pending" }, "sharpening…"));
      eMeta.textContent = "";
      const mine = (controller = new AbortController());
      const t1 = performance.now();
      job.settle({
        signal: mine.signal,
        onProgress: (p) => (eMeta.textContent = `${p.step} ${p.done}/${p.of}`),
        onStage: (s) => {
          if (!s.exact) return;
          exact.replaceChildren(render(s.value, "exact"));
          eMeta.textContent = `${s.via} · ${(performance.now() - t1).toFixed(1)} ms`;
        },
      }).catch((e) => {
        if (e.name !== "AbortError") exact.replaceChildren(h("span", { class: "err" }, e.message));
      });
    },
    error(e) {
      controller?.abort();
      quick.replaceChildren(h("span", { class: "err" }, e.message));
      exact.replaceChildren();
      qMeta.textContent = eMeta.textContent = "";
    },
  };
}

// ---------------------------------------------------------------------------
// Step-by-step player

export function player(walk) {
  let i = 0;
  const list = h("ol", { class: "steps" });
  const body = h("div", { class: "step" });
  const counter = h("span", { class: "meta" });
  const prev = h("button", { type: "button", class: "ghost", onclick: () => go(i - 1) }, "← Back");
  const next = h("button", { type: "button", onclick: () => go(i + 1) }, "Next →");
  walk.steps.forEach((s, k) => list.append(h("li", {}, h("button", { type: "button", class: "steplink", onclick: () => go(k) }, s.title))));
  const go = (k) => {
    i = Math.max(0, Math.min(walk.steps.length - 1, k));
    const s = walk.steps[i];
    const parts = [h("h4", {}, s.title)];
    if (s.text) parts.push(h("p", {}, s.text));
    if (s.lines) parts.push(h("div", { class: "lines" }, s.lines.map((l) => h("div", {}, l))));
    const mats = [];
    if (s.L) mats.push(h("figure", {}, h("figcaption", {}, "L"), matrixTable(s.L)));
    if (s.matrix) mats.push(h("figure", {}, s.L ? h("figcaption", {}, "U") : null, matrixTable(s.matrix, { bar: s.bar ?? null, pivot: s.pivot ?? null, hl: s.rows ?? [], cols: s.cols ?? [] })));
    if (s.minor) mats.push(h("figure", {}, h("figcaption", {}, "minor"), matrixTable(s.minor)));
    if (mats.length) parts.push(h("div", { class: "mats" }, mats));
    if (s.table) {
      const t = h("table", { class: "data" }, h("thead", {}, h("tr", {}, s.table.head.map((c) => h("th", {}, c)))), h("tbody", {}, s.table.rows.map((r) => h("tr", {}, r.map((c) => h("td", {}, c))))));
      parts.push(h("div", { class: "mxwrap" }, t));
    }
    body.replaceChildren(...parts);
    [...list.children].forEach((li, k) => li.classList.toggle("on", k === i));
    counter.textContent = `Step ${i + 1} of ${walk.steps.length}`;
    prev.disabled = i === 0;
    next.disabled = i === walk.steps.length - 1;
  };
  const el = h("div", { class: "player" },
    h("div", { class: "player-head" }, h("strong", {}, walk.title), counter),
    h("div", { class: "player-main" }, h("nav", { class: "steplist", "aria-label": "steps" }, list), h("div", { class: "player-body" }, body, h("div", { class: "player-nav" }, prev, next))),
  );
  go(0);
  return el;
}

// ---------------------------------------------------------------------------
// A demo: inputs, a live fast panel, and a step-by-step button

const parseMatrix = (rows) => L.Matrix.from(rows.map((r) => r.map((c) => c.trim() || "0")));
const parseVector = (s) => s.split(/[,\s]+/).filter(Boolean);
export const parsePoly = (s) => s.split(/[,\s]+/).filter(Boolean).reverse(); // highest degree first in the UI

export function demo(mount, cfg) {
  const editors = {};
  const inputsEl = h("div", { class: "inputs" });
  const values = () => {
    const v = {};
    for (const inp of cfg.inputs) {
      const raw = editors[inp.id].get();
      v[inp.id] = inp.kind === "grid" ? raw.map((r) => r.map((c) => c.trim() || "0")) : inp.kind === "matrix" ? parseMatrix(raw) : inp.kind === "vector" ? parseVector(raw) : inp.kind === "poly" ? parsePoly(raw) : raw;
    }
    return v;
  };
  const fp = cfg.fast ? fastPanel() : null;
  const teachBox = h("div", { class: "teachbox" });
  let timer = null;
  const recompute = () => {
    if (!fp) return;
    clearTimeout(timer);
    timer = setTimeout(() => {
      try {
        fp.run(cfg.fast(values()), cfg.render);
      } catch (e) {
        fp.error(e);
      }
      if (teachBox.firstChild) openTeach();
    }, 60);
  };
  const openTeach = () => {
    try {
      teachBox.replaceChildren(player(cfg.teach(values())));
    } catch (e) {
      teachBox.replaceChildren(h("p", { class: "err" }, e.message));
    }
  };
  for (const inp of cfg.inputs) {
    const isGrid = inp.kind === "matrix" || inp.kind === "grid";
    const ed = isGrid
      ? matrixEditor(inp.value, { onChange: recompute, resizable: inp.resizable !== false })
      : textEditor(inp.value, { onChange: recompute, label: inp.label });
    editors[inp.id] = ed;
    inputsEl.append(isGrid ? h("div", { class: "field" }, h("span", {}, inp.label), ed.el) : ed.el);
  }
  const actions = h("div", { class: "actions" });
  if (cfg.teach) actions.append(h("button", { type: "button", onclick: () => (teachBox.firstChild ? teachBox.replaceChildren() : openTeach()) }, cfg.teachLabel ?? "Step by step"));
  if (cfg.extra) actions.append(...cfg.extra(values, recompute));
  const el = h("div", { class: "demo-card" },
    cfg.title ? h("div", { class: "demo-title" }, cfg.title) : null,
    h("div", { class: "demo-grid" }, inputsEl, fp ? fp.el : null),
    actions,
    teachBox,
  );
  mount.replaceChildren(el);
  recompute();
  return { values, recompute, editors };
}

// ---------------------------------------------------------------------------
// Do-it-yourself row reduction

export function rowReducerUI(mount, start) {
  let rr = new L.teach.RowReducer(start);
  const state = h("div", { class: "rr-state" });
  const log = h("ol", { class: "rr-log" });
  const status = h("p", { class: "rr-status" });
  const rowSel = (id) => {
    const s = h("select", { id });
    return s;
  };
  const a = rowSel(nextId("rra"));
  const b = rowSel(nextId("rrb"));
  const c = h("input", { id: nextId("rrc"), value: "1", class: "small", "aria-label": "factor" });
  const op = h("select", { id: nextId("rrop"), "aria-label": "operation" },
    h("option", { value: "add" }, "Rᵢ ← Rᵢ + c·Rⱼ"),
    h("option", { value: "scale" }, "Rᵢ ← c·Rᵢ"),
    h("option", { value: "swap" }, "Rᵢ ↔ Rⱼ"));
  const draw = () => {
    const n = rr.matrix.length;
    for (const s of [a, b]) {
      const v = s.value;
      s.replaceChildren(...Array.from({ length: n }, (_, k) => h("option", { value: k }, `R${"₀₁₂₃₄₅₆₇₈₉"[k + 1] ?? k + 1}`)));
      if (v) s.value = v;
    }
    if (!b.value || b.value === a.value) b.value = a.value === "0" ? "1" : "0";
    state.replaceChildren(matrixTable(rr.matrix));
    log.replaceChildren(...rr.history.map((x) => h("li", {}, x.label)));
    status.textContent = rr.isReduced() ? "Reduced row echelon form. Done." : rr.isEchelon() ? "Echelon form: now clear above the pivots and make them 1." : "Not in echelon form yet.";
    status.className = `rr-status ${rr.isReduced() ? "done" : ""}`;
  };
  const msg = h("p", { class: "muted small" });
  const apply = () => {
    msg.textContent = "";
    try {
      const i = Number(a.value);
      const j = Number(b.value);
      if (op.value === "swap") rr.swap(i, j);
      else if (op.value === "scale") rr.scale(i, c.value);
      else rr.add(i, j, c.value);
    } catch (e) {
      msg.textContent = e.message;
    }
    draw();
  };
  const ed = matrixEditor(start, { onChange: () => {} });
  const el = h("div", { class: "demo-card" },
    h("div", { class: "demo-title" }, "Row reduce it yourself"),
    h("div", { class: "demo-grid" },
      h("div", { class: "field" }, h("span", {}, "Starting matrix"), ed.el, h("button", { type: "button", class: "ghost", onclick: () => { try { rr = new L.teach.RowReducer(ed.get()); msg.textContent = ""; } catch (e) { msg.textContent = e.message; } draw(); } }, "Load this matrix")),
      h("div", {}, state, status),
    ),
    h("div", { class: "rr-controls" }, op, h("label", {}, "i ", a), h("label", {}, "j ", b), h("label", {}, "c ", c),
      h("button", { type: "button", onclick: apply }, "Apply"),
      h("button", { type: "button", class: "ghost", onclick: () => { rr.undo(); draw(); } }, "Undo"),
      h("button", { type: "button", class: "ghost", onclick: () => { const hi = rr.hint(); msg.textContent = hi ? `Hint: ${hi.label}. ${hi.why}` : "Nothing left to do."; } }, "Hint"),
      h("button", { type: "button", class: "ghost", onclick: () => { const hi = rr.hint(); if (hi) hi.apply(); draw(); } }, "Do the hint")),
    msg,
    log,
  );
  mount.replaceChildren(el);
  draw();
}

export { L };
